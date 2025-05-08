import os
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
import re

# Function to find the longest common prefix among a list of strings.
def longest_common_prefix(strings):
    """
    Finds the longest common prefix among a list of strings.

    Args:
        strings: A list of strings.

    Returns:
        The longest common prefix string.
    """
    if not strings:
        return ""
    # Start with the first string as the candidate for the prefix.
    prefix = strings[0]
    # Iterate over the rest of the strings.
    for s in strings[1:]:
        # Reduce the prefix until it matches the start of s.
        while not s.startswith(prefix):
            prefix = prefix[:-1]
            if not prefix:
                return ""
    return prefix

def natural_sort_key(s):
    """
    Generate a key for natural sorting where numbers within strings are sorted numerically
    rather than lexicographically.
    """
    return [int(text) if text.isdigit() else text.lower() 
            for text in re.split(r'(\d+)', s)]

# Define the directory where the target data files are located.
ref_dir = '/home/yasuda/spt_ssp_observation/runs/2025-05/targets'

# Define the list of working groups (WG).
wg_list = ['GA', 'GE', 'CO']
wg_list_exist = []

# Initialize an empty dictionary to store data for each working group.
d = dict()

# Loop through each working group.
for wg in wg_list:

    if os.path.exists(f"{ref_dir}/{wg}/ppcList.ecsv") == False:
        print(f"File not found: {ref_dir}/{wg}/ppcList.ecsv")
        continue

    # Append the existing working group to wg_list_exist.
    wg_list_exist.append(wg)

    # Read the target data file for the current working group.
    data = ascii.read(f"{ref_dir}/{wg}/ppcList.ecsv")

    # Extract RA and Dec from the data.
    ra = data['ppc_ra']
    dec = data['ppc_dec']

    # Create SkyCoord objects for each target.
    c = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))

    # Convert RA and Dec to string format (HH:MM:SS and +/-DD:MM:SS).
    ra_str = c.ra.to_string(unit=u.hourangle, sep=':', precision=2, pad=True)
    dec_str = c.dec.to_string(sep=':', precision=1, alwayssign=True, pad=True)

    # Add the formatted RA and Dec strings to the data table.
    data['ra'] = ra_str
    data['dec'] = dec_str

    # Add the working group name to the data table.
    data['wg'] = wg

    # Select the relevant columns for each working group.
    d[wg] = data['ppc_code', 'ppc_ra', 'ppc_dec', 'ppc_pa', 'ppc_resolution', 'ppc_priority', 'ppc_exptime', 'ppc_nframes', 'ra', 'dec', 'wg']
    #if wg != 'GE':
    #    d[wg] = data['ppc_code', 'ppc_ra', 'ppc_dec', 'ppc_pa', 'ppc_resolution', 'ppc_priority', 'ppc_exptime', 'ppc_nframes', 'ra', 'dec', 'wg']
    #else:
    #    # Special handling for 'GE' working group: rename 'ppc_nframe' to 'ppc_nframes'.
    #    d[wg] = data['ppc_code', 'ppc_ra', 'ppc_dec', 'ppc_pa', 'ppc_resolution', 'ppc_priority', 'ppc_exptime', 'ppc_nframe', 'ra', 'dec', 'wg']
    #    d[wg].rename_column('ppc_nframe', 'ppc_nframes')

wg_list = wg_list_exist

# Vertically stack the data tables for all working groups.
raw_targets = vstack([d[wg] for wg in wg_list])

# Create a list of unique combinations of (ra, dec, pa, priority).
uniq_rdpp = list(set([(raw_targets[i]['ra'], raw_targets[i]['dec'], raw_targets[i]['ppc_pa'], raw_targets[i]['ppc_priority']) for i in range(len(raw_targets))]))

# Create an empty table to store the final target information.
target_table = Table(names=('wg', 'name', 'ra', 'dec', 'pa', 'nexp', 'priority'),
                     dtype=(str, str, str, str, float, int, int))

# Loop through each unique combination of (ra, dec, pa, priority).
for ra, dec, pa, priority in uniq_rdpp:
    # Select the rows in raw_targets that match the current combination.
    sub = raw_targets[(raw_targets['ra'] == ra) & (raw_targets['dec'] == dec) & (raw_targets['ppc_pa'] == pa) & (raw_targets['ppc_priority'] == priority)]
    
    # Find the longest common prefix of the 'ppc_code' for the selected rows.
    ppc_code = longest_common_prefix([sub[i]['ppc_code'] for i in range(len(sub))])

    # Special handling for 'SSP_GA' and 'EN1' prefixes.
    if ppc_code.startswith('SSP_GA'):
        ppc_code = ppc_code[:ppc_code.rfind('V')]
    elif ppc_code.startswith('EN1'):
        ppc_code = ppc_code[:ppc_code.rfind('_')]

    # Calculate the total exposure time and total number of frames.
    tot_exptime = sub['ppc_exptime'].sum()
    tot_nframes = sub['ppc_nframes'].sum()

    # Get the working group name (should be the same for all rows in 'sub').
    wg = list(set(sub['wg']))[0]

    # Calculate the number of exposures (assuming 900 seconds per exposure).
    nexp = int(tot_exptime/900)

    # Add a row to the target_table with the processed information.
    target_table.add_row([wg, ppc_code, ra, dec, pa, nexp, priority])

# Sort the target table by working group and target name.
target_table.sort(['wg', 'name'])

# First sort by working group
target_table.sort('wg')

# Then sort by name within each working group
wgs = set(target_table['wg'])
sorted_table = Table(names=target_table.colnames, dtype=target_table.dtype)

for group in sorted(wgs):
    group_rows = target_table[target_table['wg'] == group]

    if group == "CO":
        # Use natural sort for CO working group
        sorted_indices = sorted(range(len(group_rows)), key=lambda i: natural_sort_key(group_rows['name'][i]))
    else:
        # Use default sort for other working groups
        sorted_indices = sorted(range(len(group_rows)), key=lambda i: group_rows['name'][i])

    for i in sorted_indices:
        sorted_table.add_row(group_rows[i])

target_table = sorted_table

# Print the final target table.
print(target_table)

# Write the final target table to an ECSV file.
target_table.write('target_table_output_202505.ecsv', format='ascii.ecsv', overwrite=True)
