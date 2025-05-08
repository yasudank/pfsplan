import datetime
import astropy.units as u
from matplotlib.colors import LinearSegmentedColormap
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Image, Table, TableStyle, Paragraph, Spacer, PageBreak, KeepTogether
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle

# Define ANSI color codes
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def printSchedule(obsSlotList, dates_local, observer, oc, params):
    for k in range(len(dates_local)):
        xx = list()
        yy = list()
        zz = list()
        ma = list()
        print(f'Observing plan for {dates_local[k]}')
        print('Name                  Time (HST)       Airmass t_eff   RotAngs RotAnge HA     MoonSep MoonIll MoonAlt MarsSep  JupiterSep SaturnSep')
        print('--------------------- ---------------- ------- ------- ------- ------- ------ ------- ------- ------- --------- --------- ---------')
        for slot in obsSlotList:
            if slot.date == dates_local[k]:
                t = slot.target
                if t is None:
                    continue

                print(f'{t.name:21s} {(slot.obs_start+observer.utcoffset).iso[:-7]}', end=' ')

                # airmass
                if oc.airmass(slot.index, t.name) > params.airmass['limit'][t.wg]:
                    print(f'{bcolors.FAIL}{oc.airmass(slot.index, t.name):7.2f}{bcolors.ENDC}', end=' ')
                elif oc.airmass(slot.index, t.name) > params.airmass['warn']:
                    print(f'{bcolors.WARNING}{oc.airmass(slot.index, t.name):7.2f}{bcolors.ENDC}', end=' ')
                else:
                    print(f'{bcolors.OKGREEN}{oc.airmass(slot.index, t.name):7.2f}{bcolors.ENDC}', end=' ')

                # effective exposure time
                if oc.teff(slot.index, t.name) > 0.8:
                    print(f'{bcolors.OKGREEN}{oc.teff(slot.index, t.name):7.2f}{bcolors.ENDC}', end=' ')
                elif oc.teff(slot.index, t.name) > 0.6:
                    print(f'{bcolors.WARNING}{oc.teff(slot.index, t.name):7.2f}{bcolors.ENDC}', end=' ')
                else:
                    print(f'{bcolors.FAIL}{oc.teff(slot.index, t.name):7.2f}{bcolors.ENDC}', end=' ')

                rotang_warn = 150 * u.deg
                rotang_limit = 170 * u.deg
                # rotator angle
                if oc.rotang_start(slot.index, t.name) < rotang_warn and oc.rotang_start(slot.index, t.name) > -rotang_warn:
                    print(f'{bcolors.OKGREEN}{oc.rotang_start(slot.index, t.name).degree:+7.2f}{bcolors.ENDC}', end=' ')
                elif oc.rotang_start(slot.index, t.name) < rotang_limit and oc.rotang_start(slot.index, t.name) > -rotang_limit:
                    print(f'{bcolors.WARNING}{oc.rotang_start(slot.index, t.name).degree:+7.2f}{bcolors.ENDC}', end=' ')
                else:
                    print(f'{bcolors.FAIL}{oc.rotang_start(slot.index, t.name).degree:+7.2f}{bcolors.ENDC}', end=' ')

                # rotator angle
                if oc.rotang_end(slot.index, t.name) < rotang_warn and oc.rotang_end(slot.index, t.name) > -rotang_warn:
                    print(f'{bcolors.OKGREEN}{oc.rotang_end(slot.index, t.name).degree:+7.2f}{bcolors.ENDC}', end=' ')
                elif oc.rotang_end(slot.index, t.name) < rotang_limit and oc.rotang_end(slot.index, t.name) > -rotang_limit:
                    print(f'{bcolors.WARNING}{oc.rotang_end(slot.index, t.name).degree:+7.2f}{bcolors.ENDC}', end=' ')
                else:
                    print(f'{bcolors.FAIL}{oc.rotang_end(slot.index, t.name).degree:+7.2f}{bcolors.ENDC}', end=' ')


                ha_warn = (0.3 * u.hourangle).value
                ha_limit = (0.1 * u.hourangle).value
                is_north = t.coord.dec > observer.location.lat
                # hour angle
                if abs(oc.ha(slot.index, t.name).value) < ha_limit and is_north:
                    print(f'{bcolors.FAIL}{oc.ha(slot.index, t.name).value:6.2f}{bcolors.ENDC}', end=' ')
                elif abs(oc.ha(slot.index, t.name).value) < ha_warn and is_north:
                    print(f'{bcolors.WARNING}{oc.ha(slot.index, t.name).value:6.2f}{bcolors.ENDC}', end=' ')
                else:
                    print(f'{bcolors.OKGREEN}{oc.ha(slot.index, t.name).value:6.2f}{bcolors.ENDC}', end=' ')

                # moon separation
                if oc.moon_sep(slot.index, t.name) < params.moonsep['limit']:
                    print(f'{bcolors.FAIL}{oc.moon_sep(slot.index, t.name).degree:7.2f}{bcolors.ENDC}', end=' ')
                elif oc.moon_sep(slot.index, t.name) < params.moonsep['warn']:
                    print(f'{bcolors.WARNING}{oc.moon_sep(slot.index, t.name).degree:7.2f}{bcolors.ENDC}', end=' ')
                else:
                    print(f'{bcolors.OKGREEN}{oc.moon_sep(slot.index, t.name).degree:7.2f}{bcolors.ENDC}', end=' ')

                # moon illumination
                if oc.moon_ill(slot.index) < params.moonill['limit'] or oc.moon_altaz(slot.index).alt < params.moonalt['limit']:
                    print(f'{bcolors.OKGREEN}{oc.moon_ill(slot.index):7.2f}{bcolors.ENDC}', end=' ')
                elif oc.moon_ill(slot.index) < params.moonill['warn']:
                    print(f'{bcolors.WARNING}{oc.moon_ill(slot.index):7.2f}{bcolors.ENDC}', end=' ')
                else:
                    print(f'{bcolors.FAIL}{oc.moon_ill(slot.index):7.2f}{bcolors.ENDC}', end=' ')

                # Moon altitude
                if oc.moon_altaz(slot.index).alt < params.moonalt['limit']:
                    print(f'{bcolors.OKGREEN}{oc.moon_altaz(slot.index).alt.deg:7.2f}{bcolors.ENDC}', end=' ')
                elif oc.moon_altaz(slot.index).alt < params.moonalt['warn']:
                    print(f'{bcolors.WARNING}{oc.moon_altaz(slot.index).alt.deg:7.2f}{bcolors.ENDC}', end=' ')
                else:
                    print(f'{bcolors.FAIL}{oc.moon_altaz(slot.index).alt.deg:7.2f}{bcolors.ENDC}', end=' ')
                
                # Planet separation
                for planet in ["mars", "jupiter", "saturn"]:
                    if oc.planet_sep(planet, slot.index, t.name) < params.planetssep['limit']:
                        print(f'{bcolors.FAIL}{oc.planet_sep(planet, slot.index, t.name).degree:9.2f}{bcolors.ENDC}', end=' ')
                    elif oc.planet_sep(planet, slot.index, t.name) < params.planetssep['warn']:
                        print(f'{bcolors.WARNING}{oc.planet_sep(planet, slot.index, t.name).degree:9.2f}{bcolors.ENDC}', end=' ')
                    else:
                        print(f'{bcolors.OKGREEN}{oc.planet_sep(planet, slot.index, t.name).degree:9.2f}{bcolors.ENDC}', end=' ')
                print()
        print()

def createTableContents(obsSlotList, observer, oc, params):
    # Define the color sequence
    custom_colors = ['green', 'yellow', 'orange', 'red']

    # Create a custom colormap
    custom_cmap = LinearSegmentedColormap.from_list('danger_scale', custom_colors)

    texts = list()
    tables = list()
    tablestyles = list()
    for k in range(len(obsSlotList.dates)):
        texts.append(f'Observing plan for {obsSlotList.dates[k]}')
        data = list()
        l_tablestyle = list()
        #data.append(['Name', 'Time (HST)', 'Airmass', 't_eff', 'RotAng(Start)', 'RotAng(End)', 'MoonSep', 'MoonIll', 'MoonAlt', 'MarsSep', 'JupiterSep', 'SaturnSep'])
        data.append(['', '', '', '', 'Rotator Angle', '', '', '', ''])
        data.append(['Name', 'Time (HST)', 'Airmass', 't_eff', 'Start', 'End', 'HA', 'MoonSep', 'MoonIll', 'MoonAlt'])
        jj = 1
        for slot in obsSlotList.slots:
            if slot.start.iso.split(' ')[0] == obsSlotList.dates_utc[k]:
                i = slot.index
                if slot.target is None:
                    continue
                tname = slot.target.name
                data.append([tname, (slot.start+5*u.minute+observer.utcoffset).iso[:-7],
                             f'{oc.airmass(i, tname):7.2f}',
                             f'{oc.teff(i, tname):7.2f}',
                             f'{oc.rotang_start(i, tname).degree:+7.2f}',
                             f'{oc.rotang_end(i, tname).degree:+7.2f}',
                             f'{oc.ha(i, tname).value:6.2f}',
                             f'{oc.moon_sep(i, tname).degree:7.2f}', 
                             f'{oc.moon_ill(i):7.2f}', 
                             f'{oc.moon_altaz(i).alt.deg:7.2f}', 
                            # f'{oc.planet_sep('mars', i, tname).degree:7.2f}', 
                            # f'{oc.planet_sep('jupiter', i, tname).degree:7.2f}', 
                            # f'{oc.planet_sep('saturn', i, tname).degree:7.2f}'
                            ])

                if oc.airmass(i, tname) > params.airmass['limit'][slot.target.wg]:
                    l_tablestyle.append(('BACKGROUND', (2, jj+1), (2, jj+1), colors.lightcoral))
                elif oc.airmass(i, tname) > params.airmass['warn']:
                    l_tablestyle.append(('BACKGROUND', (2, jj+1), (2, jj+1), colors.lemonchiffon))               

                if oc.teff(i, tname) < 0.6:
                    l_tablestyle.append(('BACKGROUND', (3, jj+1), (3, jj+1), colors.lightcoral))
                elif oc.teff(i, tname) < 0.8:
                    l_tablestyle.append(('BACKGROUND', (3, jj+1), (3, jj+1), colors.lemonchiffon))

                rotang_warn = 150 * u.deg
                rotang_limit = 170 * u.deg
                if oc.rotang_start(i, tname) > rotang_limit or oc.rotang_start(i, tname) < -rotang_limit:
                    l_tablestyle.append(('BACKGROUND', (4, jj+1), (4, jj+1), colors.lightcoral))
                elif oc.rotang_start(i, tname) > rotang_warn or oc.rotang_start(i, tname) < -rotang_warn:
                    l_tablestyle.append(('BACKGROUND', (4, jj+1), (4, jj+1), colors.lemonchiffon))

                if oc.rotang_end(i, tname) > rotang_limit or oc.rotang_end(i, tname) < -rotang_limit:
                    l_tablestyle.append(('BACKGROUND', (5, jj+1), (5, jj+1), colors.lightcoral))
                elif oc.rotang_end(i, tname) > rotang_warn or oc.rotang_end(i, tname) < -rotang_warn:
                    l_tablestyle.append(('BACKGROUND', (5, jj+1), (5, jj+1), colors.lemonchiffon))

                ha_warn = (0.3 * u.hourangle).value
                ha_limit = (0.1 * u.hourangle).value
                is_north = slot.target.coord.dec > observer.location.lat
                if abs(oc.ha(i, tname).value) < ha_limit and is_north:
                    l_tablestyle.append(('BACKGROUND', (6, jj+1), (6, jj+1), colors.lightcoral))
                elif abs(oc.ha(i, tname).value) < ha_warn and is_north:
                    l_tablestyle.append(('BACKGROUND', (6, jj+1), (6, jj+1), colors.lemonchiffon))

                if oc.moon_sep(i, tname) < params.moonsep['limit']:
                    l_tablestyle.append(('BACKGROUND', (7, jj+1), (7, jj+1), colors.lightcoral))
                elif oc.moon_sep(i, tname) < params.moonsep['warn']:
                    l_tablestyle.append(('BACKGROUND', (7, jj+1), (7, jj+1), colors.lemonchiffon))

                if oc.moon_altaz(i).alt > params.moonalt['limit']:
                    if oc.moon_ill(i) > params.moonill['limit']:
                        l_tablestyle.append(('BACKGROUND', (8, jj+1), (8, jj+1), colors.lightcoral))
                    elif oc.moon_ill(i) > params.moonill['warn']:
                        l_tablestyle.append(('BACKGROUND', (8, jj+1), (8, jj+1), colors.lemonchiffon))

                if oc.moon_altaz(i).alt > params.moonalt['limit']:
                    l_tablestyle.append(('BACKGROUND', (9, jj+1), (9, jj+1), colors.lightcoral))
                elif oc.moon_altaz(i).alt > params.moonalt['warn']:
                    l_tablestyle.append(('BACKGROUND', (8, jj+1), (8, jj+1), colors.lemonchiffon))
                """
                if mars_sep[(i, tname)] < planetssep_limit:
                    l_tablestyle.append(('BACKGROUND', (9, jj+1), (9, jj+1), colors.lightcoral))
                elif mars_sep[(i, tname)] < planetssep_warn:
                    l_tablestyle.append(('BACKGROUND', (9, jj+1), (9, jj+1), colors.lemonchiffon))
                if jupiter_sep[(i, tname)] < planetssep_limit:
                    l_tablestyle.append(('BACKGROUND', (10, jj+1), (10, jj+1), colors.lightcoral))
                elif jupiter_sep[(i, tname)] < planetssep_warn:
                    l_tablestyle.append(('BACKGROUND', (10, jj+1), (10, jj+1), colors.lemonchiffon))
                if saturn_sep[(i, tname)] < planetssep_limit:
                    l_tablestyle.append(('BACKGROUND', (11, jj+1), (11, jj+1), colors.lightcoral))
                elif saturn_sep[(i, tname)] < planetssep_warn:
                    l_tablestyle.append(('BACKGROUND', (11, jj+1), (11, jj+1), colors.lemonchiffon))
                """

                jj += 1
            
        tables.append(data)
        tablestyles.append(l_tablestyle)

    return texts, tables, tablestyles

# Function to write text and table data to a multi-page PDF
def write_text_and_table_to_pdf(texts, tables, tablestyles, params):

    # Get the current date
    current_date = datetime.datetime.now()

    # Format the version string with current date as YYYYMMDD
    version_date = current_date.strftime("%Y%m%d")

    # Extract the base filename (without .pdf)
    base_filename = params.fname_report.replace('.pdf', '')

    # Create the new filename with version date
    fname_report_versioned = f"{base_filename}.v{version_date}.pdf"

    pdf = SimpleDocTemplate(fname_report_versioned, pagesize=A4)
    elements = []

    # Define a custom style for the text
    styles = getSampleStyleSheet()
    custom_style = ParagraphStyle(
        name='CustomStyle',
        parent=styles['Normal'],
        fontSize=16,  # Set font size
        alignment=0,  # Left alignment
        spaceAfter=12  # Space after paragraph
    )

    # Add text with the custom style
    text = f'Observing plan for {params.fname_obsdate}'
    paragraph = Paragraph(text, custom_style)
    elements.append(paragraph)

    img1 = Image('elevation.png', width=7*inch, height=4.5*inch)
    elements.append(img1)
    # Add some space between images
    elements.append(Spacer(1, 12))
    img2 = Image('observed_counts.png', width=7*inch, height=4.5*inch)
    elements.append(img2)
    # Add some space between images
    elements.append(Spacer(1, 12))

    img3 = Image('rot_angle.png', width=7*inch, height=4.5*inch)
    elements.append(img3)
    # Add some space between images
    elements.append(Spacer(1, 12))
    img4 = Image('hour_angle.png', width=7*inch, height=4.5*inch)
    elements.append(img4)
    # Add some space between images
    elements.append(Spacer(1, 12))

    # Add a page break after each page except the last one
    elements.append(PageBreak())

    for text, data, ts in zip(texts, tables, tablestyles):
        # Create a list to hold the elements that should be kept together
        keep_together_elements = []

        # Add text with the custom style
        paragraph = Paragraph(text, custom_style)
        keep_together_elements.append(paragraph)

        # Add some space between text and table
        keep_together_elements.append(Spacer(1, 12))

        # Create the table
        table = Table(data)

        # Add style to the table
        style = TableStyle([
            ('SPAN', (4, 0), (5, 0)),
            ('BACKGROUND', (0, 0), (-1, 1), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 1), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 1), 9),
            #('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            #('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ] + ts)
        table.setStyle(style)

        # Add the table to the elements list
        keep_together_elements.append(table)

        # Add some space between text and table
        keep_together_elements.append(Spacer(1, 12))

        # Add a page break after each page except the last one
        #elements.append(PageBreak())

        # Add the keep_together_elements to the main elements list using KeepTogether
        elements.append(KeepTogether(keep_together_elements))

    # Build the PDF
    pdf.build(elements)

def printSchedule_PDF(obsSlotList, observer, oc, params):
    texts, tables, tablestyles = createTableContents(obsSlotList, observer, oc, params)

    write_text_and_table_to_pdf(texts, tables, tablestyles, params)