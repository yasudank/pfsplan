from astropy.time import Time
import astropy.units as u
from astropy.coordinates import Angle, get_body, AltAz
from astroplan import moon_illumination
from Moon import MoonBrightnessModel as MBM
import math
import numpy as np
import pickle
import hashlib
import json
import base64
from tqdm import tqdm

# Configure logging
import logging
logger = logging.getLogger(__name__)

def generate_unique_id_base64(obsSlotList, targetList, observer, num_bytes=10):
    """
    SHA256 を使用し、Base64 エンコードで短縮します。

    obsSlotList, targetList, observer からユニークな ID を生成します。

    Args:
        obsSlotList: ObsSlotList オブジェクト
        targetList: TargetList オブジェクト
        observer: Observer オブジェクト

    Returns:
        str: ユニークな ID
    """

    # obsSlotList の状態を文字列化
    obs_slot_state = []
    for slot in obsSlotList.get_all_slots():
        obs_slot_state.append({
            "index": slot.index,
            "start": slot.start.iso if hasattr(slot, 'start') else None
         })
    obs_slot_str = json.dumps(obs_slot_state, sort_keys=True)

    # targetList の状態を文字列化
    target_state = []
    for target in targetList.get_all_targets():
        target_state.append({
            "name": target.name,
            "coord": f"{target.coord.ra.deg},{target.coord.dec.deg}"
        })
    target_str = json.dumps(target_state, sort_keys=True)

    # observer の状態を文字列化
    observer_str = json.dumps({
        "location": f"{observer.location.lat.deg},{observer.location.lon.deg},{observer.location.height.value}",
        "utcoffset": observer.utcoffset.value if hasattr(observer, 'utcoffset') else None
    }, sort_keys=True)

    # 全ての文字列を連結
    combined_str = obs_slot_str + target_str + observer_str

    # SHA256 でハッシュ化
    hash_object = hashlib.sha256(combined_str.encode())
    digest = hash_object.digest()

    # ハッシュ値の先頭 num_bytes バイトだけを使用
    shortened_digest = digest[:num_bytes]

    # Base64 エンコード
    unique_id = base64.urlsafe_b64encode(shortened_digest).decode('utf-8').rstrip('=')

    return unique_id

class ObservingConditions:

    def __init__(self, obsSlotList, targetList, observer, params):

        self.obsSlotList = obsSlotList
        self.targetList = targetList
        self.observer = observer
        self.params = params

        # Pre-compute and store target indices
        self._target_indices = {target.name: i for i, target in enumerate(targetList.get_all_targets())}

        # Pre-compute and store slot indices
        self._slot_indices = {slot.index: i for i, slot in enumerate(obsSlotList.get_all_slots())}

        mid_times = Time([slot.mid for slot in obsSlotList.get_all_slots()])
        start_times = Time([slot.obs_start for slot in obsSlotList.get_all_slots()])
        end_times = Time([slot.obs_end for slot in obsSlotList.get_all_slots()])
        
        target_coords = [target.coord for target in targetList.get_all_targets()]
        target_pa     = [target.pa for target in targetList.get_all_targets()]
        num_targets = len(target_coords)
        
        logger.info("Calculating airmass and hour angle")
        self._altaz = observer.altaz(mid_times, target_coords, 
                                     grid_times_targets=True)
        
        self._airmass = self._altaz.secz
        # Set airmass to a larget value for targets below 0.573 deg (=> airmass = 100)
        self._airmass[self._altaz.alt < 0.573 * u.deg] = 100.0 

        hour_angles = []
        for i, t in enumerate(target_coords):
            hour_angle = mid_times.sidereal_time('mean', longitude=observer.longitude) - t.ra
            #hour_angle = hour_angle.wrap_at(180 * u.deg)
            hour_angles.append(hour_angle)
        self._hour_angles = Angle(hour_angles)

        logger.info("Calculating rotator angle")
        parallactic_angle_at_start = observer.parallactic_angle(start_times, target_coords,
                                                                grid_times_targets=True)
        parallactic_angle_at_end   = observer.parallactic_angle(end_times,   target_coords,
                                                                grid_times_targets=True)
        self._rot_angle_at_start = Angle([parallactic_angle_at_start[i] + target_pa[i] for i in range(num_targets)]).wrap_at(180 * u.deg)
        self._rot_angle_at_end   = Angle([parallactic_angle_at_end[i]   + target_pa[i] for i in range(num_targets)]).wrap_at(180 * u.deg)
        
        logger.info("Calculating the separation from the Moon and Planets")
        moon = get_body('moon', mid_times, observer.location)
        self._moon_sep = [moon.separation(target_coords[i], origin_mismatch="ignore") for i in range(num_targets)]

        self._moon_ill = moon_illumination(mid_times)

        sun = get_body('sun', mid_times, observer.location)
        self._moon_phase = moon.separation(sun, origin_mismatch="ignore")

        frame = AltAz(obstime=mid_times, location=observer.location)
        self._moon_altaz = moon.transform_to(frame)

        self._planet_seps = {}
        for planet in ["mars", "jupiter", "saturn"]:
            planet_pos = get_body(planet, mid_times, observer.location)
            self._planet_seps[planet] = [planet_pos.separation(target_coords[i], origin_mismatch="ignore") for i in range(num_targets)]

        logger.info("Calculating effective exposure time")
        mbm = MBM()

        self._teff = []
        for i in range(num_targets):
            # Calculate the minimum zenith distance for the tareget
            zmin = abs(target_coords[i].dec - observer.location.lat)
            
            # Normalize the effective exposure time at the minimum zenith distance
            airmass0 = 1.0 / math.cos(zmin.to(u.rad).value)
            teff0 = 1.0 / (airmass0 * 10**(0.8*mbm.k['r']*(airmass0-1.0)))
            dmu = mbm.deltaMag("r",
                               self._moon_phase.deg,
                               90.-self._moon_altaz.alt.deg,
                               90.-self._altaz[i].alt.deg,
                               self._moon_sep[i].deg)
            dmu[self._moon_altaz.alt < 0] = 0.0
            self._teff.append((1.0 / (10**(-0.4*dmu) * self._airmass[i] * 10**(0.8*mbm.k['r']*(self._airmass[i]-1.0))) / teff0).value)

        self._teff = np.array(self._teff)

        uniq_id = generate_unique_id_base64(obsSlotList, targetList, observer)
        logger.info(f"Unique ID for ObservingConditions cache: {uniq_id}")

        logger.info("Calculating slew time")
        self.calc_slewTime(uniq_id)

    def airmass(self, islot, tname):
        # Use the pre-computed indices for fast lookup
        return self._airmass[self._target_indices[tname]][self._slot_indices[islot]]

    def altaz(self, islot, tname):
        return self._altaz[self._target_indices[tname]][self._slot_indices[islot]]

    def ha(self, islot, tname):
        return self._hour_angles[self._target_indices[tname]][self._slot_indices[islot]]

    def rotang_start(self, islot, tname):
        return self._rot_angle_at_start[self._target_indices[tname]][self._slot_indices[islot]]
    
    def rotang_end(self, islot, tname):
        return self._rot_angle_at_end[self._target_indices[tname]][self._slot_indices[islot]]
    
    def teff(self, islot, tname):
        if tname == 'dummy':
            return 0.0
        else:
            return self._teff[self._target_indices[tname]][self._slot_indices[islot]]
    
    def moon_sep(self, islot, tname):
        return self._moon_sep[self._target_indices[tname]][self._slot_indices[islot]]

    def moon_ill(self, islot):
        return self._moon_ill[self._slot_indices[islot]]
    
    def moon_altaz(self, islot):
        return self._moon_altaz[self._slot_indices[islot]]
    
    def moon_phase(self, islot):
        return self._moon_phase[self._slot_indices[islot]]

    def planet_sep(self, name, islot, tname):
        if name not in self._planet_seps:
            raise ValueError(f"Planet {name} not found.")
        return self._planet_seps[name][self._target_indices[tname]][self._slot_indices[islot]]

    def calc_slewTime(self, uniq_id):
        # Load the slew time from a pickle file if available
        try:
            with open(f'slew_time_{uniq_id}.pkl', 'rb') as f:
                self._slewTime = pickle.load(f)
                logger.info(f"Slew time loaded from file: slew_time_{uniq_id}.pkl")
            return
        except FileNotFoundError:
            logger.info(f"Slew time file slew_time_{uniq_id}.pkl not found. Calculating slew time.")

            pbar = tqdm(total=self.targetList.num_targets * self.targetList.num_targets * (self.obsSlotList.num_slots-1), desc="Calculating slew time")

            self._slewTime = create_3d_array(self.targetList.num_targets, self.targetList.num_targets, self.obsSlotList.num_slots-1)
        
            for j in range(self.obsSlotList.num_slots-1):
                if self.obsSlotList[j].date != self.obsSlotList[j+1].date:
                    continue
                for i1 in range(self.targetList.num_targets):
                    for i2 in range(self.targetList.num_targets):
                        cur_altaz = self._altaz[i1][j]
                        cur_rotang = self._rot_angle_at_end[i1][j]
                        tgt_altaz = self._altaz[i2][j+1]
                        tgt_rotang = self._rot_angle_at_start[i2][j+1]

                        self._slewTime[i1][i2][j] = slewTime(cur_altaz, cur_rotang, tgt_altaz, tgt_rotang, self.params)

                        pbar.update(1)
            pbar.close()

            # Save the slew time to a pickle file for future use
            try:
                with open(f'slew_time_{uniq_id}.pkl', 'wb') as f:
                    pickle.dump(self._slewTime, f)
            except Exception as e:
                logger.error(f"Error saving slew time to slew_time_{uniq_id}.pkl: {e}")

def slewTime(cur_altaz, cur_rotang, tgt_altaz, tgt_rotang, params):
    """
    Calculate the slew time based on the current and target altaz and rotator angles.
    """

    rate_az = params.slew_speed_az
    rate_el = params.slew_speed_el
    rate_rot = params.inst_rot_speed

    # Calculate the difference of the azimuth
    # If the difference is larger than 180 degrees, the telescope should rotate in the opposite direction
    az_diff = (tgt_altaz.az - cur_altaz.az).wrap_at(180 * u.deg)
    el_diff = tgt_altaz.alt - cur_altaz.alt
    rot_diff = tgt_rotang - cur_rotang

    # Calculate the slew time using a simple model
    slew_time = max(abs(az_diff) / rate_az,
                    abs(el_diff) / rate_el,
                    abs(rot_diff) / rate_rot)

    return slew_time
    
def create_3d_array(depth, height, width, initial_value=None):
    """
    指定された次元数の3次元配列を作成します。

    Args:
        depth (int): 奥行きの次元数。
        height (int): 高さの次元数。
        width (int): 幅の次元数。
        initial_value: 配列の初期値。デフォルトは None。

    Returns:
        list: 3次元配列。
    """
    array_3d = []
    for _ in range(depth):
        array_2d = []
        for _ in range(height):
            array_1d = [initial_value] * width
            array_2d.append(array_1d)
        array_3d.append(array_2d)
    return array_3d
