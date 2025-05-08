from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from ObservingConditions import slewTime
import logging
import os # for os.path.exists

# Configure logging
logger = logging.getLogger(__name__)

class ObsSlot:

    def __init__(self, index, start, end, mid, obs_start, obs_end, date, date_utc, used=False):
        self.index = index
        self.start = start
        self.end = end
        self.mid = mid
        self.obs_start = obs_start
        self.obs_end = obs_end
        self.date = date
        self.date_utc = date_utc
        self.used = used
        self.target = None

    def __repr__(self):
        return f"ObsSlot({self.index}, {self.used}, {self.date}, {self.mid}, {self.obs_start}, {self.obs_end})"
        
class ObsSlotList:

    def __init__(self):
        self.slots = []
        self.used_slots = []
        self.free_slots = []
        self.dates = []
        self.dates_utc = []
        self._slot_index_map = {} # add this

    def add_slot(self, slot):
        self.slots.append(slot)
        if not slot.date in self.dates:
            self.dates.append(slot.date)
        if not slot.date_utc in self.dates_utc:
            self.dates_utc.append(slot.date_utc)
        self._slot_index_map[slot.index] = len(self.slots) -1 # add this        
        if slot.used:
            self.used_slots.append(slot)
        else:
            self.free_slots.append(slot)

    def get_used_slots(self):
        return self.used_slots

    def get_free_slots(self):
        return self.free_slots
    
    def get_all_slots(self):
        return self.slots
    
    def get_slot_by_index(self, index):
        for slot in self.slots:
            if slot.index == index:
                return slot
        return None
    
    def get_slots_by_field(self, field):
        return [slot for slot in self.used_slots if slot.target.name.startswith(field)]
    
    def get_slots_by_date(self, date):
        return [slot for slot in self.slots if slot.date == date]
    
    def get_available_slots(self, nslot):
        _slots = []
        for date in self.dates:
            for slot in self.free_slots:
                if slot.date == date:
                    _slots.append(slot)
            if len(_slots) >= nslot:
                break
        return _slots

    def get_index(self, index):
        return self._slot_index_map.get(index) # change this        
    
    @property
    def num_slots(self):
        return len(self.slots)

    def __getitem__(self, index):
        """
        リストのようにインデックスアクセスを可能にします。

        Args:
            index (int): アクセスするスロットのインデックス。

        Returns:
            ObsSlot: 指定されたインデックスのスロット。

        Raises:
            IndexError: インデックスが範囲外の場合。
        """
        if isinstance(index, int):
            if 0 <= index < len(self.slots):
                return self.slots[index]
            else:
                raise IndexError("ObsSlotList index out of range")
        else:
            raise TypeError("Index must be an integer")
        
    def updateUsed(self, index):
        """
        スロットの使用状況を更新します。

        Args:
            index (int): 更新するスロットのインデックス。
        """
        _index = self.get_index(index)
        if _index is not None:
            self.slots[_index].used = True
            self.used_slots.append(self.slots[_index])
            self.free_slots.remove(self.slots[_index])
        else:
            raise ValueError("Invalid index for ObsSlotList")

    def updateSchedule(self, o, obs_slots, targets, targetList):
        for slot in obs_slots:
            for target in targets:
                if o[(slot.index, target.name)].varValue > 0.5:
                    slot.target = target
                    self.updateUsed(slot.index)
                    targetList.update_observed(target.name, 1)
                    break

    def updateTimeBySlew(self, oc, params):
        slew_overhead = 0.0 * u.minute
        for i in range(self.num_slots-1):
            cur_slot = self.slots[i]
            tgt_slot = self.slots[i+1]
            if cur_slot.date != tgt_slot.date:
                # logger.info("") # Removed empty print
                slew_overhead = 0.0 * u.minute
                continue
            cur_target = cur_slot.target
            tgt_target = tgt_slot.target
            if cur_target == None or tgt_target == None:
                logger.warning(f"Slot {cur_slot.index} or {tgt_slot.index} is empty during slew time calculation")
                continue
            cur_altaz = oc.altaz(cur_slot.index, cur_target.name)
            tgt_altaz = oc.altaz(tgt_slot.index, tgt_target.name)
            cur_rotang = oc.rotang_end(cur_slot.index, cur_target.name)
            tgt_rotang = oc.rotang_start(tgt_slot.index, tgt_target.name)
            slew_time = slewTime(cur_altaz, cur_rotang, tgt_altaz, tgt_rotang, params)
            logger.info(f"Slew Time: {cur_slot.index:3d} -> {tgt_slot.index:3d} ({cur_target.name:21s} -> {tgt_target.name:21s}) = {slew_time.to(u.minute).value:4.1f} min. Cumulative slew: {slew_overhead.to(u.minute).value:4.1f} min")
            slew_overhead += slew_time

            tgt_slot.start += slew_overhead
            tgt_slot.end += slew_overhead
            tgt_slot.mid += slew_overhead
            tgt_slot.obs_start += slew_overhead
            tgt_slot.obs_end += slew_overhead

    def reset(self):
        """
        スロットの使用状況をリセットします。
        """
        for slot in self.slots:
            slot.used = False
            slot.target = None
        self.used_slots = []
        self.free_slots = self.slots[:]

    def __iter__(self):
        return iter(self.slots)
    
    def __len__(self):
        return len(self.slots)
    
    def __repr__(self):
        return f"ObsSlotList({len(self.slots)} slots)"
    
    def __str__(self):
        return f"ObsSlotList with {len(self.slots)} slots"
    
    def __contains__(self, slot):
        return slot in self.slots

class ObsDate:

    # Return UTC time from the local time
    # Hours after 24:00 will be converted to the next day
    def get_time(self, date, time, utcoffset):
        hh = int(time.split(':')[0])
        if hh >= 24:
            time = f'{(hh-24):02d}'+':'+time.split(':')[1]+':00'
            date = (Time(date, format='iso') + 1 * u.day).iso.split(' ')[0]
        return Time(date+' '+time) - utcoffset

    def __init__(self, fname_obsdate, fname_obsdate_finish=None, observer=None, params=None):

        obsdate_table = Table.read(fname_obsdate, format='ascii')

        self.dates = [row['date'] for row in obsdate_table]

        self.params = params

        if fname_obsdate_finish and os.path.exists(fname_obsdate_finish):
            self.obsdate_finish_table = Table.read(fname_obsdate_finish, format='ascii')
            self.dates_finish = [row['date'] for row in self.obsdate_finish_table]
        else:
            self.dates_finish = []

        self.dates_local = []
        self.dates_utc = []

        midpt = params.t_overhead + 0.5 * (params.w_timeslot - params.t_overhead)

        self.obsSlotList = ObsSlotList()
        slot_index = 0

        for i, (date, start, end) in enumerate(obsdate_table):
            # Set time to noon at Hawaii to calculate the "next" sunset and sunrise
            time = Time(date+' 12:00:00') - observer.utcoffset
            
            if start == 'sun_set':
                start_time = observer.sun_set_time(time, which='next', horizon=params.angle_twilight)
            else:
                start_time = self.get_time(date, start, observer.utcoffset)
            if end == 'sun_rise':
                end_time = sun_rise = observer.sun_rise_time(time, which='next', horizon=params.angle_twilight)
            else:
                end_time = self.get_time(date, end, observer.utcoffset)
            # Print the time range for the observation in HST
            logger.info(f"ObsDate processing: {date} HST Start: {(start_time+observer.utcoffset).iso} HST End: {(end_time+observer.utcoffset).iso}")

            if not date in self.dates_finish:
                if not date in self.dates_local:
                    self.dates_local.append(date)
                date_utc = start_time.strftime('%Y-%m-%d')
                if not date_utc in self.dates_utc:
                    self.dates_utc.append(date_utc)

            used = False
            if date in self.dates_finish:
                used = True

            current = start_time
            while current + midpt < end_time:
                slot_index += 1
                obsSlot = ObsSlot(index=slot_index,
                                  start=current,
                                  end=current + params.w_timeslot,
                                  mid=current + midpt,    
                                  obs_start=current + params.t_overhead,
                                  obs_end=current + params.w_timeslot,
                                  date=date,
                                  date_utc=date_utc,
                                  used=used)
                self.obsSlotList.add_slot(obsSlot)
                current += params.w_timeslot

        self.dates_local.sort()
        self.dates_utc.sort()

    @property
    def nexp_max(self):
        return {w: int((self.params.frac[w] + 0.5 * self.params.frac_margin) * self.obsSlotList.num_slots) for w in self.params.frac.keys()}

    @property
    def nexp_min(self):
        return {w: int((self.params.frac[w] - 0.5 * self.params.frac_margin) * self.obsSlotList.num_slots) for w in self.params.frac.keys()}
