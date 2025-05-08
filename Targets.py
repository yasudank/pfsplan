from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import re

def natural_sort_key(s):
    """
    Generate a key for natural sorting where numbers within strings are sorted numerically
    rather than lexicographically.
    """
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]

class Target:
    def __init__(self, wg, name, coord, pa, nexp, priority, observed=0):
        self.wg = wg
        self.name = name
        self.coord = coord
        self.pa = pa
        self.nexp = nexp
        self.priority = priority
        self.observed = observed

    def __repr__(self):
        return f"Target({self.wg}, {self.name}, {self.coord}, {self.pa}, {self.nexp}, {self.priority}, {self.observed})"
    
class TargetList:
    def __init__(self):
        self.targets = []
        self.completed_targets = []
        self.observing_targets = []
        self._target_index_map = {} # add this

    def add_target(self, target):
        self.targets.append(target)
        self._target_index_map[target.name] = len(self.targets) -1 # add this
        if target.observed >= target.nexp:
            self.completed_targets.append(target)
        else:
            self.observing_targets.append(target)

    def sort_targets(self):
        self.targets.sort(key=lambda x: x.name)
        self.observing_targets.sort(key=lambda x: x.name)
        self.completed_targets.sort(key=lambda x: x.name)
        
    def get_completed_targets(self):
        return self.completed_targets

    def get_observing_targets(self):
        return self.observing_targets
    
    def get_all_targets(self):
        return self.targets
    
    def get_observing_targets_by_priority(self, priority):
        return [target for target in self.observing_targets if target.priority <= priority]
    
    def add_observed(self, target_name, nexp):
        for target in self.targets:
            if target.name == target_name:
                target.observed += nexp
                if target.observed >= target.nexp:
                    self.completed_targets.append(target)
                    self.observing_targets.remove(target)
                break

    def get_index(self, target_name):
        return self._target_index_map.get(target_name) # change this
    
    @property
    def num_targets(self):
        return len(self.targets)

    @property
    def wg_list(self):
        return sorted(list(set([target.wg for target in self.targets])))
    
    @property
    def priorities(self):
        return sorted(list(set([target.priority for target in self.targets])))
    
    @property
    def names(self):
        return [target.name for target in self.targets]
    
    @property
    def wg_objects(self):
        # Create a dictionary with working group names as keys and target names as values
        wg_objects = {wg: [] for wg in self.wg_list}
        for target in self.targets:
            wg_objects[target.wg].append(target.name)
        # Sort the target names for each working group
        for wg in wg_objects:
            if wg == 'CO':
                wg_objects[wg] = sorted(wg_objects[wg], key=lambda x:natural_sort_key(x))
            else:
                wg_objects[wg].sort()
        # Return the dictionary
        return wg_objects

    @property
    def nexp_wg_finished(self):
        return {wg: sum(target.observed 
                        for target in self.targets 
                        if target.wg == wg) for wg in self.wg_list}
    
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
                return self.targets[index]
            else:
                raise IndexError("TargetList index out of range")
        else:
            raise TypeError("Index must be an integer")
        
    def update_observed(self, target_name, nexp):
        """
        ターゲットの観測数を更新します。

        Args:
            target_name (str): 更新するターゲットの名前。
            nexp (int): 更新する観測数。
        """
        index = self.get_index(target_name)
        if index is not None:
            target = self.targets[index]
            target.observed += nexp
            if target.observed >= target.nexp:
                self.completed_targets.append(target)
                self.observing_targets.remove(target)
        else:
            print(f"Target {target_name} not found in the target list.")


class TargetManager:
    def __init__(self, fname_targets, fname_targets_finish=None):
        
        target_table = Table.read(fname_targets, format='ascii')

        self.targetList = TargetList()

        for row in target_table:
            wg = row['wg']
            name = row['name']
            coord = SkyCoord(row['ra']+' '+row['dec'], unit=(u.hourangle, u.deg))
            pa = row['pa'] * u.deg
            nexp = row['nexp']
            priority = row['priority']
            observed = row['observed'] if 'observed' in row.colnames else 0
            target = Target(wg, name, coord, pa, nexp, priority, observed)
            self.targetList.add_target(target)

        if fname_targets_finish and os.path.exists(fname_targets_finish):
            target_table_finish = Table.read(fname_targets_finish, format='ascii')
            for row in target_table_finish:
                name = row['name']
                nexp = int(row['exptime'] / 900)
                if name.startswith('SSP_GA'):
                    nexp = min(nexp, 2)
                self.targetList.add_observed(name, nexp)
