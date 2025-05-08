import yaml
import astropy.units as u

class Params:
    def __init__(self, yaml_file):
        with open(yaml_file, 'r') as file:
            self.params = yaml.safe_load(file)

    def __getitem__(self, key):
        return self.params.get(key, None)
    
    @property
    def fname_obsdate(self):
        return self.params.get('fname_obsdate', None)
    
    @property
    def fname_obsdate_finish(self):
        return self.params.get('fname_obsdate_finish', None)
    
    @property
    def fname_targets(self):
        return self.params.get('fname_targets', None)
    
    @property
    def fname_targets_finish(self):
        return self.params.get('fname_targets_finish', None)
    
    @property
    def fname_report(self):
        return self.params.get('fname_report', None)
    
    @property
    def frac(self):
        return self.params.get('frac', None)
    
    @property
    def frac_margin(self):
        return self.params.get('frac_margin', None)
    
    @property
    def GA_last(self):
        return self.params.get('GA_last', None)
    
    @property
    def n_continuous(self):
        return self.params.get('n_continuous', None)
    
    @property
    def w_timeslot(self):
        return self.params.get('w_timeslot', None) * u.minute
    
    @property
    def t_overhead(self):
        return self.params.get('t_overhead', None) * u.minute
    
    @property
    def angle_twilight(self):
        return self.params.get('angle_twilight', None) * u.degree
    
    @property
    def airmass(self):
        return self.params.get('airmass', None) 
    
    @property
    def meridian(self):
        _ = self.params.get('meridian', None)
        return {key: value * u.hourangle for key, value in _.items()}
    
    @property
    def moonsep(self):
        _ = self.params.get('moonsep', None)
        return {key: value * u.degree for key, value in _.items()}
    
    @property
    def moonill(self):
        return self.params.get('moonill', None)
    
    @property
    def moonalt(self):
        _ = self.params.get('moonalt', None)
        return {key: value * u.degree for key, value in _.items()}

    @property
    def planetssep(self):
        _ = self.params.get('planetssep', None)
        return {key: value * u.degree for key, value in _.items()}
    
    @property
    def weight_comp(self):
        return self.params.get('weight_comp', None)
    
    @property
    def weight_pri(self):
        return self.params.get('weight_pri', None)

    @property
    def slew_speed_az(self):
        return self.params.get('slew_speed_az', None) * u.degree / u.second
    
    @property
    def slew_speed_el(self):
        return self.params.get('slew_speed_el', None) * u.degree / u.second
    
    @property
    def inst_rot_speed(self):
        return self.params.get('inst_rot_speed', None) * u.degree / u.second