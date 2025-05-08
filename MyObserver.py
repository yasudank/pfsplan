from astroplan import Observer
import astropy.units as u
from datetime import datetime
from astropy.coordinates import EarthLocation

class MyObserver(Observer):

    """
    A custom observer class that extends astroplan.Observer to include 
    the UTC offset calculation support.
    This class allows for the specification of a timezone, which is used to calculate
    the UTC offset for the observer's location.
    """

    def __init__(self, *args, **kwargs):
        """
        Initializes MyObserver, inheriting from astroplan.Observer.

        Parameters
        ----------
        *args : tuple
            Positional arguments passed to astroplan.Observer.__init__.
        timezone : str, optional
            Timezone name (e.g., 'America/Los_Angeles'). If None, defaults to UTC.
        **kwargs : dict
            Keyword arguments passed to astroplan.Observer.__init__.
        """
        super().__init__(*args, **kwargs)  # Initialize the astroplan.Observer part

        _timezone = kwargs.pop('timezone', None)  # Extract timezone from kwargs
        
        self._timezone = _timezone  # Store the timezone

        if _timezone is not None:
            # Calculate the UTC offset in hours
            self._utcoffset = (datetime.now().astimezone(self.timezone).utcoffset().total_seconds() * u.s).to(u.hour)
        else:
            self._utcoffset = 0 * u.hour  # Default to UTC

    @property
    def utcoffset(self):
        """
        The UTC offset of the observer's timezone in hours.
        """
        return self._utcoffset
