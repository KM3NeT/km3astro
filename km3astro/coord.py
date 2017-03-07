from astropy.time import Time
from astropy.units import rad, deg  # noqa
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, Longitude,
                                 Latitude)
from datetime import datetime
import numpy as np


def second_from_interval(start, stop, n=1):
    """Sample random times from an interval (in seconds)."""
    sec = np.timedelta64(1, 's')
    n_seconds = (stop - start) / sec
    samples = np.random.randint(low=0, high=n_seconds, size=n) * sec
    return start + samples


def equidistant_from_interval(start, stop, step=np.timedelta64(2, 'm')):
    start = np.datetime64(start)
    stop = np.datetime64(stop)
    duration = stop - start
    n_steps = np.ceil(duration/step)
    samples = np.arange(n_steps) * step
    return start + samples


def random_date(year=2015, **randargs):
    """Create random dates in the given year."""
    start = np.datetime64('{}-01-01'.format(year))
    stop = np.datetime64('{}-01-01'.format(year+1))
    return second_from_interval(start, stop, **randargs)


def orca_event(time, altitude, azimuth):
    """Create astropy events from detector coordinates."""
    from km3astro.constants import orca_longitude, orca_latitude, orca_height
    orca_loc = EarthLocation.from_geodetic(
        Longitude(orca_longitude * deg),
        Latitude(orca_latitude * deg),
        height=orca_height
    )
    time = Time(np.atleast_1d(time).tolist())
    orca_frame = AltAz(obstime=time, location=orca_loc)
    altitude *= rad
    azimuth *= rad
    event = SkyCoord(alt=altitude, az=azimuth, frame=orca_frame)
    return event


def gc_dist(event):
    """Returns distance to galactic center, in ICRS coordinates.

    If no events (SkyCoord
    """
    events = event.icrs
    gc = SkyCoord(0*deg, 0*deg, frame='galactic').icrs
    return events.separation(gc)


def np_to_datetime(time):
    return datetime.utcfromtimestamp((time - np.datetime64(
        '1970-01-01T00:00:00Z')) / np.timedelta64(1, 's'))
