from astropy.time import Time
from astropy.units import rad, deg  # noqa
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, Longitude,
                                 Latitude)
import numpy as np

from km3astro.constants import orca_longitude, orca_latitude, orca_height
from km3astro.time import np_to_datetime, random_date


def orca_event(azimuth, time, zenith):
    """Create astropy events from detector coordinates."""
    time = np.atleast_1d(time)
    zenith = np.atleast_1d(zenith)
    azimuth = np.atleast_1d(azimuth)

    orca_loc = EarthLocation.from_geodetic(
        Longitude(orca_longitude * deg),
        Latitude(orca_latitude * deg),
        height=orca_height
    )
    time = Time(np_to_datetime(time))
    orca_frame = AltAz(obstime=time, location=orca_loc)

    altitude = zenith - np.pi/2
    event = SkyCoord(alt=altitude*rad, az=azimuth*rad, frame=orca_frame)
    return event


def gc_dist(event):
    """Returns distance to galactic center, in ICRS coordinates.

    If no events (SkyCoord
    """
    events = event.icrs
    gc = SkyCoord(0*deg, 0*deg, frame='galactic').icrs
    return events.separation(gc).radian


def orca_gc_dist(azimuth, time, zenith):
    evt = orca_event(azimuth, time, zenith)
    dist = gc_dist(evt)
    return dist


def random_azimuth(n=1, unit='rad'):
    azi = (np.pi * np.random.random_sample(size=n))
    if unit == 'rad':
        return azi
    elif unit == 'deg':
        return azi / np.pi * 180
    else:
        raise KeyError("Unknown unit '{}'".format(unit))


def gc_dist_random(zenith, time='random', azimuth='random'):
    n_evts = len(zenith)
    if time == 'random':
        time = random_date(n=n_evts)
    if azimuth == 'random':
        azimuth = random_azimuth(n=n_evts)
    dist = orca_gc_dist(azimuth, time, zenith)
    return dist
