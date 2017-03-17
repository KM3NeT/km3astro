"""Coordinate transformations.
"""
from astropy.time import Time
from astropy.units import rad, deg, hourangle  # noqa
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, Longitude,
                                 Latitude, get_sun, get_moon)
import numpy as np

from km3astro.constants import (
    arca_longitude, arca_latitude, arca_height,
    orca_longitude, orca_latitude, orca_height,
)
from km3astro.time import np_to_datetime, random_date, np_to_astrotime


ARCA_LOC = EarthLocation.from_geodetic(
    Longitude(arca_longitude * deg),
    Latitude(arca_latitude * deg),
    height=arca_height
)

ORCA_LOC = EarthLocation.from_geodetic(
    Longitude(orca_longitude * deg),
    Latitude(orca_latitude * deg),
    height=orca_height
)


def transform_to_orca(event, time):
    time = np_to_astrotime(time)
    orca_frame = AltAz(obstime=time, location=ORCA_LOC)
    return event.transform_to(orca_frame)


def orca_event(azimuth, time, zenith):
    """Create astropy events from detector coordinates."""
    zenith = np.atleast_1d(zenith)
    azimuth = np.atleast_1d(azimuth)

    time = np_to_astrotime(time)
    orca_frame = AltAz(obstime=time, location=ORCA_LOC)

    altitude = zenith - np.pi/2
    event = SkyCoord(alt=altitude*rad, az=azimuth*rad, frame=orca_frame)
    return event


def random_azimuth(n=1, unit='rad'):
    """Draw azimuth, uniformly distributed."""
    azi = (2 * np.pi * np.random.random_sample(size=n))
    if unit == 'rad':
        return azi
    elif unit == 'deg':
        return azi / np.pi * 180
    else:
        raise KeyError("Unknown unit '{}'".format(unit))


def random_zenith(n=1, unit='rad'):
    """Draw zenith, uniformly distributed in cos(zen)."""
    coszen = (2 * np.random.random_sample(size=n) - 1)
    zen = np.arccos(coszen)
    if unit == 'rad':
        return zen
    elif unit == 'deg':
        return zen / np.pi * 180
    else:
        raise KeyError("Unknown unit '{}'".format(unit))


def sun_in_orca(time):
    time = np_to_astrotime(time)
    orca_frame = AltAz(obstime=time, location=ORCA_LOC)
    sun = get_sun(time)
    sun_orca = sun.transform_to(orca_frame)
    return sun_orca


def gc_in_orca(time):
    time = np_to_astrotime(time)
    orca_frame = AltAz(obstime=time, location=ORCA_LOC)
    gc = SkyCoord(0*deg, 0*deg, frame='galactic')
    gc_orca = gc.transform_to(orca_frame)
    return gc_orca
