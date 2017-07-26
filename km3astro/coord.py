"""Coordinate transformations.

Galactic:
    GC at (0, 0),
    gal. longitude, latitude (l, b)

Horizontal / altaz (km3):
    centered at detector position
    altitude, azimuth (altitude = 90deg - zenith)

EquatorialJ200 / FK5 / ICRS / GCRS
    (right ascension, declination)

    Equatorial is the same as FK5. FK5 is superseded by the ICRS, so use
    this instead. Note that FK5/ICRS are _barycentric_ implementations,
    so if you are looking for *geocentric* equatorial (i.e.
    for solar system bodies), use GCRS.


A note on maing conventions:
``phi`` and ``theta`` refer to neutrino directions, ``azimuth`` and
``zenith`` to source directions (i.e. the inversed neutrino direction).
The former says where the neutrino points to, the latter says where it comes
from.
"""
from astropy.units import rad, deg  # noqa
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, Longitude,
                                 Latitude, get_sun)
import numpy as np

from km3astro.constants import (
    arca_longitude, arca_latitude, arca_height,
    orca_longitude, orca_latitude, orca_height,
    antares_longitude, antares_latitude, antares_height,
)
from km3astro.time import np_to_astrotime
from km3astro.random import random_date, random_azimuth


LOCATIONS = {
    'arca': EarthLocation.from_geodetic(
        lon=Longitude(arca_longitude * deg),
        lat=Latitude(arca_latitude * deg),
        height=arca_height
    ),
    'orca': EarthLocation.from_geodetic(
        lon=Longitude(orca_longitude * deg),
        lat=Latitude(orca_latitude * deg),
        height=orca_height
    ),
    'antares': EarthLocation.from_geodetic(
        lon=Longitude(antares_longitude * deg),
        lat=Latitude(antares_latitude * deg),
        height=antares_height
    ),
}


def get_location(location='orca'):
    try:
        loc = LOCATIONS[location]
    except KeyError:
        raise KeyError("Invalid location, valid are 'orca', 'arca', 'antares'")
    return loc


def neutrino_to_source_direction(phi, theta, radian=True):
    """Flip the direction.

    Parameters
    ==========
    phi, theta: neutrino direction
    radian: bool [default=True]
        receive + return angles in radian? (if false, use degree)
    """
    phi = np.atleast_1d(phi)
    theta = np.atleast_1d(theta)
    if not radian:
        phi *= np.pi / 180
        theta *= np.pi / 180
    azimuth = (phi + np.pi) % (2 * np.pi)
    zenith = np.pi - theta
    if not radian:
        azimuth *= 180 / np.pi
        zenith *= 180 / np.pi
    return azimuth, zenith


def local_frame(time, location='orca'):
    loc = get_location(location)
    frame = AltAz(obstime=time, location=loc)
    return frame


def local_event(azimuth, time, zenith, radian=True,
                location='orca', **kwargs):
    """Create astropy events from detector coordinates."""
    time = np_to_astrotime(time)
    zenith = np.atleast_1d(zenith)
    azimuth = np.atleast_1d(azimuth)
    if not radian:
        azimuth *= np.pi / 180
        zenith *= np.pi / 180
    altitude = zenith - np.pi / 2
    frame = local_frame(time, location=location)
    event = SkyCoord(alt=altitude * rad, az=azimuth * rad, frame=frame,
                     **kwargs)
    return event


def sun_local(time, loc='orca'):
    time = np_to_astrotime(time)
    frame = local_frame(time, location='orca')
    sun = get_sun(time)
    sun_local = sun.transform_to(frame)
    return sun_local


def galactic_center():
    return SkyCoord(0 * deg, 0 * deg, frame='galactic')


def gc_in_local(time, loc='orca'):
    time = np_to_astrotime(time)
    frame = local_frame(time, location='orca')
    gc = galactic_center()
    gc_local = gc.transform_to(frame)
    return gc_local


def orca_gc_dist(azimuth, time, zenith, frame='detector'):
    """Return angular distance of event to GC.

    Parameters
    ==========
    frame: str, [default: 'detector']
        valid are 'detector', 'galactic', 'icrs', 'gcrs'
    """
    evt = local_event(azimuth, time, zenith)
    galcen = gc_in_local(time, loc='orca')
    if frame == 'detector':
        pass
    elif frame in ('galactic', 'icrs', 'gcrs'):
        evt = evt.transform_to(frame)
        galcen = galcen.transform_to(frame)
    return evt.separation(galcen).radian


def orca_sun_dist(azimuth, time, zenith):
    """Return distance of event to sun, in detector coordinates."""
    evt = local_event(azimuth, time, zenith)
    sun = sun_local(time, loc='orca')
    dist = evt.separation(sun).radian
    return dist


def gc_dist_random(zenith, frame='detector'):
    """Generate random (time, azimuth) events and get distance to GC."""
    n_evts = len(zenith)
    time = random_date(n=n_evts)
    azimuth = random_azimuth(n=n_evts)
    dist = orca_gc_dist(azimuth, time, zenith, frame=frame)
    return dist


def sun_dist_random(zenith):
    """Generate random (time, azimuth) events and get distance to GC."""
    n_evts = len(zenith)
    time = random_date(n=n_evts)
    azimuth = random_azimuth(n=n_evts)
    dist = orca_sun_dist(azimuth, time, zenith)
    return dist


def hsin(theta):
    """haversine"""
    return (1.0 - np.cos(theta)) / 2.


def space_angle(zen_1, zen_2, azi_1, azi_2):
    """Space angle between two directions specified by zenith and azimuth.
    """
    return hsin(azi_2 - azi_1) + numpy.cos(azi_1) * np.cos(azi_2) * hsin(zen_2 - zen_1)
