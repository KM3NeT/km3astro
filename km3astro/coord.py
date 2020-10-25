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

Also radian is the default. Degree can be used, but generally the default is
to assume radian.
"""
from astropy import units as u
from astropy.units import rad, deg, hourangle  # noqa
from astropy.coordinates import (
    EarthLocation,
    SkyCoord,
    AltAz,
    Longitude,
    Latitude,
    get_sun,
    get_moon,
)
import astropy.time
import numpy as np

from km3astro.constants import (
    arca_longitude,
    arca_latitude,
    arca_height,
    orca_longitude,
    orca_latitude,
    orca_height,
    antares_longitude,
    antares_latitude,
    antares_height,
)
from km3astro.time import np_to_astrotime
from km3astro.random import random_date, random_azimuth
from km3astro.sources import GALACTIC_CENTER


LOCATIONS = {
    "arca": EarthLocation.from_geodetic(
        lon=Longitude(arca_longitude * deg),
        lat=Latitude(arca_latitude * deg),
        height=arca_height,
    ),
    "orca": EarthLocation.from_geodetic(
        lon=Longitude(orca_longitude * deg),
        lat=Latitude(orca_latitude * deg),
        height=orca_height,
    ),
    "antares": EarthLocation.from_geodetic(
        lon=Longitude(antares_longitude * deg),
        lat=Latitude(antares_latitude * deg),
        height=antares_height,
    ),
}


def neutrino_to_source_direction(phi, theta, radian=True):
    """Flip the direction.

    Parameters
    ----------
    phi, theta: neutrino direction
    radian: bool [default=True]
        receive + return angles in radian? (if false, use degree)

    """
    phi = np.atleast_1d(phi).copy()
    theta = np.atleast_1d(theta).copy()
    if not radian:
        phi *= np.pi / 180
        theta *= np.pi / 180
    assert np.all(phi <= 2 * np.pi)
    assert np.all(theta <= np.pi)
    azimuth = (phi + np.pi) % (2 * np.pi)
    zenith = np.pi - theta
    if not radian:
        azimuth *= 180 / np.pi
        zenith *= 180 / np.pi
    return azimuth, zenith


def source_to_neutrino_direction(azimuth, zenith, radian=True):
    """Flip the direction.

    Parameters
    ----------
    zenith : float
        neutrino origin
    azimuth: float
        neutrino origin
    radian: bool [default=True]
        receive + return angles in radian? (if false, use degree)

    """
    azimuth = np.atleast_1d(azimuth).copy()
    zenith = np.atleast_1d(zenith).copy()
    if not radian:
        azimuth *= np.pi / 180
        zenith *= np.pi / 180
    phi = (azimuth - np.pi) % (2 * np.pi)
    theta = np.pi - zenith
    if not radian:
        phi *= 180 / np.pi
        theta *= 180 / np.pi
    return phi, theta


def get_location(location="orca"):
    try:
        loc = LOCATIONS[location]
    except KeyError:
        raise KeyError("Invalid location, valid are 'orca', 'arca', 'antares'")
    return loc


def Sun(time):
    """Wrapper around astropy's get_sun, accepting numpy/pandas time objects."""
    if not isinstance(time, astropy.time.Time):
        # if np.datetime64, convert to astro time
        time = np_to_astrotime(time)
    return get_sun(time)


def Moon(time):
    """Wrapper around astropy's get_moon, accepting numpy/pandas time objects."""
    if not isinstance(time, astropy.time.Time):
        # if np.datetime64, convert to astro time
        time = np_to_astrotime(time)
    return get_moon(time)


def local_frame(time, location="orca"):
    """Get the (horizontal) coordinate frame of your detector."""
    if not isinstance(time, astropy.time.Time):
        # if np.datetime64, convert to astro time
        time = np_to_astrotime(time)
    loc = get_location(location)
    frame = AltAz(obstime=time, location=loc)
    return frame


def local_event(azimuth, time, zenith, radian=True, location="orca", **kwargs):
    """Create astropy events from detector coordinates."""
    zenith = np.atleast_1d(zenith).copy()
    azimuth = np.atleast_1d(azimuth).copy()
    if not radian:
        azimuth *= np.pi / 180
        zenith *= np.pi / 180
    altitude = zenith - np.pi / 2
    # neutrino telescopes call the co-azimuth "azimuth"
    true_azimuth = (np.pi / 2 - azimuth) % np.pi
    frame = local_frame(time, location=location)
    event = SkyCoord(alt=altitude * rad, az=true_azimuth * rad, frame=frame, **kwargs)
    return event


def sun_local(time, loc="orca"):
    """Sun position in local coordinates."""
    frame = local_frame(time, location=loc)
    sun = Sun(time)
    sun_local = sun.transform_to(frame)
    return sun_local


def moon_local(time, loc="orca"):
    """Moon position in local coordinates."""
    frame = local_frame(time, location=loc)
    moon = Moon(time)
    moon_local = moon.transform_to(frame)
    return moon_local


def gc_in_local(time, loc="orca"):
    """Galactic center position in local coordinates."""
    frame = local_frame(time, location=loc)
    gc = GALACTIC_CENTER
    gc_local = gc.transform_to(frame)
    return gc_local


def orca_gc_dist(azimuth, time, zenith, frame="detector"):
    """Return angular distance of event to GC.

    Parameters
    ==========
    frame: str, [default: 'detector']
        valid are 'detector', 'galactic', 'icrs', 'gcrs'
    """
    evt = local_event(azimuth, time, zenith)
    galcen = gc_in_local(time, loc="orca")
    if frame == "detector":
        pass
    elif frame in ("galactic", "icrs", "gcrs"):
        evt = evt.transform_to(frame)
        galcen = galcen.transform_to(frame)
    return evt.separation(galcen).radian


def orca_sun_dist(azimuth, time, zenith):
    """Return distance of event to sun, in detector coordinates."""
    evt = local_event(azimuth, time, zenith)
    sun = sun_local(time, loc="orca")
    dist = evt.separation(sun).radian
    return dist


def gc_dist_random(zenith, frame="detector"):
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


class Event(object):
    def __init__(self, zenith, azimuth, time, location="orca"):
        self.zenith = zenith
        self.azimuth = azimuth
        self.time = time

    @classmethod
    def from_zenith(cls, zenith, **initargs):
        zenith = np.atleast_1d(zenith)
        n_evts = zenith.shape[0]
        azimuth = random_azimuth(n_evts)
        time = random_date(n_evts)
        return cls(zenith, azimuth, time, **initargs)


def convergence_angle(lon, lat):
    """Calculate the converge angle on the UTM grid.

    Parameters
    ----------
    lon : number
        Longitude in rad
    lat : number
        Latitude in rad

    """
    latitude_deg = lat * u.deg

    if latitude_deg > 84 * u.deg or latitude_deg < -80 * u.deg:
        raise ValueError(
            "UTM coordinate system is only defined between -80deg S and 84deg N."
        )

    # detector position, longitude and latitude in rad
    # lambda  = longitude
    phi = lat

    # find UTM zone and central meridian

    # longitude of the central meridian of UTM zone in rad
    lambda0 = longitude_of_central_meridian(utm_zone(lon))
    omega = lon - lambda0

    # parameters of the Earth ellipsoid
    sma = 6378137  # semi-major axis in meters (WGS84)
    ecc = 0.0066943800  # eccentricity (WGS84)

    rho = sma * (1 - ecc) / pow(1 - ecc * np.sin(phi) ** 2, 3 / 2)
    nu = sma / np.sqrt(1 - ecc * np.sin(phi) ** 2)
    psi = nu / rho
    t = np.tan(phi)

    angle = (
        np.sin(phi) * omega
        - np.sin(phi) * omega ** 3 / 3 * pow(np.cos(phi), 2) * (2 * psi ** 2 - psi)
        - np.sin(phi)
        * omega ** 5
        / 15
        * pow(np.cos(phi), 4)
        * (
            psi ** 4 * (11 - 24 * t ** 2)
            - psi ** 3 * (11 - 36 * t ** 2)
            + 2 * psi ** 2 * (1 - 7 * t ** 2)
            + psi * t ** 2
        )
        - np.sin(phi)
        * omega ** 7
        / 315
        * pow(np.cos(phi), 6)
        * (17 - 26 * t ** 2 + 2 * t ** 4)
    )

    return angle


def utm_zone(lat):
    """The UTM zone for a given latitude

    Parameters
    ----------
    lat : number
        Latitude in rad

    """
    return 1 + int((np.pi + lat) / (6 * np.pi / 180))


def longitude_of_central_meridian(utmzone):
    """The longitude of the central meridian for a given UTM zone.

    Parameters
    ----------
    utmzone : number
        The UTM zone.

    """
    zone_width = 6 * np.pi / 180
    return -np.pi + (utmzone - 1) * zone_width + zone_width / 2
