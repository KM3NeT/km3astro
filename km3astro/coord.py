from astropy.time import Time
from astropy.units import rad, deg  # noqa
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, Longitude,
                                 Latitude)
import numpy as np
from km3astro.constants import orca_longitude, orca_latitude, orca_height


def orca_event(time, altitude, azimuth):
    """Create astropy events from detector coordinates."""
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
