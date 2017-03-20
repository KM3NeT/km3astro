"""
===============================
Sun and GC in local coordinates
===============================

Show off some coordinate transformations.
"""

# Author: Moritz Lotze <mlotze@km3net.de>
# License: BSD-3

from astropy.units import rad, deg
from astropy.coordinates import get_sun, AltAz, SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from km3astro.time import random_date, np_to_astrotime
from km3astro.coord import (random_azimuth, random_zenith, ORCA_LOC,
                            orca_gc_dist)


#############################
# generate some random events

n_evts = 1e3
zen = random_zenith(n=n_evts)
tim = random_date(n=n_evts)
azi = random_azimuth(n=n_evts)
astrotim = np_to_astrotime(tim)


#############################

orca_frame = AltAz(obstime=astrotim, location=ORCA_LOC)
gc = SkyCoord(0 * rad, 0 * rad, frame='galactic')
sun = get_sun(astrotim)


sun_orca = sun.transform_to(orca_frame)
gc_orca = gc.transform_to(orca_frame)


#############################

sun_azi = sun_orca.az.rad
sun_zen = (90 * deg - sun_orca.alt).rad
gc_azi = gc_orca.az.rad
gc_zen = (90 * deg - gc_orca.alt).rad


############################

plt.hexbin(sun_zen, sun_azi)


############################

plt.hexbin(np.cos(sun_zen), sun_azi)


############################

plt.hexbin(gc_zen, gc_azi)


############################

plt.hexbin(np.cos(gc_zen), gc_azi)


# In[10]:

df = pd.DataFrame({
    'dist_gc': orca_gc_dist(azi, tim, zen, frame='gcrs'),
    'dist_sun': orca_gc_dist(azi, tim, zen, frame='gcrs'),
})

############################

df.dist_gc.hist()


############################

df.dist_sun.hist()
