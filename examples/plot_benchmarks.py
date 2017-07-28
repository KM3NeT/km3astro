"""
======================================================
Benchmark Calculations for Coordinates Transformations
======================================================

Run the benchmarks similar to
`http://antares.in2p3.fr/internal/dokuwiki/doku.php?id=benchmarks_astro`

BIG CAVEAT:
What IceCube and ANTARES/KM3NeT call "azimuth" is actually the
co-azimuth, or ``(90 - azimuth) % 360``!
"""

# Author: Moritz Lotze <mlotze@km3net.de>
# License: BSD-3

from collections import OrderedDict

import numpy as np
import pandas as pd

from km3astro.coord import (neutrino_to_source_direction,
                            local_event, local_frame)
from km3astro.sources import SIRIUS, CANOPUS, ARCTURUS, ANTARES


time = pd.to_datetime([
    '2007-10-04 03:03:03.00',
    '2007-10-04 03:03:03.00',
    '2007-10-04 03:03:03.00',
    '2007-10-04 03:03:03.00',
    '2007-10-04 03:03:03.00',
])
theta_deg = np.array([
    135.00,
    11.97,
    22.97,
    33.97,
    85.23,
])
phi_deg = np.array([
    97.07,
    23.46,
    97.07,
    192.50,
    333.33,
])
theta = theta_deg * np.pi / 180
phi = phi_deg * np.pi / 180
azimuth, zenith = neutrino_to_source_direction(phi, theta)

zenith = np.pi - zenith
azimuth = np.pi + azimuth

evt = local_event(azimuth, time, zenith, location='antares')
equat = evt.fk5
#dec = equat.dec.degree
#ra = equat.ra.degree
dec = equat.dec
ra = equat.ra
gal = evt.galactic
l = gal.l
b = gal.b

data = OrderedDict([
    ('time', time),
    ('theta [deg]', theta_deg),
    ('phi [deg]', phi_deg),
    ('theta', theta),
    ('phi', phi),
    ('zenith [deg]', zenith * 180 / np.pi),
    ('azimuth [deg]', azimuth * 180 / np.pi),
    ('zenith', zenith),
    ('azimuth', azimuth),
    ('declication', dec),
    ('right ascension', ra),
    ('gal_longitude', l),
    ('gal_latitude', b),
])
bench = pd.DataFrame(data=data)
print(bench[:])

########################################################
# bla bla

frame = local_frame(time[0], location='antares')
sirius_local = SIRIUS.transform_to(frame)
sirius_alt = sirius_local.alt.degree
sirius_az = sirius_local.az.degree
sirius_zen = 90 - sirius_alt
print()
print(sirius_az)
print(sirius_zen)

canopus_local = CANOPUS.transform_to(frame)
canopus_alt = canopus_local.alt.degree
canopus_az = canopus_local.az.degree
canopus_zen = 90 - canopus_alt
print()
print(canopus_az)
print(canopus_zen)
