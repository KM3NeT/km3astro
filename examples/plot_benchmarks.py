"""
======================================================
Benchmark Calculations for Coordinates Transformations
======================================================

Run the benchmarks similar to
`http://antares.in2p3.fr/internal/dokuwiki/doku.php?id=benchmarks_astro`
"""

# Author: Moritz Lotze <mlotze@km3net.de>
# License: BSD-3

import numpy as np
import pandas as pd

from km3astro.coord import neutrino_to_source_direction


time = pd.to_datetime([
    '2007-10-04 03:03:03.00',
    '2007-10-04 03:03:03.00',
    '2007-10-04 03:03:03.00',
    '2007-10-04 03:03:03.00',
    '2007-10-04 03:03:03.00',
])
theta = np.array([
    135.00,
    11.97,
    22.97,
    33.97,
    85.23,
])
phi = np.array([
    97.07,
    23.46,
    97.07,
    192.50,
    333.33,
])
azimuth, zenith = neutrino_to_source_direction(phi, theta, radian=False)

bench = pd.DataFrame({
    'time': time,
    'theta': theta,
    'phi': phi,
})
