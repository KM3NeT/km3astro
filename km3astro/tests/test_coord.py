from unittest import TestCase

import numpy as np
import pandas as pd

from km3astro.coord import gc_dist, orca_event, orca_gc_dist


class TestGalCenSep(TestCase):
    def setUp(self):
        rand_t = pd.Series(np.array(['2015-08-09T23:18:18', '2015-06-23T03:53:35',
                           '2015-11-19T22:44:57', '2015-10-30T10:01:46',
                           '2015-11-08T02:49:24'], dtype='datetime64[s]'))
        self.X = pd.DataFrame({
            'azimuth': {0: 4.7519763176762861, 1: 2.4658551815113445,
                        2: 0.16061711736367185, 3: 2.8638258582673122,
                        4: 4.327783960628202},
            'time': rand_t,
            'zenith': {0: 1.5058163860685134, 1: 1.395836281144522,
                       2: 0.8759613758638688, 3: 2.0804470026610935,
                       4: 0.54816919871933478}
        })
        self.event = orca_event(self.X.azimuth, self.X.time, self.X.zenith)

    def test_distance(self):
        expected = np.array([ 0.5633849 ,  1.69012499,  1.28103244,
                             0.42668377,  1.21924223])
        dist = gc_dist(self.event)
        assert np.allclose(expected, dist)
