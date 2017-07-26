from unittest import TestCase

import numpy as np
from numpy.testing import assert_allclose

from km3astro.coord import neutrino_to_source_direction


class TestCoord(TestCase):
    def setUp(self):
        self.n_evts = 100
        self.n_evts_funny = 1e2

    def test_neutrino_flip(self):
        phi_deg = np.array([97.07, 23.46, 97.07, 192.5, 333.33])
        theta_deg = np.array([135., 11.97, 22.97, 33.97, 85.23])
        azimuth_exp_deg = np.array([277.07, 203.46, 277.07, 12.5, 153.33])
        zenith_exp_deg = np.array([45., 168.03, 157.03, 146.03, 94.77])
        azi_deg, zen_deg = neutrino_to_source_direction(phi_deg, theta_deg,
                                                        radian=False)
        assert_allclose(zen_deg, zenith_exp_deg)
        assert_allclose(azi_deg, azimuth_exp_deg)
        phi = phi_deg * np.pi / 180
        theta = theta_deg * np.pi / 180
        azimuth_exp = azimuth_exp_deg * np.pi / 180
        zenith_exp = zenith_exp_deg * np.pi / 180
        azi, zen = neutrino_to_source_direction(phi, theta, radian=True)
        assert_allclose(zen, zenith_exp)
        assert_allclose(azi, azimuth_exp)


class TestEvent(TestCase):
    pass
