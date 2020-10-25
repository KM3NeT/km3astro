from unittest import TestCase

import numpy as np
from numpy.testing import assert_allclose

from km3astro.coord import (
    neutrino_to_source_direction,
    sun_local,
    convergence_angle,
    utm_zone,
    longitude_of_central_meridian,
)
from km3astro.random import random_date


class TestCoord(TestCase):
    def setUp(self):
        self.n_evts = 100
        self.n_evts_funny = 1e2

    def test_neutrino_flip_degree(self):
        phi = np.array([97.07, 23.46, 97.07, 192.5, 333.33])
        theta = np.array([135.0, 11.97, 22.97, 33.97, 85.23])
        azi_exp = np.array([277.07, 203.46, 277.07, 12.5, 153.33])
        zen_exp = np.array([45.0, 168.03, 157.03, 146.03, 94.77])
        azi, zen = neutrino_to_source_direction(phi, theta, radian=False)
        assert_allclose(azi, azi_exp)
        assert_allclose(zen, zen_exp)

    def test_neutrino_flip_radian(self):
        phi = np.array([97.07, 23.46, 97.07, 192.5, 333.33]) * np.pi / 180
        theta = np.array([135.0, 11.97, 22.97, 33.97, 85.23]) * np.pi / 180
        azi_exp = np.array([277.07, 203.46, 277.07, 12.5, 153.33]) * np.pi / 180
        zen_exp = np.array([45.0, 168.03, 157.03, 146.03, 94.77]) * np.pi / 180
        azi, zen = neutrino_to_source_direction(phi, theta, radian=True)

        assert_allclose(azi, azi_exp)
        assert_allclose(zen, zen_exp)


class TestCoordRandom(TestCase):
    def test_sun(self):
        date = random_date(n=100)
        sun = sun_local(date)


class TestConvergenceAngle(TestCase):
    def test_convergence_angle(self):
        ca = convergence_angle(1.3, 1.5)
        self.assertAlmostEqual(-0.00897440033130838, ca)


class TestUTMStuff(TestCase):
    def test_utm_zone(self):
        assert 38 == utm_zone(np.pi / 180 * 42.8871)

    def test_longitude_of_central_meridian(self):
        self.assertAlmostEqual(0.785398163397448, longitude_of_central_meridian(38))
