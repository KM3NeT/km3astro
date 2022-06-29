from unittest import TestCase

from km3net_testdata import data_path
import km3astro.toolbox as ktb
import km3astro.testing_tools as ktt
import km3astro.coord as kc


class TestToolBox(TestCase):
    def setUp(self):
        self.det = "antares"
        self.date = "2007-10-04"
        self.time = "03:03:03"
        self.phi = 97.07
        self.theta = 135.0
        self.az = 277.07
        self.ze = 45.00
        self.ra = 70.613
        self.dec = -1.852
        self.l = 198.7
        self.b = -29.298
        self.threshold = 0.02

    def test_tool_box_eq_to_utm(self):

        SC_check = ktb.eq_to_utm(self.ra, self.dec, self.date, self.time, self.det)
        SC_true = kc.build_event(
            "UTM", self.date, self.time, self.az, self.ze, "deg", self.det
        )
        sep = ktt.test_Skycoord_separation(SC_true, SC_check)
        if self.threshold < sep:
            raise AssertionError(
                "Error: Maximum angle separation = "
                + str(sep)
                + " > "
                + str(self.threshold)
            )

    def test_tool_box_eq_to_loc(self):

        SC_check = ktb.eq_to_loc(self.ra, self.dec, self.date, self.time, self.det)
        SC_true = kc.build_event(
            "ParticleFrame", self.date, self.time, self.theta, self.phi, "deg", self.det
        )
        sep = ktt.test_Skycoord_separation(SC_true, SC_check)
        if self.threshold < sep:
            raise AssertionError(
                "Error: Maximum angle separation = "
                + str(sep)
                + " > "
                + str(self.threshold)
            )

    def test_tool_box_eq_to_gal(self):

        SC_check = ktb.eq_to_gal(self.ra, self.dec, self.date, self.time)
        SC_true = kc.build_event(
            "galactic", self.date, self.time, self.l, self.b, "deg"
        )
        sep = ktt.test_Skycoord_separation(SC_true, SC_check)
        if self.threshold < sep:
            raise AssertionError(
                "Error: Maximum angle separation = "
                + str(sep)
                + " > "
                + str(self.threshold)
            )

    def test_tool_box_loc_to_eq(self):

        SC_check = ktb.loc_to_eq(self.phi, self.theta, self.date, self.time, self.det)
        SC_true = kc.build_event(
            "equatorial", self.date, self.time, self.ra, self.dec, "deg"
        )
        sep = ktt.test_Skycoord_separation(SC_true, SC_check)
        if self.threshold < sep:
            raise AssertionError(
                "Error: Maximum angle separation = "
                + str(sep)
                + " > "
                + str(self.threshold)
            )

    def test_tool_box_loc_to_utm(self):

        SC_check = ktb.loc_to_utm(self.phi, self.theta, self.date, self.time, self.det)
        SC_true = kc.build_event(
            "UTM", self.date, self.time, self.az, self.ze, "deg", self.det
        )
        sep = ktt.test_Skycoord_separation(SC_true, SC_check)
        if self.threshold < sep:
            raise AssertionError(
                "Error: Maximum angle separation = "
                + str(sep)
                + " > "
                + str(self.threshold)
            )

    def test_tool_box_loc_to_gal(self):

        SC_check = ktb.loc_to_gal(self.phi, self.theta, self.date, self.time, self.det)
        SC_true = kc.build_event(
            "galactic", self.date, self.time, self.l, self.b, "deg"
        )
        sep = ktt.test_Skycoord_separation(SC_true, SC_check)
        if self.threshold < sep:
            raise AssertionError(
                "Error: Maximum angle separation = "
                + str(sep)
                + " > "
                + str(self.threshold)
            )
