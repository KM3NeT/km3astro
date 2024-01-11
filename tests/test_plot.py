from unittest import TestCase, skipIf

import pandas as pd
from km3net_testdata import data_path
import km3astro.toolbox as ktb
import km3astro.testing_tools as ktt
import km3astro.coord as kc
import km3astro.plot as kp

import sys


class TestPlotSkymap(TestCase):
    @skipIf(sys.version_info < (3, 8), "ligo.skymap requires Python 3.8+")
    def test_skymap_list(self):

        _ = kp.skymap_list(
            frame="equatorial",
            frame_input="UTM",
            detector="antares",
        )

        table_read = pd.read_csv(
            data_path("astro/antares_coordinate_systems_benchmark.csv"), comment="#"
        )
        alert_type = [
            "GRB",
            "Transient",
            "Neutrino",
            "NuEM",
            "GRB",
            "GRB",
            "Transient",
            "Neutrino",
            "GRB",
            "GRB",
            "Neutrino",
            "Neutrino",
            "GRB",
            "GRB",
            "Transient",
            "GRB",
            "NuEM",
        ]
        table_read["Alert_type"] = alert_type

        _ = kp.skymap_list(
            dataframe=table_read,
            frame="equatorial",
            frame_input="galactic",
            detector="antares",
        )

        _ = kp.skymap_list(
            dataframe=table_read,
            frame="galactic",
            frame_input="UTM",
            detector="orca",
        )

    @skipIf(sys.version_info < (3, 8), "ligo.skymap requires Python 3.8+")
    def test_skymap_alert(self):
        _ = kp.skymap_alert(
            ra=80,
            dec=-20,
            obstime="2022-07-18T03:03:03",
            frame="equatorial",
            detector="orca",
        )
        _ = kp.skymap_alert(
            ra=80,
            dec=-20,
            error_radius=5,
            obstime="2022-07-18T03:03:03",
            frame="galactic",
            detector="antares",
        )
