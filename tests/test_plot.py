from unittest import TestCase

import pandas as pd
from km3net_testdata import data_path
import km3astro.toolbox as ktb
import km3astro.testing_tools as ktt
import km3astro.coord as kc
import km3astro.plot as kp


class TestPlotSkymap(TestCase):
    def test_skymap_list(self):

        _ = kp.skymap_list(
            file0=data_path("astro/antares_coordinate_systems_benchmark.csv"),
            frame="UTM",
            detector="antares",
            plot_frame="galactic",
            detector_to="antares",
            save=True,
        )

        table_read = pd.read_csv(
            data_path("astro/antares_coordinate_systems_benchmark.csv"), comment="#"
        )
        alert_type = [
            "GRB",
            "GW",
            "Neutrino",
            "NuEM",
            "SK_SN",
            "SNEWS",
            "Transient",
            "Random",
            "GRB",
            "GW",
            "Neutrino",
            "NuEM",
            "SK_SN",
            "SNEWS",
            "Transient",
            "Random",
            "Random",
        ]

        table_read["Alert_type"] = alert_type

        _ = kp.skymap_list(
            dataframe=table_read,
            frame="UTM",
            detector="antares",
            plot_frame="equatorial",
            detector_to="antares",
            title="test_title",
            save=True,
            name="test_dataframe_input",
        )

        _ = kp.skymap_list(
            file0=data_path("astro/ORCA_coordinate_systems_benchmark.csv"),
            frame="UTM",
            detector="orca",
            plot_frame="equatorial",
            detector_to="orca",
        )

        _ = kp.skymap_list(
            file0=data_path("astro/ORCA_coordinate_systems_benchmark.csv"),
            frame="UTM",
            detector="orca",
            plot_frame="galactic",
            detector_to="orca",
        )

        _ = kp.skymap_list(
            file0=data_path("astro/ARCA_coordinate_systems_benchmark.csv"),
            frame="UTM",
            detector="arca",
            plot_frame="equatorial",
            detector_to="arca",
        )

        _ = kp.skymap_list(
            file0=data_path("astro/ARCA_coordinate_systems_benchmark.csv"),
            frame="UTM",
            detector="arca",
            plot_frame="galactic",
            detector_to="arca",
        )

    def test_skymap_alert(self):

        _ = kp.skymap_alert(
            file0=data_path("astro/antares_coordinate_systems_benchmark.csv"),
            frame="UTM",
            detector="antares",
            plot_frame="ParticleFrame",
            detector_to="antares",
            save=True,
        )

        _ = kp.skymap_alert(
            ra=80,
            dec=-20,
            obstime="2022-07-18T03:03:03",
            plot_frame="galactic",
            detector="orca",
            detector_to="orca",
        )

        _ = kp.skymap_alert(
            ra=80,
            dec=-20,
            obstime="2022-07-18T03:03:03",
            plot_frame="equatorial",
            detector="dummy",
            detector_to="orca",
        )

    def test_skymap_hpx(self):

        _ = kp.skymap_hpx(
            file0="https://gracedb.ligo.org/api/superevents/MS230522k/files/bayestar.multiorder.fits,1",
            save=True,
        )
