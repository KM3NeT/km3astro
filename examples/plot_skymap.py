"""
========================
Multiple exemple of how to use km3astro.plot
========================

"""

# Author: Tedjditi Hichem <htedjditi@km3net.de>


import numpy as np
import pandas as pd

import km3astro.plot as kp
from km3net_testdata import data_path


def main():

    fig = kp.skymap_list(
        data_path("astro/antares_coordinate_systems_benchmark.csv"),
        frame="UTM",
        detector="antares",
        plot_frame="equatorial",
        detector_to="antares",
        save=True,
    )

    fig.savefig("test_save_fig.png")

    table_read = pd.read_csv(
        data_path("astro/antares_coordinate_systems_benchmark.csv"), comment="#"
    )

    _ = kp.skymap_list(
        dataframe=table_read,
        frame="UTM",
        detector="antares",
        plot_frame="equatorial",
        detector_to="antares",
        title="test_title_input",
        save=True,
        name="test_dataframe_input",
    )

    _ = kp.skymap_list(
        data_path("astro/antares_coordinate_systems_benchmark.csv"),
        frame="UTM",
        detector="antares",
        plot_frame="galactic",
        detector_to="antares",
    )

    _ = kp.skymap_list(
        data_path("astro/ORCA_coordinate_systems_benchmark.csv"),
        frame="UTM",
        detector="orca",
        plot_frame="equatorial",
        detector_to="orca",
    )

    _ = kp.skymap_list(
        data_path("astro/ORCA_coordinate_systems_benchmark.csv"),
        frame="UTM",
        detector="orca",
        plot_frame="galactic",
        detector_to="orca",
    )

    _ = kp.skymap_list(
        data_path("astro/ARCA_coordinate_systems_benchmark.csv"),
        frame="UTM",
        detector="arca",
        plot_frame="equatorial",
        detector_to="arca",
    )

    _ = kp.skymap_list(
        data_path("astro/ARCA_coordinate_systems_benchmark.csv"),
        frame="UTM",
        detector="arca",
        plot_frame="galactic",
        detector_to="arca",
    )

    _ = kp.skymap_alert(
        file0=data_path("astro/antares_coordinate_systems_benchmark.csv"),
        frame="UTM",
        detector="antares",
        plot_frame="ParticleFrame",
        detector_to="antares",
    )

    _ = kp.skymap_alert(
        ra=80,
        dec=-20,
        obstime="2022-07-18T03:03:03",
        plot_frame="galactic",
        detector="dummy",
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


if __name__ == "__main__":
    main()
