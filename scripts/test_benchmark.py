import numpy as np
import pandas as pd

from km3net_testdata import data_path

import km3astro.coord as kc 
import km3astro.testing_tools as ktt


def main():

    # testing_angle_separation of :

    # coordinate system benchmark

    ktt.test_angle_separation(
        data_path("astro/antares_coordinate_systems_benchmark.csv"),"antares","antares"
    )
    ktt.test_angle_separation(
        data_path("astro/ARCA_coordinate_systems_benchmark.csv"), "arca", "arca"
    )
    ktt.test_angle_separation(
        data_path("astro/ORCA_coordinate_systems_benchmark.csv"), "orca", "orca"
    )

    # astro object benchmark
    ktt.test_angle_separation(
        data_path("astro/antares_astro_objects_benchmark.csv"), "antares", "antares"
    )
    ktt.test_angle_separation(
        data_path("astro/ARCA_astro_objects_benchmark.csv"), "arca", "arca"
    )
    ktt.test_angle_separation(
        data_path("astro/ORCA_astro_objects_benchmark.csv"), "orca", "orca"
    )

    # moon sun postion benchmark
    ktt.test_angle_separation(
        data_path("astro/antares_moon_sun_position_benchmark.csv"), "antares", "antares"
    )

    # need new benchmark table:
    # test_angle_separation(data_path("astro/ARCA_moon_sun_position_benchmark.csv"), "arca", "arca")
    # test_angle_separation(data_path("astro/ORCA_moon_sun_position_benchmark.csv"), "orca", "orca")


if __name__ == "__main__":
    main()
