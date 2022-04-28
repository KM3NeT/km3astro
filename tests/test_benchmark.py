import numpy as np
import pandas as pd

from km3net_testdata import data_path
from km3astro.coord import *


def test_Skycoord_separation(SC_true, SC_check):

    sep = SC_true.separation(SC_check)
    sep_deg = sep.deg
    return sep_deg


def test_benchmark_conversion(table_true, table_check):

    data_ = [table_true["SkyCoord_base"], table_check["SkyCoord_new"]]
    table_ = pd.concat(data_, axis=1)

    sep_table = table_.apply(
        lambda x: test_Skycoord_separation(x.SkyCoord_base, x.SkyCoord_new),
        axis=1,
        result_type="expand",
    )

    return sep_table


def test_angle_separation(
    file0="", detector_="antares", detector_to_="antares", to_print=False
):

    print("Starting Angle Separation Test...")
    if file0 == "":
        file0 = data_path("astro/antares_coordinate_systems_benchmark.csv")

    print("Read " + file0 + " file. \n")

    table_read = reader_from_file(file0)
    print(table_read)

    angle_treshold = 0.02
    # Calculate value for all frame from initial benchmark value

    if set(["phi", "theta"]).issubset(table_read.columns):

        table_loc_to_utm = transform_to_new_frame(
            table_read, "ParticleFrame", "UTM", detector_, detector_to_
        )
        table_loc_to_eq = transform_to_new_frame(
            table_read, "ParticleFrame", "equatorial", detector_, detector_to_
        )
        table_loc_to_gal = transform_to_new_frame(
            table_read, "ParticleFrame", "galactic", detector_, detector_to_
        )

    if set(["azimuth", "zenith"]).issubset(table_read.columns):

        table_utm_to_loc = transform_to_new_frame(
            table_read, "UTM", "ParticleFrame", detector_, detector_to_
        )
        table_utm_to_eq = transform_to_new_frame(
            table_read, "UTM", "equatorial", detector_, detector_to_
        )
        table_utm_to_gal = transform_to_new_frame(
            table_read, "UTM", "galactic", detector_, detector_to_
        )

    if set(["RA-J2000", "DEC-J2000"]).issubset(table_read.columns):

        table_eq_to_utm = transform_to_new_frame(
            table_read, "equatorial", "UTM", detector_, detector_to_
        )
        table_eq_to_loc = transform_to_new_frame(
            table_read, "equatorial", "ParticleFrame", detector_, detector_to_
        )
        table_eq_to_gal = transform_to_new_frame(
            table_read, "equatorial", "galactic", detector_, detector_to_
        )

    if set(["gal_lon", "gal_lat"]).issubset(table_read.columns):

        table_gal_to_loc = transform_to_new_frame(
            table_read, "galactic", "ParticleFrame", detector_, detector_to_
        )
        table_gal_to_eq = transform_to_new_frame(
            table_read, "galactic", "equatorial", detector_, detector_to_
        )
        table_gal_to_utm = transform_to_new_frame(
            table_read, "galactic", "UTM", detector_, detector_to_
        )

    # Start the Angle separation calculation

    # Checking Phi/Theta
    if set(["phi", "theta"]).issubset(table_read.columns):

        if set(["azimuth", "zenith"]).issubset(table_read.columns):
            sep_utm_to_loc = test_benchmark_conversion(
                table_loc_to_utm, table_utm_to_loc
            )

            mean_ = sep_utm_to_loc.mean()
            min_ = sep_utm_to_loc.min()
            max_ = sep_utm_to_loc.max()

            if to_print == True:
                print("UTM TO LOC")
                print(sep_utm_to_loc)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise Exception(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

        if set(["RA-J2000", "DEC-J2000"]).issubset(table_read.columns):
            sep_eq_to_loc = test_benchmark_conversion(table_loc_to_eq, table_eq_to_loc)

            mean_ = sep_eq_to_loc.mean()
            min_ = sep_eq_to_loc.min()
            max_ = sep_eq_to_loc.max()

            if to_print == True:
                print("EQ TO LOC")
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )
                print(sep_eq_to_loc)

            if max_ > angle_treshold:
                raise Exception(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

        if set(["gal_lon", "gal_lat"]).issubset(table_read.columns):

            sep_gal_to_loc = test_benchmark_conversion(
                table_loc_to_gal, table_gal_to_loc
            )

            mean_ = sep_gal_to_loc.mean()
            min_ = sep_gal_to_loc.min()
            max_ = sep_gal_to_loc.max()

            if to_print == True:
                print("GAL TO LOC")
                print(sep_gal_to_loc)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise Exception(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

    # Checking Az/Ze
    if set(["azimuth", "zenith"]).issubset(table_read.columns):

        if set(["phi", "theta"]).issubset(table_read.columns):
            sep_loc_to_utm = test_benchmark_conversion(
                table_utm_to_loc, table_loc_to_utm
            )

            mean_ = sep_loc_to_utm.mean()
            min_ = sep_loc_to_utm.min()
            max_ = sep_loc_to_utm.max()

            if to_print == True:
                print("LOC TO UTM")
                print(sep_loc_to_utm)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise Exception(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

        if set(["RA-J2000", "DEC-J2000"]).issubset(table_read.columns):

            sep_eq_to_utm = test_benchmark_conversion(table_utm_to_eq, table_eq_to_utm)

            mean_ = sep_eq_to_utm.mean()
            min_ = sep_eq_to_utm.min()
            max_ = sep_eq_to_utm.max()

            if to_print == True:
                print("EQ TO UTM")
                print(sep_eq_to_utm)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise Exception(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

        if set(["gal_lon", "gal_lat"]).issubset(table_read.columns):
            sep_gal_to_utm = test_benchmark_conversion(
                table_utm_to_gal, table_gal_to_utm
            )

            mean_ = sep_gal_to_utm.mean()
            min_ = sep_gal_to_utm.min()
            max_ = sep_gal_to_utm.max()

            if to_print == True:
                print("GAL TO UTM")
                print(sep_gal_to_utm)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise Exception(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

    # Checking RA/DEC
    if set(["RA-J2000", "DEC-J2000"]).issubset(table_read.columns):

        if set(["phi", "theta"]).issubset(table_read.columns):
            sep_loc_to_eq = test_benchmark_conversion(table_eq_to_loc, table_loc_to_eq)

            mean_ = sep_loc_to_eq.mean()
            min_ = sep_loc_to_eq.min()
            max_ = sep_loc_to_eq.max()

            if to_print == True:
                print("LOC TO EQ")
                print(sep_loc_to_eq)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise Exception(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

        if set(["azimuth", "zenith"]).issubset(table_read.columns):
            sep_utm_to_eq = test_benchmark_conversion(table_eq_to_utm, table_utm_to_eq)

            mean_ = sep_utm_to_eq.mean()
            min_ = sep_utm_to_eq.min()
            max_ = sep_utm_to_eq.max()

            if to_print == True:
                print("UTM TO EQ")
                print(sep_utm_to_eq)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise Exception(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

        if set(["gal_lon", "gal_lat"]).issubset(table_read.columns):
            sep_gal_to_eq = test_benchmark_conversion(table_eq_to_gal, table_gal_to_eq)

            mean_ = sep_gal_to_eq.mean()
            min_ = sep_gal_to_eq.min()
            max_ = sep_gal_to_eq.max()

            if to_print == True:
                print("GAL TO EQ")
                print(sep_gal_to_eq)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise Exception(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

    # Checking GAL L/B
    if set(["gal_lon", "gal_lat"]).issubset(table_read.columns):

        if set(["phi", "theta"]).issubset(table_read.columns):
            sep_loc_to_gal = test_benchmark_conversion(
                table_gal_to_loc, table_loc_to_gal
            )

            mean_ = sep_loc_to_gal.mean()
            min_ = sep_loc_to_gal.min()
            max_ = sep_loc_to_gal.max()

            if to_print == True:
                print("LOC TO GAL")
                print(sep_loc_to_gal)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise ValueError(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

        if set(["azimuth", "zenith"]).issubset(table_read.columns):
            sep_utm_to_gal = test_benchmark_conversion(
                table_gal_to_utm, table_utm_to_gal
            )

            mean_ = sep_utm_to_gal.mean()
            min_ = sep_utm_to_gal.min()
            max_ = sep_utm_to_gal.max()

            if to_print == True:
                print("UTM TO GAL")
                print(sep_utm_to_gal)
                print(
                    "Mean = "
                    + str(mean_)
                    + " Max = "
                    + str(max_)
                    + " Min = "
                    + str(min_)
                )

            if max_ > angle_treshold:
                raise ValueError(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

        if set(["RA-J2000", "DEC-J2000"]).issubset(table_read.columns):
            sep_eq_to_gal = test_benchmark_conversion(table_gal_to_eq, table_eq_to_gal)

            mean_ = sep_eq_to_gal.mean()
            min_ = sep_eq_to_gal.min()
            max_ = sep_eq_to_gal.max()

            if to_print == True:
                print("EQ TO GAL")
                print(sep_eq_to_gal)
                print(f"Mean = {mean_} Max = {max_} Min = {min_}")

            if max_ > angle_treshold:
                raise ValueError(
                    "Error: Maximum angle separation = "
                    + str(max_)
                    + " > "
                    + str(angle_treshold)
                )

    print("... End of Angle Separation Test")
    return 0


def main():

    # testing_angle_separation of :

    # coordinate system benchmark

    test_angle_separation(to_print=False)
    test_angle_separation(
        data_path("astro/ARCA_coordinate_systems_benchmark.csv"), "arca", "arca"
    )
    test_angle_separation(
        data_path("astro/ORCA_coordinate_systems_benchmark.csv"), "orca", "orca"
    )

    # astro object benchmark
    test_angle_separation(
        data_path("astro/antares_astro_objects_benchmark.csv"), "antares", "antares"
    )
    test_angle_separation(
        data_path("astro/ARCA_astro_objects_benchmark.csv"), "arca", "arca"
    )
    test_angle_separation(
        data_path("astro/ORCA_astro_objects_benchmark.csv"), "orca", "orca"
    )

    # moon sun postion benchmark
    test_angle_separation(
        data_path("astro/antares_moon_sun_position_benchmark.csv"), "antares", "antares"
    )

    # need new benchmark table:
    # test_angle_separation(data_path("astro/ARCA_moon_sun_position_benchmark.csv"), "arca", "arca")
    # test_angle_separation(data_path("astro/ORCA_moon_sun_position_benchmark.csv"), "orca", "orca")


if __name__ == "__main__":
    main()
