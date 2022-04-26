import numpy as np
import pandas as pd

from km3net_testdata import data_path
from km3astro.coord import *


def test_coordinate_transformator(
    file0="",
    detector_="antares",
    detector_to_="antares",
    type_name="coordinate_systems",
):

    print("Starting Test Coordinate Transformator...")
    # take a benchmark and create new benchmark table

    if file0 == "":
        file0 = data_path("astro/antares_coordinate_systems_benchmark.csv")

    table_read = reader_from_file(file0)

    print("Read " + file0 + " file. \n")
    print(table_read)
    path = "/tmp/"
    g_comment = "# benchmark table from https://git.km3net.de/km3py/km3net-testdata/-/tree/master/km3net_testdata/data/astro\n# values were obtained with the corresponding benchmark "

    has_object = False
    date_ = table_read["date"]
    time_ = table_read["time"]
    if "object" in table_read.columns:
        object_ = table_read["object"]
        has_object = True

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

        phi_ = table_read["phi"]
        theta_ = table_read["theta"]

        az_, ze_ = zip(
            *table_loc_to_utm.apply(
                lambda x: get_az_zenith(x.SkyCoord_new, detector_to_), axis=1
            )
        )

        ra_, dec_ = zip(
            *table_loc_to_eq.apply(lambda x: get_ra_dec(x.SkyCoord_new), axis=1)
        )
        l_, b_ = zip(*table_loc_to_gal.apply(lambda x: get_l_b(x.SkyCoord_new), axis=1))

        data_ = {
            "date": date_,
            "time": time_,
            "phi": phi_,
            "theta": theta_,
            "azimuth": az_,
            "zenith": ze_,
            "RA-ICRS": ra_,
            "DEC-ICRS": dec_,
            "gal_lon": l_,
            "gal_lat": b_,
        }
        if has_object == True:
            data_ = {
                "date": date_,
                "time": time_,
                "object": object_,
                "phi": phi_,
                "theta": theta_,
                "azimuth": az_,
                "zenith": ze_,
                "RA-ICRS": ra_,
                "DEC-ICRS": dec_,
                "gal_lon": l_,
                "gal_lat": b_,
            }

        df_ = pd.DataFrame(data=data_)

        comment = (
            g_comment
            + "for 'phi' and 'theta'\n# other value are obtained with km3py_astro.\n# units: [deg]\n# theta, phi: track in horiz. UTM\n# zenith, azimuth: Source in horiz. UTM\n# DEC/RA: source in equatorial coos.\n# gal_lon, gal_lat: source in galactic coos.\n"
        )
        name = path + detector_ + "_" + type_name + "_benchmark_km3py_from_track.csv"
        f = open(name, "w")
        f.write(comment)

        df_.to_csv(f, index=False)
        f.close()

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

        az_ = table_read["azimuth"]
        ze_ = table_read["zenith"]

        phi_, theta_ = zip(
            *table_utm_to_loc.apply(
                lambda x: get_phi_theta(x.SkyCoord_new, detector_to_), axis=1
            )
        )

        ra_, dec_ = zip(
            *table_utm_to_eq.apply(lambda x: get_ra_dec(x.SkyCoord_new), axis=1)
        )
        l_, b_ = zip(*table_utm_to_gal.apply(lambda x: get_l_b(x.SkyCoord_new), axis=1))

        data_ = {
            "date": date_,
            "time": time_,
            "phi": phi_,
            "theta": theta_,
            "azimuth": az_,
            "zenith": ze_,
            "RA-ICRS": ra_,
            "DEC-ICRS": dec_,
            "gal_lon": l_,
            "gal_lat": b_,
        }
        if has_object == True:
            data_ = {
                "date": date_,
                "time": time_,
                "object": object_,
                "phi": phi_,
                "theta": theta_,
                "azimuth": az_,
                "zenith": ze_,
                "RA-ICRS": ra_,
                "DEC-ICRS": dec_,
                "gal_lon": l_,
                "gal_lat": b_,
            }

        df_ = pd.DataFrame(data=data_)

        comment = (
            g_comment
            + "for 'azimuth' and 'zenith'\n# other value are obtained with km3py_astro.\n# units: [deg]\n# theta, phi: track in horiz. UTM\n# zenith, azimuth: Source in horiz. UTM\n# DEC/RA: source in equatorial coos.\n# gal_lon, gal_lat: source in galactic coos.\n"
        )
        name = path + detector_ + "_" + type_name + "_benchmark_km3py_from_utm.csv"
        f = open(name, "w")
        f.write(comment)

        df_.to_csv(f, index=False)
        f.close()

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

        ra_ = table_read["RA-J2000"]
        dec_ = table_read["DEC-J2000"]

        phi_, theta_ = zip(
            *table_eq_to_loc.apply(
                lambda x: get_phi_theta(x.SkyCoord_new, detector_to_), axis=1
            )
        )

        az_, ze_ = zip(
            *table_eq_to_utm.apply(
                lambda x: get_az_zenith(x.SkyCoord_new, detector_to_), axis=1
            )
        )

        l_, b_ = zip(*table_eq_to_gal.apply(lambda x: get_l_b(x.SkyCoord_new), axis=1))

        data_ = {
            "date": date_,
            "time": time_,
            "phi": phi_,
            "theta": theta_,
            "azimuth": az_,
            "zenith": ze_,
            "RA-J2000": ra_,
            "DEC-J2000": dec_,
            "gal_lon": l_,
            "gal_lat": b_,
        }
        if has_object == True:
            data_ = {
                "date": date_,
                "time": time_,
                "object": object_,
                "phi": phi_,
                "theta": theta_,
                "azimuth": az_,
                "zenith": ze_,
                "RA-J2000": ra_,
                "DEC-J2000": dec_,
                "gal_lon": l_,
                "gal_lat": b_,
            }

        df_ = pd.DataFrame(data=data_)

        comment = (
            g_comment
            + "for 'RA-J2000' and 'RA-J2000'\n# other value are obtained with km3py_astro.\n# units: [deg]\n# theta, phi: track in horiz. UTM\n# zenith, azimuth: Source in horiz. UTM\n# DEC/RA: source in equatorial coos.\n# gal_lon, gal_lat: source in galactic coos.\n"
        )
        name = (
            path + detector_ + "_" + type_name + "_benchmark_km3py_from_equatorial.csv"
        )
        f = open(name, "w")
        f.write(comment)

        df_.to_csv(f, index=False)
        f.close()

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

        l_ = table_read["gal_lon"]
        b_ = table_read["gal_lat"]

        phi_, theta_ = zip(
            *table_utm_to_loc.apply(
                lambda x: get_phi_theta(x.SkyCoord_new, detector_to_), axis=1
            )
        )
        az_, ze_ = zip(
            *table_eq_to_utm.apply(
                lambda x: get_az_zenith(x.SkyCoord_new, detector_to_), axis=1
            )
        )

        ra_, dec_ = zip(
            *table_utm_to_eq.apply(lambda x: get_ra_dec(x.SkyCoord_new), axis=1)
        )

        data_ = {
            "date": date_,
            "time": time_,
            "phi": phi_,
            "theta": theta_,
            "azimuth": az_,
            "zenith": ze_,
            "RA-ICRS": ra_,
            "DEC-ICRS": dec_,
            "gal_lon": l_,
            "gal_lat": b_,
        }
        if has_object == True:
            data_ = {
                "date": date_,
                "time": time_,
                "object": object_,
                "phi": phi_,
                "theta": theta_,
                "azimuth": az_,
                "zenith": ze_,
                "RA-ICRS": ra_,
                "DEC-ICRS": dec_,
                "gal_lon": l_,
                "gal_lat": b_,
            }

        df_ = pd.DataFrame(data=data_)

        comment = (
            g_comment
            + "for 'gal_lon' and 'gal_lat'\n# other value are obtained with km3py_astro.\n# units: [deg]\n# theta, phi: track in horiz. UTM\n# zenith, azimuth: Source in horiz. UTM\n# DEC/RA: source in equatorial coos.\n# gal_lon, gal_lat: source in galactic coos.\n"
        )
        name = path + detector_ + "_" + type_name + "_benchmark_km3py_from_galactic.csv"
        f = open(name, "w")
        f.write(comment)

        df_.to_csv(f, index=False)
        f.close()

    print("... End of Test Coordinate Transformator")

    return 0


def main():

    test_coordinate_transformator()
    # test_coordinate_transformator(data_path("astro/ARCA_coordinate_systems_benchmark.csv"), "arca", "arca")
    # test_coordinate_transformator(data_path("astro/ORCA_coordinate_systems_benchmark.csv"), "orca", "orca")
    # test_coordinate_transformator(data_path("astro/ARCA_astro_objects_benchmark.csv"), "arca", "arca", type_name = "astro_objects")
    # test_coordinate_transformator(data_path("astro/ORCA_moon_sun_position_benchmark.csv"), "orca", "orca", type_name = "moon_sun_position")


if __name__ == "__main__":
    main()
