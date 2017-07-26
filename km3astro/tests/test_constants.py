from unittest import TestCase

from km3astro.constants import (
    orca_easting, orca_northing, orca_utm_zone_letter, orca_utm_zone_number,
    orca_latitude, orca_longitude,
    arca_easting, arca_northing, arca_utm_zone_letter, arca_utm_zone_number,
    arca_latitude, arca_longitude,
)

# taken from loi
arca_latitude_naive = 36 + (16 / 60)  # degree
arca_longitude_naive = 16 + (6 / 60)  # degree
arca_lat_delta = 0.1
arca_lon_delta = 0.2

# taken from antares (should be the same)
orca_utm_zone_number_ref = 32
orca_utm_zone_letter_ref = 'T'

antares_easting = 268221.6
antares_northing = 4742381.9

orca_antares_delta_east = 12000    # m
orca_antares_delta_north = 600    # m


class TestUTM(TestCase):
    def test_ref_vs_computed(self):
        self.assertEqual(orca_utm_zone_number, orca_utm_zone_number_ref)
        self.assertEqual(orca_utm_zone_letter, orca_utm_zone_letter_ref)
        self.assertAlmostEqual(orca_easting, antares_easting,
                               delta=orca_antares_delta_east)      # 10 km
        self.assertAlmostEqual(orca_northing, antares_northing,
                               delta=orca_antares_delta_north)
        self.assertAlmostEqual(arca_longitude_naive, arca_longitude,
                               delta=arca_lon_delta)
        self.assertAlmostEqual(arca_latitude_naive, arca_latitude,
                               delta=arca_lat_delta)
