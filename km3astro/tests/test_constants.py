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

# taken from antares (should be the same)
orca_utm_zone_number_ref = 32
orca_utm_zone_letter_ref = 'T'


class TestUTM(TestCase):
    def test_ref_vs_computed(self):
        self.assertEqual(orca_utm_zone_number, orca_utm_zone_number_ref)
        self.assertEqual(orca_utm_zone_letter, orca_utm_zone_letter_ref)
        self.assertAlmostEqual(arca_longitude_naive, arca_longitude, delta=0.2)
        self.assertAlmostEqual(arca_latitude_naive, arca_latitude, delta=0.1)
