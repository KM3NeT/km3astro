"""Constants, like geographical positions."""

import utm

HEMISPHERE = 'north'
DATUM = 'WGS84'

orca_latitude = 42 + (48 / 60)  # degree
orca_longitude = 6 + (2 / 60)  # degree
orca_height = -2450  # m

# no longer needed (computet from latlon)
#orca_utm_zone_number = 32
#orca_utm_zone_letter = 'T'
orca_utm_x, orca_utm_y, orca_utm_zone_number, orca_utm_zone_letter = \
    utm.from_latlon(orca_latitude, orca_longitude)
orca_utm_zone = '{num}{let}'.format(
    num=orca_utm_zone_number, let=orca_utm_zone_letter)

arca_latitude = 36 + (16 / 60)  # degree
arca_longitude = 10 + (6 / 60)  # degree
arca_height = -3500  # m

# no longer needed (computet from latlon)
#arca_utm_zone_number = 33
#arca_utm_zone_letter = 'S'
arca_utm_x, arca_utm_y, arca_utm_zone_number, arca_utm_zone_letter = \
    utm.from_latlon(arca_latitude, arca_longitude)
arca_utm_zone = '{num}{let}'.farmat(
    num=arca_utm_zone_number, let=arca_utm_zone_letter)
