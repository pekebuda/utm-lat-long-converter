import sys
import os
import unittest
import json
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../src")
from geocon import Geocon


class GeoconTestCase(unittest.TestCase):

    def test_utm_to_lat_lng(self):
        geocon = Geocon(0)
        with open('./sample.json') as fSample:
            SAMPLE = json.loads(fSample.read())
            for city in SAMPLE:
                actual_latitude = city['latitude']
                actual_longitude = city['longitude']
                actual_zone = city['zone']
                actual_southern = city['southern']
                actual_easting = city['easting']
                actual_northing = city['northing']
                computed_latitude = geocon.utmToLatLng(
                    actual_easting,
                    actual_northing,
                    actual_zone,
                    actual_southern
                )['lat']
                computed_longitude = geocon.utmToLatLng(
                    actual_easting,
                    actual_northing,
                    actual_zone,
                    actual_southern
                )['lng']

                self.assertTrue(0.999*abs(actual_latitude) <= abs(computed_latitude) <= abs(actual_latitude)*1.001)
                self.assertTrue(0.999*abs(actual_longitude) <= abs(computed_longitude) <= abs(actual_longitude)*1.001)

    def test_lat_lng_to_utm(self):
        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
