# UtmLatLongConverter

A UTM -> Lat/Long (or vice versa) converter.

## How to Use

Call setDatum with the index of the datum to be used and then call the various
conversion functions latLngToUtm(), utmToLatLng() or natoToLatLng(), utmToNato(),
natoToUtm() to convert between the various coordinate systems, hopefully accurately!

## Caveat

**NOTE**: no attempt is made to compensate for the irregular grid in the area around
the southern coast of Norway and Svalbard (zones 32V and 31X, 33X, 35X and 37X)
because of this results returned for NATO coordinates for lat/lng or
UTM values located in these regions will not be correct.

## Acknowledgments

This converter has been translated into python from the original javascript
file develop at [Montana University](http://www.rcn.montana.edu/Resources/Converter.aspx), 
which in turn was adapted from the script used at http://www.uwgb.edu/dutchs/UsefulData/ConvertUTMNoOZ.HTM
