
# Array description

L-shaped array with 3 receivers in each direction.
Length of one degree of latitude is roughly 110 km,
whereas at 89.8N length of one degree of longitude
is about 0.39 km. This means that the longitude array
needs to be much more spaced out (in the stations file,
where the locations are given in degrees).

To space Lon receivers 12 km from one another, they
must be placed with about 30 degrees between them.
This is for a domain radius of roughly 27 km.

For a domain radius of roughly 1 degree, that is
110 km, the length of one degree of latitude is
about 112 km.The length of one degree of longitude
is 1.56 km. The Lon array below has 5 receivers placed
at 23.4 km (15 degree) intervals.

To check the length of one degree of longitude quickly
without having to run own code, use
[this online calculator](http://www.csgnetwork.com/degreelenllavcalc.html)
in case domain size has been changed and wish to expand
the Lat receivers to cover a larger area.

## Receiver locations for r = 27km

|Receiver name|Latitude   |Longitude   |
|-------------|-----------|------------|
|LAT0         |90.        |0.          |
|LAT1         |89.95      |0.          |
|LAT2         |89.9       |0.          |
|LAT3         |89.85      |0.          |
|LAT4         |89.8       |0.          |
|LON1         |89.8       |15          |
|LON2         |89.8       |30          |
|LON3         |89.8       |45          |
|LON4         |89.8       |60          |

## Receiver locations for r = 110km

|Receiver name|Latitude   |Longitude   |
|-------------|-----------|------------|
|LAT0         |90.        |0.          |
|LAT1         |89.8       |0.          |
|LAT2         |89.6       |0.          |
|LAT3         |89.4       |0.          |
|LAT4         |89.2       |0.          |
|LON1         |89.2       |15          |
|LON2         |89.2       |30          |
|LON3         |89.2       |45          |
|LON4         |89.2       |60          |
