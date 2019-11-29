# topography-transformation
Transformation of data with geographic coordinates into Cartesian & cylindrical coordinates

# Overview

- [x] Transform to goecentric coordinates
- [x] Transform to spherical coordinate system
- [x] Transform to Cartesian
- [x] Transform to cylindrical (from spherical)
- [ ] Check if transofrmation valid - transform back and plot both before and 'after after'  
- [ ] Rotate domain to be centred about the N pole 

1. Obtain coordinates for spherical coordinate system
1. Transformation to Cartesian & cylindrical coordinates
1. Test the validity of transformation
   1. Plot data in geographic coordinates
   1. Plot data in spherical coordinates
   1. Plot data in Cartesian/cylindrical

## Step 1
Spherical polar coordinates (radius, colatitude, azimuth) for the bathymetry data obtained via transforming latitudes from geographic to geocentric coordinates, calculating the radius of the reference ellipsoid WGS84 at each of the latitudes, and then calculating the distance from the centre to each of the topography data points. Colatitude is found as 90 minus geocentric latitude. 
This gives the bathymetry (distance from reference ellipsoid) and geographic latitude-longitude as radius, colatitude and azimuth (=longitude). 

## Step 2
Transform the data from spherical polar to Cartesian and cylindrical coordinates, e.g. [Wiki_page](https://en.wikipedia.org/wiki/Spherical_coordinate_system). 
For the transformation operation, r, colatitude, and azimuth must all be arrays of identical shape. 
Rotate domain to be centred about the N pole, using the rotation matrix in 

## Step 3
Validify the transformation by plotting in: 
1. geographic coordinates
2. spherical coordinates
3. Cartesian/cylindrical coordinates

