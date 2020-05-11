# Which simulations

Which simulation should be done for bechmark and how
to display the results of each simulation.

- [Which simulations](#which-simulations)
  - [General plots](#general-plots)
  - [Test curvature](#test-curvature)
    - [Which plots](#which-plots)
  - [Test effect of bathymetry](#test-effect-of-bathymetry)
  - [Test sensitivity to bathymetry](#test-sensitivity-to-bathymetry)
  - [Test sensitivity to mesh resolution](#test-sensitivity-to-mesh-resolution)

## General plots

1. domain bathymetry
  From 35.5 - 39.5 N, 14 - 19 E. Lots of bathymetric variation in the SE part of the region - hence the cross-shaped array was chosen to cross that and contrast the SW part where the bathymetry smoothly deepens.
  ![bathymetry](thesis_plots/../../../thesis_plots/Feb_dom_bathymetry_more.png)
1. stations and bathymetry
  ![stations](thesis_plots/../../../thesis_plots/stations.png)
1. bathymetry and wavelengths
  To compare the horizontal lengthscale of bathymetric variations to the seismic wavelengths - should show that at longer wavelengths the sensitivity to bathymetry is lost.

## Test curvature

This involves three simulations:

- Cartesian with all interfaces being flat
- Spherical with all interfaces curved (calculated for a spherical Earth with average radius of the domain)
- Geographic with all interfaces curved as defined by the WGS84 reference ellipsoid

### Which plots

Mainly I'm thinking just a simple plot with the three
results overlain at three different distances (one at
source, one relatively close, and one as far as
possible).
To display the plot, save it as a png. pdf files are not loaded in Markdown preview.
![curvature](thesis_plots/../../../thesis_plots/which_curvature.pdf)

## Test effect of bathymetry

Plot with vs without bathymetry as a move out plot,
more information about this in notebook Feb notes

Bathymetry of the domain with the different wave lengths for different source periods:
at 36.5N the length of one degree of longitude is about 90 km.

| Source period | Wavelength in water | Wavelength in solid | Rel to 1degree |
| ------------- |:-------------------:|:-------------------:|:--------------:|
| 2.5 sec       | 2.5 km | 7.5 km | 1/36(0.03); 1/12 |
| 5 sec         | 5 km   | 15 km  | 1/18(0.06); 1/6  |
| 10 sec        | 10 km  | 30 km  | 1/9; 1/3         |
| 20 sec        | 20 km  | 60 km  | 0.22; 2/3        |

![bathy_moveout](../../thesis_plots/bathy_moveout_5s.png)

## Test sensitivity to bathymetry

Test the effect of different source periods, and see
at which wavelengths the solution with bathymetry
starts to converge with the solution without
bathymetry, i.e. flat ocean bottom.

Plot again as a move out with the different wavelengths
overlain atop one another.

## Test sensitivity to mesh resolution

Run same as in last one but on a finer mesh and see
the differences of results.
