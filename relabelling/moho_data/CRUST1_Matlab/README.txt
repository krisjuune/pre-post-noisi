README.txt

The CRUST1.0 database has been reorganized into a matlab structure called C1 that has been saved
in a matlab data file called CRUST1.mat.   This file and two low-level access functions called
findCrust1cell.m and getCrust1.m comprise the contents of the folder CRUST1 which constitutes a 
matlab toolbox. This folder should be put with the user’s other toolboxes, and it should be added
to the path using pathtool.m or addpath.m

An example script called queryCRUST1.m utilizes this toolbox. Normally this script would be placed
in a project folder rather than in the toolbox. When you run script queryCRUST1, it will to prompt
you to input the latitude (lat) and longitude (lon) of a point or station of interest, and then it
write out the vertical profile for with each CRUST1.0 cell or tile associated with the station.

Please send bug reports to Michael Bevis at mbevis@gmail.com

