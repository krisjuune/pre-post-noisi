%Script to transform bathymetry data from spherical to Cartesian
%coordinates using the Jacobian. 

%Load the data
lon = load('lon14.5_22W.txt'); %degrees E
lat = load('lat35.5_41.4N.txt'); %degrees N
bathy = load('bathy22W41.4N.txt'); %m, elevation below sealevel -ve

%% Format data

%Elevation must be expressed as a radius
r_ref = 6370.074; %km, reference radius at 38N
    %This assumes Earth's radius is constant over this domain which is not
    %strictly true, could use the equation in method notes to calculate actual
    %radius for each of the latitudes in the domain. 
bathy_r = bathy./1000 + r_ref; %km, distance from Earth's centre


%Azimuth & elevation (lon, lat) must be in radians
lon_rad = lon.*(pi/180); %rad, counterclockwise anle in x-y from +ve x
lat_rad = lat.*(pi/180); %rad, elevation angle from x-y plane
%Need to mesh the lat, lon for compatibility with bathy size
[Lat_rad, Lon_rad] = meshgrid(lat_rad,lon_rad); 

%Transform (0,0) to OBS coordinates (rotate receiver to N pole)????

%% Transformation into Cartesian

% [x,y,z] = sph2cart(Lon_rad, Lat_rad, bathy_r); 
% 
% h = surf(x, y, z);
% set(h,'LineStyle','none')
% fontsize = 14; 
% scheisse

[x,y,z] = axesm([35.5 41.4], [-22 -14.5], [-7000 0]); 
[x,y,z] = mfwdtran(Lat, Lon, bathy); 
bathy_doubel = double(bathy); 

geoshow(lat, lon, bathy, 'surface')

