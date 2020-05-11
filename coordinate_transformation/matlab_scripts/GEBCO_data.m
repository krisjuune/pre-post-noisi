%Script to read the bathymetry .nc data set



%Import data from the GEBCO_2019 netCDF file
%The file has three variables lon(86400, double, degrees_east, 'X'), 
%lat(43200, double, degrees_north, 'Y'), elevation(86400x43200, single, 
%height above ref ellipsoid, i.e. elevation rel to sea level, m)

%To see the information about the variables and dimensions of the data, use
%ncdisp('GEBCO_2019.nc')

% raw_lon = ncread('GEBCO_2019.nc','lon'); 
% raw_lat = ncread('GEBCO_2019.nc','lat');
% raw_elevation = ncread('GEBCO_2019.nc','elevation');

% %Save raw lat, lon, and elevation as text files
% fid = fopen('raw_lon.txt','w');
%     fprintf(fid,'%6.2f  %12.8f\n',raw_lon');
%     fclose(fid);
%     
% fid = fopen('raw_lat.txt','w');
%     fprintf(fid,'%6.2f  %12.8f\n',raw_lat');
%     fclose(fid);
%     
% fid = fopen('raw_elevation.txt','w');
%     fprintf(fid,'%6.2f  %12.8f\n',raw_elevation');
%     fclose(fid);

%%
%Select depths for the domain specified
%Dimensions of the domain (~500x500km): 35,5...41,4N -22...-14,5E

% raw_lon = read('raw_lon.txt'); 
% raw_lat = read('raw_lat.txt');
% raw_elevation = read('raw_elevation.txt');

%Find and store indices of the domain longitudes
indx_lon = find(raw_lon<=-14.5 & raw_lon>=-22);
lon = raw_lon(indx_lon); 

%Find and store indices of the domain latitudes
indx_lat = find(raw_lat<=41.4 & raw_lat>=35.5); 
lat = raw_lat(indx_lat); 

%Create an array with the domain elevation data
bathy = raw_elevation(indx_lon,indx_lat); %add relevant elevation values

%Save truncated topography files
fid = fopen('lon_Prt.txt','w');
    fprintf(fid,'%6.2f  %12.8f\n',lon);
    fclose(fid);
    
fid = fopen('lat_Prt.txt','w');
    fprintf(fid,'%6.2f  %12.8f\n',lat);
    fclose(fid);
    
fid = fopen('bathy_Prt.txt','w');
    fprintf(fid,'%6.2f  %12.8f\n',bathy);
    fclose(fid);

%%
%Plot the bathymetry
[Lat, Lon] = meshgrid(lat, lon);
h = surf(Lon, Lat, bathy);
set(h,'LineStyle','none')
fontsize = 14; 
title('Bathymetry 35.5-41.4N 14.5-22W', 'Fontsize', fontsize)
xlabel('Longitude E', 'Fontsize', fontsize)
ylabel('Latitude N', 'Fontsize', fontsize)
zlabel('Depth (m)', 'Fontsize', fontsize)

