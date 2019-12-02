%Script to test if .nc file created in Python for topography works 
%29.11.2019

ncdisp('topography_coord.nc')

%Output: 
% Source:
%            /Users/kristiinajoon/Desktop/4th_year/4th_year_project/Project/topography_coord.nc
% Format:
%            netcdf4
% Global Attributes:
%            description = 'Togography data in sph, Cart, and cyl coordinates'
%            history     = 'Created 29/11/19'
% Dimensions:
%            colat = 1416
%            lon   = 1800
%            x     = 1416
%            y     = 1416
%            z     = 1416
%            s     = 1416
%            phi   = 1800
%            z_cyn = 1416
% Variables:
%     radial_distance    
%            Size:       1800x1416
%            Dimensions: lon,colat
%            Datatype:   single
%            Attributes:
%                        units = 'km'
%     colatitude         
%            Size:       1416x1
%            Dimensions: colat
%            Datatype:   single
%            Attributes:
%                        units = 'degrees north'
%     longitude          
%            Size:       1800x1
%            Dimensions: lon
%            Datatype:   single
%            Attributes:
%                        units = 'degrees east'
%     x_value            
%            Size:       1416x1
%            Dimensions: x
%            Datatype:   single
%            Attributes:
%                        units = 'km'
%     y_value            
%            Size:       1416x1
%            Dimensions: y
%            Datatype:   single
%            Attributes:
%                        units = 'km'
%     z_value            
%            Size:       1416x1
%            Dimensions: z
%            Datatype:   single
%            Attributes:
%                        units = 'km'
%     cyl_radial_distance
%            Size:       1416x1
%            Dimensions: s
%            Datatype:   single
%            Attributes:
%                        units = 'km'
%     azimuth            
%            Size:       1800x1
%            Dimensions: phi
%            Datatype:   single
%            Attributes:
%                        units = 'degrees east'
%     height             
%            Size:       1416x1
%            Dimensions: z_cyn
%            Datatype:   single
%            Attributes:
%                        units = 'km'
