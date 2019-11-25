function [i,j]= findCrust1cell(lat,lon)
%findCrust1cell find the CRUST1.0 cell(s) containing the point (lat,lon)
% CRUST1.0 was constructed one 1 x 1 degree cell at a time. If a station
% or point (lat,lon) lies within a cell, the cell indices (i,j) of that 
% cell are returned. If the station sits on the border of two cells, then 
% the indicess of both cells are returned. If the station sites
% on a vertex shared by four cells, then the indices of all
% four cells are returned.
%
% USAGE:        [i,j]= findCrust1cell(lat,lon);
%
% INPUT:
%   lat   scalar containing the latitude in degrees (-90>=lat<=+90)
%   lon   scalar containing the longitude in degrees (-180>=lat<=+180)
% 
% OUTPUT:
%    i    row index in the various matrices of the C1 structure
%    j    column index in the various matrices of the C1 structure
%
% Notes: 
% (1) In this matlab implementation of CRUST1.0, i is the latitude- 
% or y- related index, and j is the longitude- or x-related index.
%
% (2) The user can count the number of cells associated with (lat,lon)
% using the statement  ncell=length(i)
%
% See also functions getCrust1.m and data file Crust1.mat

%  Version 1.0             Michael Bevis            20 June 2017
if nargin~=2
    error('findnearestCrust1gridpoints takes 2 input arguments')
end
if numel(lat)~=1 || numel(lon)~=1
    error('input arguments lat and lon must be scalar')
end
if lon< -180
    lon=lon+360;
elseif lon> 180
    lon=lon-360;
end

if lat==90
    lat=89.999; 
elseif lat==-90
    lat=-89.999;
end
flat=fix(lat);
if rem(lat,1)~=0
    if lat>0
        i=flat+91;
    else
        i=flat+90;
    end
else
    i=[flat+90 flat+91];
end

flon=fix(lon);
if rem(lon,1)~=0
    if lon>0
        j=flon+181;
    else
        j=flon+180;
    end
else
    if flon== 180 || flon== -180
        j=[1 360];
    else
        j=[flon+180 flon+181];
    end
end
i=i(:); j=j(:);
ni=length(i);
nj=length(j);
if ni>1 || nj>1
    if ni==1 & nj==2
        I=[i i]'; J=j;
    elseif ni==2 & nj==1
        I=i; J=[j j]';
    else
        I=[i(1) i(1) i(2) i(2)]';
        J=[j(1) j(2) j(1) j(2)]';
    end
    i=I(:);
    j=J(:);
end

