function [b,p,s,r] = getCrust1(i,j,C1)
%getCrust1 get vertical profile from the (i,j) cell or tile of CRUST1.0
%
% See the CRUST1.0 Readme file on the CRUST1.0 website at URL
% https://igppweb.ucsd.edu/~gabi/crust1.html#download for more details
% about the CRUST1.0 model
%
% USAGE:
%           [b,p,s,r] = getCrust1(i,j);
%           [b,p,s,r] = getCrust1(i,j,C1);
%
% INPUT:
%     i    scalar. The integer index associated with cell latitude
%                  1 <= i <= 180
%     j    scalar. The integer index associated with cell longitude
%                  1 <= j <= 360
%   C1  a structure containing the CRUST1.0 grid (this can be found in the
%         matlab data file CRUST1.mat). If getCrust1 is to be called many
%         times, providing C1 directly speeds things up since it is not
%         then necessary to load CRUST1.mat many times. If C1 is not 
%         provided via the input list, then it we be obtained by loading
%         Crust1.mat
% OUTPUT:
%    b   boundary topography vector with 9 elements 
%    p   Vp vector with 9 elements giving the layer P-wave velocties
%          p(9) is VPn at the top of the mantle
%    s   Vs vector with 9 elements giving the layer S-wave velocties
%          s(9) is VSn at the top of the mantle
%    r    vector with 9 elements giving the layer densities
%
% See also function findCrust1cell.m and dataset CRUST1.mat

%  version 1.0             Michael Bevis        20 June 2017

if nargin<2 || nargin>3
    error('getCrust1 takes 2 or 3 input arguments')
end
if nargin==2
    load CRUST1
end
if numel(i)~=1 || numel(j)~=1
    error('input arguments i nd j must be scalar')
end
if rem(i,1)~=0 || rem(j,1)~=0
    error('lat and lon values must be half degrees')
end
if i<1 || i>180
    error('i has an illegal values')
end
if j<1 || j>360
    error('j has an illegal value')
end
b=squeeze(C1.B(i,j,:));
p=squeeze(C1.P(i,j,:));
s=squeeze(C1.S(i,j,:));
r=squeeze(C1.R(i,j,:));

