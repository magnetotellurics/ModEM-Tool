function [x,y,z] = latlon2cart(lat,long,varargin)
% Convert lat, long, z to x, y, z Cartesian centered at the centre of Earth
%
% INPUT
% lat   - latitude
% long  - longitude
% z     - depth
%
% OUTPUT
% x, y, z coordinates with origo at the centre of Earth

lat = deg2rad(lat);
long = deg2rad(lat);
R = 6378.1; % radius of Earth
r = R;%-z;
x = r .* cos(lat) .* cos(long);
y = r .* cos(lat) .* sin(long);
z = r .* sin(lat);