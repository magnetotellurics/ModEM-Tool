function [lat,lon] = xformToLL(x,y,Origins)
%  Usage:   [lat,lon] = xformToLL(x,y,Origins)

%   I am assuming x is east-west, y is North-south

degPerKm = 180/(6370.8*pi);
lat = (y-Origins.y0)*degPerKm+Origins.lat0;
%lon = (x-Origins.x0)*degPerKm/cos(Origins.lat0*pi/180)+Origins.lon0;
lon = (x-Origins.x0)*degPerKm./cos(lat*pi/180)+Origins.lon0;
