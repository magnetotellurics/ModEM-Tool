function [iSites,ix,iy] = findSites2D(siteLoc,x,y,rho,rhoMax,rhoMin)

% function iSites = findSites2D(siteLoc,x,y,rho,rhoMax,rhoMin)
%
% Finds sites in siteLoc size nSites*2 for which
% rhoMin <= rho <= rhoMax, where rho is 2D size nx*ny.
% Can be used to find sites in the ocean after regridding,
% by setting rhoMin = 0 and rhoMax = 0.35
% (usual resistivity of ocean water is 0.3 Ohm).

if nargin < 5
    rhoMax = 0.35;
end

if nargin < 6
    rhoMin = 0.;
end

iSites = [];
ix = [];
iy = [];
for i = 1:size(siteLoc,1)
   [xval,xind]=max(x(x<=siteLoc(i,1)));
   [yval,yind]=max(y(y<=siteLoc(i,2)));
   if rho(xind,yind) <= rhoMax && rho(xind,yind) >= rhoMin
       iSites = [iSites i];
       ix = [ix xind];
       iy = [iy yind];
   end
end