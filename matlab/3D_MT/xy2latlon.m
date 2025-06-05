function [lat,lon] = xy2latlon(x,y,lat0,lon0,units,projlat)

% [lat,lon] = xy2latlon(x,y,lat0,lon0,units,projlat)
% x = R (th-th0) & y = R cos(th) (ph-ph0)
% Earth's equator radius 6378.137 km makes it
% compatible with m_map (m_lldist etc)
%
% optionally, specify the "projection latitude"
% to use one (average) latitude for all mappings
% (for example, projlat = 42 for Earthscope models) 

% avglon = repmat(mean(lon),1,length(lat));
% avglat = repmat(mean(lat),1,length(lon));
% m_proj('equidistant');
% dx = m_lldist(avglon,lat);
% dy = m_lldist(lon,avglat);
% km = [sum(dx)/2 sum(dy)/2 0];

if nargin > 4
    if strcmp(units,'km')
        x = x*1000;
        y = y*1000;
    end
end

th0 = repmat(lat0,size(x)) * (pi/180);
ph0 = repmat(lon0,size(y)) * (pi/180);
R = 6378137; %6731500;

th  = (x/R) + th0;
if nargin > 5
    ph = (y./(R*cos(projlat*pi/180))) + ph0;
else
    ph = (y./(R*cos(th))) + ph0; 
end
    
lat = th * (180/pi);
lon = ph * (180/pi);


% lat = lat0 + km2deg(x/1000);
% lon = lon0 + km2deg(y/1000).*(cosd(lat))^2;

