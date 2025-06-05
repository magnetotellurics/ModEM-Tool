function [x,y] = latlon2xy(lat,lon,lat0,lon0,units,projlat)

% [x,y] = latlon2xy(lat,lon,lat0,lon0,units,projlat)
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

% lat_km = distance(lat0,lon0,lat,lon0,[6371 0]);
% lon_km = distance(lat0,lon0,lat0,lon,[6371 0]);
% x = 1000*lat_km;
% y = 1000*lon_km;

if isscalar(lat)
    lat = lat * ones(size(lon));
end

if isscalar(lon)
    lon = lon * ones(size(lat));
end

th  = lat * (pi/180);
th0 = repmat(lat0,size(lat)) * (pi/180);
ph  = lon * (pi/180);
ph0 = repmat(lon0,size(lon)) * (pi/180);
R = 6378137; %6731500;

x = R * (th - th0);
if nargin > 5
    y = R * cos(projlat*pi/180) .* (ph - ph0);
else
    y = R * cos(th) .* (ph - ph0);
end

if nargin > 4
    if strcmp(units,'km')
        x = x/1000;
        y = y/1000;
    end
end
