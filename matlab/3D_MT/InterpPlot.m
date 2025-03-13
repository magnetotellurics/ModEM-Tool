function [h] = InterpPlot(lon,lat,data,cblabel,htitle,cbaxis,msize)

% [h] = InterpPlot(lon,lat,data,cblabel,htitle,cbaxis,msize)
%
% Interpolates and plots values such as apparent resistivity,
% Niblett-Bostick resistivity, or phase, at a set of sites.
% Uses m_map [ http://www.eos.ubc.ca/~rich/map.html ].
%
% Political boundaries may be obtained from
% http://www.maproom.psu.edu/dcw/
% using "download points" option.
% For these to work, -180 < lon <= 180 should be used.
%
% (c) Anna Kelbert, 2009-10-01
%
% Usage example:
%   InterpPlot(lon,lat,log10(rhoNB),'log_{10} (\rho_{NB}) [\Omega m]','30 km');
%   colormap(flipud(jet));
%   print('-r300','-djpeg',['EMscope_rhoNB_30km']);

i = find(lon>180); lon(i) = lon(i) - 360;
gddel = min((max(lon)-min(lon))/100,(max(lat)-min(lat))/100);
gdlon = min(lon):gddel:max(lon);
gdlat = min(lat):gddel:max(lat);
[gdlat,gdlon] = meshgrid(gdlat,gdlon);
gddata = griddata(lat,lon,data,gdlat,gdlon);

if nargout > 0
    h = clf;
end
% 'lambert' or 'miller'
gdpad = (max(lat)-min(lat))/10;
m_proj('miller',...
    'long',[min(lon)-gdpad max(lon)+gdpad],...
    'lat', [min(lat)-gdpad max(lat)+gdpad]);
m_pcolor(gdlon,gdlat,gddata); shading interp;
cb=colorbar('FontSize',14); colormap jet;
title(cb,cblabel,'FontSize',14,'FontWeight','demi');
set(gca,'FontSize',18,'FontWeight','demi');
if nargin > 5
    caxis([cbaxis(1) cbaxis(2)]);
end
if nargin <= 6
    msize = 5;
end
hold on; m_plot(lon,lat,'k*','MarkerSize',msize,'LineWidth',1.5)
%m_gshhs_i('color','k');              % Coastline...
%m_gshhs_i('speckle','color','k');    % with speckle added
% Plot US state boundaries on an existing map with m_map
% To use, install m_map and add a setenv line to startup, e.g.
% setenv('MAPPATH','/home/mt/anya/ModEM/matlab/Toolbox/m_map/data/');
cbndry  = [0.5 0.5 0.5];
%m_gshhs_i('color',cbndry);              % Coastline... higher resolution
%m_gshhs_i('speckle','color',cbndry);    % with speckle added
m_coast('color',cbndry);              % Coastline...
m_coast('speckle','color',cbndry);    % with speckle added
states = m_shaperead('state_bounds');
for i = 1:length(states.ncst)
    stateslon = states.ncst{i}(:,1);
    stateslat = states.ncst{i}(:,2);
    m_line(stateslon,stateslat,'linewi',2,'color',cbndry);
end
m_grid('box','fancy','tickdir','in');
%ylabel('Latitude (degrees)','FontSize',16,'FontWeight','demi');
%xlabel('Longitude (degrees)','FontSize',16,'FontWeight','demi');
if nargin > 4
    title(htitle,'FontSize',20,'FontWeight','bold');
end

