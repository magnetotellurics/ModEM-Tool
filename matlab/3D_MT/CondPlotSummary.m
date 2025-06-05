function h1 = CondPlotSummary(h0,d1,d2,pts1,m1,pts2,m2)

% h1 = CondPlotSummary(h0,d1,d2,pts1,m1,pts2,m2)
%
% Make a summary plot using CondPlot figure
% Plots two depths d1 and d2
% and two cross-sections input by the user manually;
% optionally annotates the cross-sections with pts
% marked with the strings m
%
% Yellowstone usage example:
% CondPlotSummary(h0,25,50,[44.6 -110.3],'Y')

P  = get(h0,'UserData');
ax = get(h0,'CurrentAxes');

[yline1,xline1] = getline(ax);
[yline2,xline2] = getline(ax);

% compute distances to the point of intersection in km (or meters)
[yint,xint]=polyxpoly(yline1,xline1,yline2,xline2)
if P.PlotLatLon % convert to lat/lon and back to km (or meters)
    [lon,lat] = m_xy2ll(yline1,xline1);
    [x1,y1] = latlon2xy(lat,lon,P.lat0,P.lon0,P.units);
    [lon,lat] = m_xy2ll(yline2,xline2);
    [x2,y2] = latlon2xy(lat,lon,P.lat0,P.lon0,P.units);
    [lon,lat] = m_xy2ll(yint,xint);
    [xint,yint] = latlon2xy(lat,lon,P.lat0,P.lon0,P.units);
end
dist1 = sqrt((xint-x1(1))^2 + (yint-y1(1))^2); mark1 = 'x';
dist2 = sqrt((xint-x2(1))^2 + (yint-y2(1))^2); mark2 = 'x';

% also compute distances to the optional points, if specified
if nargin > 3
    if P.PlotLatLon % convert to lat/lon and back to km (or meters)
        [xpts,ypts] = latlon2xy(pts1(:,1),pts1(:,2),P.lat0,P.lon0,P.units);
    else
        xpts = pts1(:,1); ypts = pts1(:,2);
    end
    distpts1 = sqrt((xpts-x1(1)).^2 + (ypts-y1(1)).^2);
    dist1 = [dist1; distpts1];
    mark1 = [mark1; m1];
end
if nargin > 5
    if P.PlotLatLon % convert to lat/lon and back to km (or meters)
        [xpts,ypts] = latlon2xy(pts2(:,1),pts2(:,2),P.lat0,P.lon0,P.units);
    else
        xpts = pts2(:,1); ypts = pts2(:,2);
    end
    distpts2 = sqrt((xpts-x2(1)).^2 + (ypts-y2(1)).^2);
    dist2 = [dist2; distpts2];
    mark2 = [mark2; m2];
end

h1 = figure('Position',[100,100,1400,800],...
    'Units','normalized','Color','w');

s1 = subplot(3,2,[1 3])
i1 = find(P.Z <= d1, 1, 'last');
ptitle = [num2str(P.Z(i1)) ' - ' num2str(P.Z(i1+1)) ' km'];
P.PlotColorBar = 0;
P.Slice = 'Z';
P.Np = i1+1;
P.ax = s1;
CondReplot
line(yline1,xline1,'Color',[0.7 0.7 0.7],'LineWidth',2);
line(yline2,xline2,'Color',[0.7 0.7 0.7],'LineWidth',2);
title(ptitle,'FontWeight','demi','FontSize',14);
hold off

s2 = subplot(3,2,[2 4])
i2 = find(P.Z <= d2, 1, 'last');
ptitle = [num2str(P.Z(i2)) ' - ' num2str(P.Z(i2+1)) ' km'];
P.PlotColorBar = 1;
P.Slice = 'Z';
P.Np = i2+1;
P.ax = s2;
CondReplot
line(yline1,xline1,'Color',[0.7 0.7 0.7],'LineWidth',2);
line(yline2,xline2,'Color',[0.7 0.7 0.7],'LineWidth',2);
title(ptitle,'FontWeight','demi','FontSize',14);
hold off

s3 = subplot(3,2,5)
P.PlotColorBar = 0;
CondPlotSection(P,s3,yline1,xline1,dist1,mark1)
hold off
 
s4 = subplot(3,2,6)
P.PlotColorBar = 1;
CondPlotSection(P,s4,yline2,xline2,dist2,mark2)
hold off