function h1 = CondPlotSummary2(h0,d1,d2,d3,d4,d5,d6,pts1,m1,pts2,m2)

% h1 = CondPlotSummary2(h0,d1,d2,d3,d4,d5,d6,pts1,m1,pts2,m2)
%
% Make a full page summary plot using CondPlot figure
% Plots six depths d1-d6
% and one cross-sections input by the user manually;
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
dist1 = sqrt((xint-x1(1))^2 + (yint-y1(1))^2); mark1 = ' x ';
dist2 = sqrt((xint-x2(1))^2 + (yint-y2(1))^2); mark2 = ' x ';

% also compute distances to the optional points, if specified
if nargin > 7
    if P.PlotLatLon % convert to lat/lon and back to km (or meters)
        [xpts,ypts] = latlon2xy(pts1(:,1),pts1(:,2),P.lat0,P.lon0,P.units);
    else
        xpts = pts1(:,1); ypts = pts1(:,2);
    end
    distpts1 = sqrt((xpts-x1(1)).^2 + (ypts-y1(1)).^2);
    dist1 = [dist1; distpts1];
    mark1 = [mark1; [' ' m1 ' ']];
end
if nargin > 9
    if P.PlotLatLon % convert to lat/lon and back to km (or meters)
        [xpts,ypts] = latlon2xy(pts2(:,1),pts2(:,2),P.lat0,P.lon0,P.units);
    else
        xpts = pts2(:,1); ypts = pts2(:,2);
    end
    distpts2 = sqrt((xpts-x2(1)).^2 + (ypts-y2(1)).^2);
    dist2 = [dist2; distpts2];
    mark2 = [mark2; [' ' m2 ' ']];
end

h1 = figure('Position',[100,50,1400,1200],...
    'Units','normalized','Color','w');

s1 = subplot(6,2,[1 3])
i1 = find(P.Z <= d1, 1, 'last');
ptitle = [num2str(P.Z(i1)) ' - ' num2str(P.Z(i1+1)) ' km'];
P.PlotColorBar = 0;
P.Label = 0;
P.Slice = 'Z';
P.Np = i1+1;
P.ax = s1;
CondReplot
line(yline1,xline1,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline1(1),xline1(1),'A','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline1(2),xline1(2),'A''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
line(yline2,xline2,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline2(1),xline2(1),'B','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline2(2),xline2(2),'B''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
title(ptitle,'FontWeight','demi','FontSize',14);
hold off

s2 = subplot(6,2,[5 7])
i2 = find(P.Z <= d2, 1, 'last');
ptitle = [num2str(P.Z(i2)) ' - ' num2str(P.Z(i2+1)) ' km'];
P.PlotColorBar = 0;
P.Label = 0;
P.Slice = 'Z';
P.Np = i2+1;
P.ax = s2;
CondReplot
line(yline1,xline1,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline1(1),xline1(1),'A','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline1(2),xline1(2),'A''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
line(yline2,xline2,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline2(1),xline2(1),'B','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline2(2),xline2(2),'B''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
title(ptitle,'FontWeight','demi','FontSize',14);
hold off

s3 = subplot(6,2,[9 11])
i3 = find(P.Z <= d3, 1, 'last');
ptitle = [num2str(P.Z(i3)) ' - ' num2str(P.Z(i3+1)) ' km'];
P.PlotColorBar = 0;
P.Label = 1;
P.Slice = 'Z';
P.Np = i3+1;
P.ax = s3;
CondReplot
line(yline1,xline1,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline1(1),xline1(1),'A','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline1(2),xline1(2),'A''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
line(yline2,xline2,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline2(1),xline2(1),'B','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline2(2),xline2(2),'B''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
title(ptitle,'FontWeight','demi','FontSize',14);
hold off

s4 = subplot(6,2,[2 4])
posn = get(s4,'OuterPosition');
set(s4,'OuterPosition',[0.9*posn(1) posn(2) posn(3) posn(4)]);
i4 = find(P.Z <= d4, 1, 'last');
ptitle = [num2str(P.Z(i4)) ' - ' num2str(P.Z(i4+1)) ' km'];
P.PlotColorBar = 1;
P.Label = 0;
P.Slice = 'Z';
P.Np = i4+1;
P.ax = s4;
CondReplot
line(yline1,xline1,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline1(1),xline1(1),'A','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline1(2),xline1(2),'A''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
line(yline2,xline2,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline2(1),xline2(1),'B','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline2(2),xline2(2),'B''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
title(ptitle,'FontWeight','demi','FontSize',14);
hold off

s5 = subplot(6,2,[6 8])
posn = get(s5,'OuterPosition');
set(s5,'OuterPosition',[0.9*posn(1) posn(2) posn(3) posn(4)]);
i5 = find(P.Z <= d5, 1, 'last');
ptitle = [num2str(P.Z(i5)) ' - ' num2str(P.Z(i5+1)) ' km'];
P.PlotColorBar = 0;
P.Label = 0;
P.Slice = 'Z';
P.Np = i5+1;
P.ax = s5;
CondReplot
line(yline1,xline1,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline1(1),xline1(1),'A','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline1(2),xline1(2),'A''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
line(yline2,xline2,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline2(1),xline2(1),'B','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline2(2),xline2(2),'B''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
title(ptitle,'FontWeight','demi','FontSize',14);
hold off

s6 = subplot(6,2,[10 12])
posn = get(s6,'OuterPosition');
set(s6,'OuterPosition',[0.9*posn(1) posn(2) posn(3) posn(4)]);
i6 = find(P.Z <= d6, 1, 'last');
ptitle = [num2str(P.Z(i6)) ' - ' num2str(P.Z(i6+1)) ' km'];
P.PlotColorBar = 0;
P.Label = 0;
P.Slice = 'Z';
P.Np = i6+1;
P.ax = s6;
CondReplot
line(yline1,xline1,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline1(1),xline1(1),'A','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline1(2),xline1(2),'A''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
line(yline2,xline2,'Color',[0.7 0.7 0.7],'LineWidth',2);
text(yline2(1),xline2(1),'B','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
text(yline2(2),xline2(2),'B''','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16);
title(ptitle,'FontWeight','demi','FontSize',14);
hold off

h2 = figure('Position',[500,100,900,600],...
    'Units','normalized','Color','w');

s21 = subplot(2,1,1)
P.PlotColorBar = 0;
P.Label = 0;
[lon,lat] = m_xy2ll(yline1,xline1);
[xlen,ylen] = latlon2xy(lat,lon,P.lat0,P.lon0,P.units);
linelength = sqrt((ylen(2)-ylen(1))^2+(xlen(2)-xlen(1))^2);
CondPlotSection(P,s21,yline1,xline1,[dist1; 0; linelength-30],[mark1; ' A '; ' A''']);
hold off
 
s22 = subplot(2,1,2)
P.PlotColorBar = 1;
P.Label = 1;
[lon,lat] = m_xy2ll(yline2,xline2);
[xlen,ylen] = latlon2xy(lat,lon,P.lat0,P.lon0,P.units);
linelength = sqrt((ylen(2)-ylen(1))^2+(xlen(2)-xlen(1))^2);
CondPlotSection(P,s22,yline2,xline2,[dist2; 0; linelength-30],[mark2; ' B '; ' B''']);
hold off