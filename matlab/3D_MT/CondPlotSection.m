function [h1,cb] = CondPlotSection(h0,h1,yline,xline,pts,marks)

% [h1,cb] = CondPlotSection(h0,h1,yline,xline,pts,marks)
%
% Plot a cross section at an angle using CondPlot figure
% can also be used to plot a predefined cross-section on any figure
% by passing the structure P instead of h0
% pts is an optional array of distances from zero which should be marked
% marks are one-character text strings

if ishandle(h0)
    P  = get(h0,'UserData');
    ax = get(h0,'CurrentAxes');
    [yline,xline] = getline(ax);
else
    P = h0;
    if nargin < 4
        error('Line coordinates not provided')
    end
end


if P.PlotLatLon % convert to lat/lon and back to km (or meters)
    [lon,lat] = m_xy2ll(yline,xline);
    [xline,yline] = latlon2xy(lat,lon,P.lat0,P.lon0,P.units);
else % just get lat/lon for text annotation
    [lat,lon] = xy2latlon(xline,yline,P.lat0,P.lon0,P.units);
end

dy = (yline(2)-yline(1))/100;
newy = yline(1):dy:yline(2);

dx = (xline(2)-xline(1))/100;
newx = xline(1):dx:xline(2);

% oldx = oldgrid.origin(1)+cumsum([0; oldgrid.dx]);
% oldy = oldgrid.origin(2)+cumsum([0; oldgrid.dy]);
% oldz = oldgrid.origin(3)+cumsum([0; oldgrid.dz]);
oldx = P.X(P.i1:P.i2);
oldy = P.Y(P.j1:P.j2);
oldz = P.Z(P.k1:P.k2);
oldv = P.Cond(P.i1:P.i2,P.j1:P.j2,P.k1:P.k2);

% [oldY,oldX,oldZ] = meshgrid(oldy(1:end-1)+diff(oldy)/2,...
%     oldx(1:end-1)+diff(oldx)/2,...
%     oldz(1:end-1)+diff(oldz)/2);
[oldY,oldX,oldZ] = meshgrid(oldy,oldx,oldz);

newz = oldz(1:end-1)+diff(oldz)/2;


newv(length(newx),length(newz)) = 0;
for i = 1:length(newx)
    temp = interp3(oldY,oldX,oldZ,oldv,newy(i),newx(i),newz);
    newv(i,:) = squeeze(temp)';
end

[Y,X] = meshgrid(newz,sqrt((newx-xline(1)).^2+(newy-yline(1)).^2));

if nargin < 2
    h1 = figure('Position',[200,400,1400,600],...
        'Units','normalized','Color','w');
else
    axes(h1);
end
pcolor(X,Y,newv);caxis(P.Clims);colormap(P.cmap);
if P.PlotColorBar
    cb=colorbar('FontSize',14);
    ylabel(cb,P.cblabel,'FontWeight','demi','FontSize',14);
end
if P.flipud; colormap(flipud(P.cmap)); end
axis ij; shading flat
if nargin > 4
    hold on;
    text(pts,-10*ones(size(pts)),marks,'FontWeight','demi','FontSize',14);
end
hold off
set(gca,'FontWeight','demi','FontSize',16)
if P.Label
    xlabel(['Distance along the profile (' P.units ')'])
    ylabel(['Depth (' P.units ')'])
end

if nargin < 2
    pt1 = ['( ' num2str(lat(1)) ', ' num2str(lon(1)) ' )'];
    text(-0.05,1.05,pt1,'FontWeight','demi','FontSize',10,...
        'Color','b','units','normalized');
    pt2 = ['( ' num2str(lat(2)) ', ' num2str(lon(2)) ' )'];
    text(0.95,1.05,pt2,'FontWeight','demi','FontSize',10,...
        'Color','b','units','normalized');
end