function h1 = CondPlotMap(h0,map,rho)

% h1 = CondPlotMap(h0,map,rho)
%
% Plot a map of, say, conductance using CondPlot figure
% (can be easily modified to plot something else too)

P  = get(h0,'UserData');

if findstr(P.cblabel,'log_{10} \rho')
    P.Cond = - P.Cond; % convert to conductivity locally
    P.Clims = fliplr(- P.Clims);
end

if strcmp(map,'S')
    dz = diff(P.Z)*1000;
    S = zeros(length(P.X),length(P.Y));
    for k = P.k1:P.k2-1
        S = S + squeeze((10.^P.Cond(:,:,k))*dz(k));
    end
    S = S(P.i1:P.i2,P.j1:P.j2); S = log10(S); 
    Slims = P.Clims + log10(sum(dz(P.k1:P.k2-1)));
    Sticks = 0.5:0.5:5;
    Stitle = ['Total Conductance from ' num2str(P.Z(P.k1)) ' to ' num2str(P.Z(P.k2)) ' km depth'];
    cblabel = 'log_{10} (S)';
    cmap = P.cmap;
elseif strcmp(map,'D')
    if nargin<3
        rho = 40;
    end
    D = NaN*zeros(length(P.X),length(P.Y));
    for i = P.i1:P.i2
        for j = P.j1:P.j2
            k = find(P.Cond(i,j,P.k1:P.k2)>log10(1/rho),1,'first');
            if isempty(k)
                D(i,j) = NaN;
            else
                D(i,j) = P.Z(P.k1+k);
            end
        end
    end
    S = D(P.i1:P.i2,P.j1:P.j2);
    Slims = [nanmin(nanmin(D)) nanmax(nanmax(D))];
    Sticks = Slims(1):5:Slims(2);
    Stitle = ['Depth to the ' num2str(rho) ' \Omega m interface (' P.units ')'];
    cblabel = ''; %P.units;
    cmap = flipud(P.cmap);
end

if P.PlotLatLon
  S = S.';
  xLab = 'Longitude (degrees)';
  yLab = 'Latitude (degrees)';
  avgX = mean(P.X(P.i1:P.i2));
  avgY = mean(P.Y(P.j1:P.j2));
  [dummy,X] = xy2latlon(avgX,P.Y(P.j1:P.j2),P.lat0,P.lon0,P.units);
  [Y,dummy] = xy2latlon(P.X(P.i1:P.i2),avgY,P.lat0,P.lon0,P.units);
  [Y,X] = meshgrid(Y,X);
  [SiteY,SiteX] = xy2latlon(P.SiteX,P.SiteY,P.lat0,P.lon0,P.units);
  AxMode = 'xy';
else
  xLab = 'Y';
  yLab = 'X';
  X =  P.Y(P.j1:P.j2);
  Y = P.X(P.i1:P.i2);
  SiteY = P.SiteX;
  SiteX = P.SiteY;
  AxMode = 'xy';
end

h1 = figure('Position',[200,400,800,600],...
    'Units','normalized','Color','w');
if P.PlotLatLon
    m_proj('miller',...
        'long',[min(X(:,1)) max(X(:,1))],...
        'lat', [min(Y(1,:)) max(Y(1,:))]);
    m_pcolor(X,Y,S);caxis(Slims); 
    if ~P.nointerp; shading flat; end
    cb=contourcmap(Sticks,'jet','colorbar','on','location','horizontal','FontSize',14);
    %cb=colorbar('FontSize',14); colormap(cmap);
    hold on;
    plotStates
    if P.PlotSites && (P.Slice == 'Z')
        m_plot(SiteX,SiteY,'k^','MarkerSize',10,'LineWidth',2)
    end
    m_grid('box','fancy','tickdir','in'); %m_grid
else
    pcolor(X,Y,S);caxis(Slims);colormap(cmap);
    cb=colorbar('FontSize',14);
    if ~P.nointerp; shading flat; end
    hold on
    if P.PlotSites && (P.Slice == 'Z')
        plot(SiteX,SiteY,'k^','MarkerSize',10,'LineWidth',2);
    end
end
title(Stitle,'FontWeight','demi','FontSize',14);
hold off
set(gca,'FontWeight','demi','FontSize',16)
xlabel(xLab)
ylabel(yLab)
axis(AxMode)
title(cb,cblabel,'FontWeight','demi','FontSize',14);