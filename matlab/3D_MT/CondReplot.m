if(P.Slice == 'Z') && P.PlotLatLon
  C = squeeze(P.Cond(P.i1:P.i2,P.j1:P.j2,P.Np)); C = C.';
  xLab = 'Longitude (degrees)';
  yLab = 'Latitude (degrees)';
  avgX = mean(P.X(P.i1:P.i2));
  avgY = mean(P.Y(P.j1:P.j2));
  [dummy,X] = xy2latlon(avgX,P.Y(P.j1:P.j2),P.lat0,P.lon0,P.units);
  [Y,dummy] = xy2latlon(P.X(P.i1:P.i2),avgY,P.lat0,P.lon0,P.units);
  [Y,X] = meshgrid(Y,X);
  [SiteY,SiteX] = xy2latlon(P.SiteX,P.SiteY,P.lat0,P.lon0,P.units);
  AxMode = 'xy';
elseif (P.Slice == 'Z')
  C = squeeze(P.Cond(P.i1:P.i2,P.j1:P.j2,P.Np));
  xLab = 'Y';
  yLab = 'X';
  X =  P.Y(P.j1:P.j2);
  Y = P.X(P.i1:P.i2);
  SiteY = P.SiteX;
  SiteX = P.SiteY;
  AxMode = 'xy';
elseif (P.Slice == 'X')
  C = squeeze(P.Cond(P.Np,P.j1:P.j2,P.k1:P.k2)).';
  if P.PlotLatLon
    xLab = 'Longitude (degrees)';
    [n,X] = xy2latlon(P.X(P.Np),P.Y(P.j1:P.j2).',P.lat0,P.lon0,P.units);
  else
    xLab = 'Y';
    X = P.Y(P.j1:P.j2);
  end
  yLab = 'Z';
  Y = P.Z(P.k1:P.k2);
  AxMode = 'ij';
elseif(P.Slice == 'Y')
  C = squeeze(P.Cond(P.i1:P.i2,P.Np,P.k1:P.k2)).';
  if P.PlotLatLon
    xLab = 'Latitude (degrees)';
    [X,n] = xy2latlon(P.X(P.i1:P.i2).',P.Y(P.Np),P.lat0,P.lon0,P.units);
  else
    xLab = 'X';
    X = P.X(P.i1:P.i2);
  end
  yLab = 'Z';
  Y = P.Z(P.k1:P.k2);
  AxMode = 'ij';
end
axes(P.ax);
p=get(P.ax,'position');
set(P.ax,'FontWeight','demi','FontSize',20)
if(P.Slice == 'Z') && P.PlotLatLon
    m_proj(P.proj,...
        'long',[min(X(:,1)) max(X(:,1))],...
        'lat', [min(Y(1,:)) max(Y(1,:))]);
    m_pcolor(X,Y,C);caxis(P.Clims); 
    if ~P.nointerp; shading flat; end
    if P.PlotColorBar
        cb=colorbar('FontSize',14); colormap(P.cmap);
    end
    set(P.ax,'position',p)
    hold on;
    plotStates
    if P.PlotSites && (P.Slice == 'Z')
        m_plot(SiteX,SiteY,'k^','MarkerSize',10,'LineWidth',2)
    end
    m_grid('box','fancy','tickdir','in','FontWeight','demi','FontSize',20);
else
    pcolor(X,Y,C);caxis(P.Clims);colormap(P.cmap);
    if P.PlotColorBar
        cb=colorbar('FontSize',14);
    end
    if ~P.nointerp; shading flat; end
    hold on
    if P.PlotSites && (P.Slice == 'Z')
        plot(SiteX,SiteY,'k^','MarkerSize',10,'LineWidth',2);
    end
end
hold off
set(P.ax,'FontWeight','demi','FontSize',24)
if P.Label
    xlabel(xLab)
    ylabel(yLab)
end
axis(AxMode)
if P.PlotColorBar
    title(cb,P.cblabel,'FontWeight','demi','FontSize',14);
end

