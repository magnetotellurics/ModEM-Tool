%   script that plots either real and imaginary parts of array C
%   or amplitude and phase ... or for impedances apparent resistivity
%   and phase
if P.AmpPhase
   %  plot log_10(rho)/phi instead of amp/phi unless this is Hz TF
   P.Rho = P.Comp < 3;
   amp = abs(C).^2;
   %   phase  : need to account for assumption about +- wt in modeling code
   %if(P.Comp ~= P.Mode & P.Comp <3)
   %   phi = atan(-imag(C)./real(C))*180/pi;
   %else
   %   phi = atan2(-imag(C),real(C))*180/pi;
   %end
   phi = unwrap(atan2(-imag(C),real(C)))*180/pi;
   if(P.Rho)
     if(P.Slice == 'T')
        rho = P.T(P.Np)*amp/5;
     else
        rho = amp;
        for  k = 1:nPer
           rho(k,:) = P.T(k)*rho(k,:)/5.;
        end
     end
   end
end

%  plot real part (or amp/rho) in top panel
axes(h1);
if Mplot && (P.Slice == 'T'); ax = worldmap(P.Xlim,P.Ylim); end
if P.AmpPhase
   if P.Rho
      if Mplot && (P.Slice == 'T'); pcolorm(P.X,P.Y,log10(rho)); else pcolor(X,Y,log10(rho)); shading interp; end
      ctitle = 'Apparent Resistivity';
   else
      if Mplot && (P.Slice == 'T'); pcolorm(P.X,P.Y,amp); else pcolor(X,Y,amp); shading interp; end
      ctitle = 'Amplitude';
   end
else
   if Mplot && (P.Slice == 'T'); pcolorm(P.X,P.Y,real(C)); else pcolor(X,Y,real(C)); shading interp; end
   ctitle = 'Real Part';
end
%shading interp
if(P.CLmode == 'manu')
  set(h1,'Clim',P.CLreal);
else
  set(h1,'CLimMode','auto');
end
hold on
if P.PlotSites && (P.Slice == 'T')
    if Mplot
        plotm(P.SiteX,P.SiteY,'k*','MarkerSize',10,'LineWidth',1.8)
    else 
        plot(P.SiteY,P.SiteX,'k*','MarkerSize',10,'LineWidth',1.8)
    end
end
hold on
if Mplot && (P.Slice == 'T')
    S = shaperead('landareas','UseGeoCoords',true);
    geoshow([S.Lat], [S.Lon],'Color','black');
    latlim = getm(ax, 'MapLatLimit');
    lonlim = getm(ax, 'MapLonLimit');
    states = shaperead('usastatehi',...
        'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
    stateslat = [states.LabelLat];
    stateslon = [states.LabelLon];
    geoshow([states.Lat], [states.Lon], 'Color', 'black')
    tf = ingeoquad(stateslat, stateslon, latlim, lonlim);
    textm(stateslat(tf), stateslon(tf), {states(tf).Name}, ...
        'Color','white','FontWeight','Demi')
end
hold off
set(gca,'FontWeight','demi','FontSize',12)
xlabel(xLab)
ylabel(yLab)
text(.1,.9,ctitle,'Units','Normalized',...
	'FontWeight','demi','FontSize',14)
%title(ctitle)
axis(AxMode)
hCB = colorbar;
axes(hCB);
set(gca,'FontWeight','Demi','FontSize',12)

%  plot imag part (or phase) in lower panel
axes(h2);
if Mplot && (P.Slice == 'T'); ax = worldmap(P.Xlim,P.Ylim); end;
if P.AmpPhase
   if Mplot && (P.Slice == 'T'); pcolorm(P.X,P.Y,phi); else pcolor(X,Y,phi); shading interp; end
   ctitle = 'Phase';
else
   if Mplot && (P.Slice == 'T'); pcolorm(P.X,P.Y,imag(C)); else pcolor(X,Y,imag(C)); shading interp; end
   ctitle = 'Imag Part';
end
%shading interp
if(P.CLmode == 'manu')
  set(h2,'Clim',P.CLimag);
else
  set(h2,'CLimMode','auto');
end
hold on
if P.PlotSites && (P.Slice == 'T')
    if Mplot
        plotm(P.SiteX,P.SiteY,'k*','MarkerSize',10,'LineWidth',1.8)
    else 
        plot(P.SiteY,P.SiteX,'k*','MarkerSize',10,'LineWidth',1.8)
    end
end
hold on
if Mplot && (P.Slice == 'T')
    S = shaperead('landareas','UseGeoCoords',true);
    geoshow([S.Lat], [S.Lon],'Color','black');
    latlim = getm(ax, 'MapLatLimit');
    lonlim = getm(ax, 'MapLonLimit');
    states = shaperead('usastatehi',...
        'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
    stateslat = [states.LabelLat];
    stateslon = [states.LabelLon];
    geoshow([states.Lat], [states.Lon], 'Color', 'black')
    tf = ingeoquad(stateslat, stateslon, latlim, lonlim);
    textm(stateslat(tf), stateslon(tf), {states(tf).Name}, ...
        'Color','white','FontWeight','Demi')
end
hold off
%caxis(Clims);
set(gca,'FontWeight','demi','FontSize',12)
xlabel(xLab)
ylabel(yLab)
text(.1,.9,ctitle,'Units','Normalized',...
	'FontWeight','demi','FontSize',14)
%title(ctitle)
axis(AxMode)
hCB = colorbar;
axes(hCB);
set(gca,'FontWeight','Demi','FontSize',12)
cmap = jet(64);
cmap = cmap(:,end:-1:1);
colormap(cmap)

