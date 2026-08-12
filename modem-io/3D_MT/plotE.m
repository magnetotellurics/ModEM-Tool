%   script that plots real and imaginary parts of array C
%  plot real part
axes(h1)
pcolor(X,Y,real(C));
shading interp;
if(P.CLmode == 'manu')
  set(h1,'Clim',P.CLreal);
else
  set(h1,'CLimMode','auto');
end
%caxis(Clims);
set(gca,'FontWeight','demi','FontSize',12)
xlabel(xLab)
ylabel(yLab)
axis(AxMode)
hCB = colorbar;
set(hCB,'FontWeight','Demi','FontSize',12)

%  plot imag part
axes(h2)
pcolor(X,Y,imag(C));
shading interp
if(P.CLmode == 'manu')
  set(h2,'Clim',P.CLimag);
else
  set(h2,'CLimMode','auto');
end
%caxis(Clims);
set(gca,'FontWeight','demi','FontSize',12)
xlabel(xLab)
ylabel(yLab)
axis(AxMode)
hCB = colorbar;
set(hCB,'FontWeight','Demi','FontSize',12)
