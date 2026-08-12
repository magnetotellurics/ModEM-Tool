fdir = '.';
filt = [fdir '/*.soln'];

%  then get solution file
[filename, pathname] = uigetfile(filt, 'Solution File');
cfile = [pathname filename];
% Ex mode
fNum = 1; mNum=2;
[E1,T,Grid,Modes] = rdExyz(cfile, fNum, mNum);
% Ey mode
fNum = 1; mNum=1;
[E2] = rdExyz(cfile, fNum, mNum);

% compute magnetic fields H in A/m
mu0 = 4*pi*10^(-7);
H1 = CurlE(E1,Grid); 
H1.x = H1.x/(-1i*E1.omega*mu0);
H1.y = H1.y/(-1i*E1.omega*mu0);
H1.z = H1.z/(-1i*E1.omega*mu0);
H2 = CurlE(E2,Grid); 
H2.x = H2.x/(-1i*E1.omega*mu0);
H2.y = H2.y/(-1i*E1.omega*mu0);
H2.z = H2.z/(-1i*E1.omega*mu0);

% compute cell centers
NzAir = Grid.NzAir;
xEdge = Grid.origin(1)+[0 ; cumsum(Grid.dx)];
yEdge = Grid.origin(2)+[0 ; cumsum(Grid.dy)];
zEdge = Grid.origin(3)+[0 ; cumsum(Grid.dz)];
zSurface = zEdge(NzAir);
zEdge = zEdge-zSurface;
del = [ Grid.dx(1); Grid.dx(1:end-1)+Grid.dx(2:end) ]/2;
xCenter = Grid.origin(1)+cumsum(del);
del = [ Grid.dy(1); Grid.dy(1:end-1)+Grid.dy(2:end) ]/2;
yCenter = Grid.origin(2)+cumsum(del);
del = [ Grid.dz(1); Grid.dz(1:end-1)+Grid.dz(2:end) ]/2;
zCenter = Grid.origin(3)+cumsum(del)-zSurface;

% take a Z-slice at the surface
Np = Grid.NzAir + 1;
Surface.E1.x = squeeze(E1.x(:,:,Np));
Surface.E1.y = squeeze(E1.y(:,:,Np));
Surface.E1.z = squeeze(E1.z(:,:,Np));
Surface.H1.x = squeeze(H1.x(:,:,Np));
Surface.H1.y = squeeze(H1.y(:,:,Np));
Surface.H1.z = squeeze(H1.z(:,:,Np));
Surface.E2.x = squeeze(E2.x(:,:,Np));
Surface.E2.y = squeeze(E2.y(:,:,Np));
Surface.E2.z = squeeze(E2.z(:,:,Np));
Surface.H2.x = squeeze(H2.x(:,:,Np));
Surface.H2.y = squeeze(H2.y(:,:,Np));
Surface.H2.z = squeeze(H2.z(:,:,Np));

% now, map to cell centers
[X,Y] = meshgrid(xCenter,yCenter);
[Xold,Yold] = meshgrid(xEdge,yCenter);
Surface.H1.x = interp2(Xold,Yold,Surface.H1.x.',X',Y');
Surface.H2.x = interp2(Xold,Yold,Surface.H2.x.',X',Y');
Surface.E1.y = interp2(Xold,Yold,Surface.E1.y.',X',Y');
Surface.E2.y = interp2(Xold,Yold,Surface.E2.y.',X',Y');
[Xold,Yold] = meshgrid(xCenter,yEdge);
Surface.H1.y = interp2(Xold,Yold,Surface.H1.y.',X',Y');
Surface.H2.y = interp2(Xold,Yold,Surface.H2.y.',X',Y');
Surface.E1.x = interp2(Xold,Yold,Surface.E1.x.',X',Y');
Surface.E2.x = interp2(Xold,Yold,Surface.E2.x.',X',Y');
Xnew = X; Ynew = Y;

% compute impedances  (invert horizontal H matrix using Kramer's rule)
det = Surface.H1.x.*Surface.H2.y-Surface.H2.x.*Surface.H1.y;
Inv.H1.x =  Surface.H2.y./det;
Inv.H2.y =  Surface.H1.x./det;
Inv.H2.x = -Surface.H2.x./det;
Inv.H1.y = -Surface.H1.y./det;

F = ImpUnits('[V/m]/[A/m]','[mV/km]/[nT]');
Zxx = F * (Surface.E1.x.*Inv.H1.x + Surface.E2.x.*Inv.H1.y);
Zxy = F * (Surface.E1.x.*Inv.H2.x + Surface.E2.x.*Inv.H2.y);
Zyx = F * (Surface.E1.y.*Inv.H1.x + Surface.E2.y.*Inv.H1.y);
Zyy = F * (Surface.E1.y.*Inv.H2.x + Surface.E2.y.*Inv.H2.y);
%figure; pcolor(Ynew,Xnew,real(Zxy)'); colorbar

AppResXY = T(fNum)*(abs(Zxy).^2)/5;
PhaseXY = (180/pi)*atan(-imag(Zxy)./real(Zxy));
AppResYX = T(fNum)*(abs(Zyx).^2)/5;
PhaseYX = (180/pi)*atan(-imag(Zyx)./real(Zyx));

% plot XY
figure; 
subplot(2,1,1);
pcolor(Ynew,Xnew,log10(AppResXY')); colorbar;
title('Zxy: Apparent Resistivity (log_{10} \rho)','fontsize',16);
cmap = jet(64);
cmap = cmap(:,end:-1:1);
colormap(cmap)
subplot(2,1,2);
pcolor(Ynew,Xnew,PhaseXY'); colorbar;
title('Zxy: Phase (\phi)','fontsize',16);
% plot YX
figure; 
subplot(2,1,1);
pcolor(Ynew,Xnew,log10(AppResYX')); colorbar;
title('Zyx: Apparent Resistivity (log_{10} \rho)','fontsize',16);
cmap = jet(64);
cmap = cmap(:,end:-1:1);
colormap(cmap)
subplot(2,1,2);
pcolor(Ynew,Xnew,PhaseYX'); colorbar;
title('Zyx: Phase (\phi)','fontsize',16);
