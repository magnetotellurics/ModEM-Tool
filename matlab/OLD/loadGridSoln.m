gridFile = 'CondTest';
solnFile = 'EsolTest.out';

%  load in grid (including conductivity)
[Dy,Dz,Cond,periods,Nza] = rdGridCond(gridFile);
y = cumsum([0;Dy]);
z = cumsum([0;Dz]);
zCenter = z-z(Nza+1);
Nmid = fix((length(y)+1)/2);
yCenter = y-y(Nmid);
Cond = [Cond Cond(:,end)];
Cond = [Cond ;Cond(end,:)];
figure
pcolor(yCenter,zCenter,log10(Cond')); shading flat; axis('ij');

caxis([-4,-.5])
colorbar

%  load in solution (so far just for one period)
[E] = readVec2D(solnFile);
surfaceAmp = mean(abs(E.v(:,Nza+1)));
Z = [0,0]; Y = [min(yCenter),max(yCenter)];
figure
pcolor(yCenter,zCenter,real(E.v)'/surfaceAmp); shading flat; axis('ij');
hold on;  plot(Y,Z,'lineWidth',2,'Color','k');
caxis([-.5,2])
colorbar
figure
pcolor(yCenter,zCenter,imag(E.v/surfaceAmp)'); shading flat; axis('ij');
hold on;  plot(Y,Z,'lineWidth',2,'Color','k');
caxis([0,1])
colorbar

%  find magnetic field
omega = 2*pi/periods(1);
B = TEmag(E.v,Dz,omega);
BsurfaceAmp = mean(mean(abs(B(:,Nza:Nza+1))));
zEdge = (zCenter(1:end-1)+zCenter(2:end))/2;
figure
pcolor(yCenter,zEdge,real(B)'/BsurfaceAmp); shading flat; axis('ij');
hold on;  plot(Y,Z,'lineWidth',2,'Color','k');
colorbar

figure
pcolor(yCenter,zEdge,imag(B)'/BsurfaceAmp); shading flat; axis('ij');
hold on;  plot(Y,Z,'lineWidth',2,'Color','k');
colorbar

Esurface = E.v(:,Nza+1); 
Bsurface = (B(:,Nza)+B(:,Nza+1))/2;
Zsurface = Esurface./Bsurface;
%  change to standard impedance units
Zsurface = Zsurface/1000;

Esites = [ Esurface(50) (Esurface(50)+Esurface(51))/2 ...
	Esurface(51) (2*Esurface(51)+Esurface(52))/3 Esurface(52)];
Bsites = [ Bsurface(50) (Bsurface(50)+Bsurface(51))/2 ...
	Bsurface(51) (2*Bsurface(51)+Bsurface(52))/3 Bsurface(52)];
Zsites = Esites./Bsites;
