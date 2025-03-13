gridFile = 'Test1.grd';
condFile = 'Test1.cpr';
TEZFile = 'TestOut.imp';
solnFile = 'Test1a.sol';
PLOTfields = 0;

%  load in grid 
[gridDef] = readGrid2D(gridFile);
Ny = gridDef.Ny;
Nz = gridDef.Nz;
Nza = gridDef.Nza;
%  conductivity
[cond] = readCond2D(condFile);

y = cumsum([0 gridDef.Dy])';
z = cumsum([0 gridDef.Dz])';
zCenter = z-z(Nza+1);
Nmid = fix((length(y)+1)/2);
yCenter = y-y(Nmid);
Cond = ones(Ny,Nz)*cond.AirCond;
Cond(:,Nza+1:end) = cond.v;
Cond = [Cond Cond(:,end)];
Cond = [Cond ;Cond(end,:)];
figure
pcolor(yCenter,zCenter,log10(Cond')); shading flat; axis('ij');

caxis([-4,-.5])
colorbar

%  load in solution for all periods
[VEC,header] = readAllVec2D(solnFile);
nPer = length(VEC)

% load in measured impedances for all periods/sites
[allData] = readTEZ(TEZFile);

for iPer = 1:nPer
   E = VEC{iPer};
   surfaceAmp = mean(abs(E.v(:,Nza+1)));
   Z = [0,0]; Y = [min(yCenter),max(yCenter)];
   if PLOTfields
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
   end

   %  find magnetic field
   omega = 2*pi/allData{iPer}.T;
   B = TEmag(E.v,gridDef.Dz,omega);
   BsurfaceAmp = mean(mean(abs(B(:,Nza:Nza+1))));
   zEdge = (zCenter(1:end-1)+zCenter(2:end))/2;
   if PLOTfields
      figure
      pcolor(yCenter,zEdge,real(B)'/BsurfaceAmp); shading flat; axis('ij');
      hold on;  plot(Y,Z,'lineWidth',2,'Color','k');
      colorbar
   
      figure
      pcolor(yCenter,zEdge,imag(B)'/BsurfaceAmp); shading flat; axis('ij');
      hold on;  plot(Y,Z,'lineWidth',2,'Color','k');
      colorbar
   end

   Esurface = E.v(:,Nza+1); 
   Bsurface = (B(:,Nza)+B(:,Nza+1))/2;
   Zsurface = Esurface./Bsurface;
   %  change to standard impedance units
   Zsurface = Zsurface/1000;
   figure
   plot(y,real(Zsurface),'b-');
   hold on; plot(y,-imag(Zsurface),'r--')
   ySites = allData{iPer}.siteLoc(1,:);
   Z = allData{iPer}.Z/1000; 
   plot(ySites,real(Z),'b*')
   plot(ySites,-imag(Z),'ro')

   %  plot apparent resistivity and phase
   figure
   rhoSurface = allData{iPer}.T*(abs(Zsurface).^2)/5;
   rho = allData{iPer}.T*(abs(Z).^2)/5;
   plot(y,rhoSurface,'b-'); hold on; plot(ySites,rho,'b*')
   phi = atan(-imag(Z)./real(Z))*180/pi;
   phiSurface = atan(-imag(Zsurface)./real(Zsurface))*180/pi;
   figure
   plot(y,phiSurface,'b-'); hold on; plot(ySites,phi,'b*')
end 
