function [rhoOut,maskOut] = regridModel(gridIn,gridOut,rhoIn,maskIn);

seaRes = 0.3;
%   refine output grid to improve accuracy of interpolation
nRef = 3;
minFracOcean = .4;

%   coordinates of input grid cell edges
x = gridIn.origin(1)+[0; cumsum(gridIn.dx)];
y = gridIn.origin(2)+[0; cumsum(gridIn.dy)];
z = gridIn.origin(3)+[0; cumsum(gridIn.dz)];
%   coordinates of output grid cell edges
x1 = gridOut.origin(1)+[0; cumsum(gridOut.dx)];
y1 = gridOut.origin(2)+[0; cumsum(gridOut.dy)];
z1 = gridOut.origin(3)+[0; cumsum(gridOut.dz)];
%   coordinates of input grid cell centers
xc = x(1:end-1)+ gridIn.dx/2;
yc = y(1:end-1)+ gridIn.dy/2;
zc = z(1:end-1)+ gridIn.dz/2;

temp = [gridOut.dx/nRef gridOut.dx/nRef gridOut.dx/nRef];
[n1,n2] = size(temp);
dx2 = reshape(temp',n1*n2,1);
temp = [gridOut.dy/nRef gridOut.dy/nRef gridOut.dy/nRef];
[n1,n2] = size(temp);
dy2 = reshape(temp',n1*n2,1);
temp = [gridOut.dz/nRef gridOut.dz/nRef gridOut.dz/nRef];
[n1,n2] = size(temp);
dz2 = reshape(temp',n1*n2,1);

%   coordinates of refined output grid cell edges
x2 = gridOut.origin(1)+[0; cumsum(dx2)];
y2 = gridOut.origin(2)+[0; cumsum(dy2)];
z2 = gridOut.origin(3)+[0; cumsum(dz2)];
%   coordinates of refined output grid cell centers
xc2 = x2(1:end-1)+ dx2/2;
yc2 = y2(1:end-1)+ dy2/2;
zc2 = z2(1:end-1)+ dz2/2;

% interpolate mask
[XC,YC,ZC] = ndgrid(xc,yc,zc);
[XC2,YC2,ZC2] = ndgrid(xc2,yc2,zc2);
mask2 = interpn(XC,YC,ZC,maskIn,XC2,YC2,ZC2,'nearest');
rho2 = interpn(XC,YC,ZC,log10(rhoIn),XC2,YC2,ZC2,'linear');

ix = [1:nRef:length(xc2)-nRef+1]-1;
iy = [1:nRef:length(yc2)-nRef+1]-1;
iz = [1:nRef:length(zc2)-nRef+1]-1;
nx = length(ix);
ny = length(iy);
nz = length(iz);

rhoOut = zeros(nx,ny,nz);
maskOut = zeros(nx,ny,nz);
n = zeros(nx,ny,nz);
n2 = ones(size(rho2));
n2(isnan(rho2)) = 0;
rho2(isnan(rho2)) = 0;

for j = 1:nRef
   for k = 1:nRef
      for l = 1:nRef
         rhoOut = rhoOut+rho2(ix+j,iy+k,iz+l)...
		.*(1-mask2(ix+j,iy+k,iz+l));
         maskOut = maskOut+mask2(ix+j,iy+k,iz+l);
         n = n + n2(ix+j,iy+k,iz+l) ...
		.*(1-mask2(ix+j,iy+k,iz+l));
      end
   end
end
rhoOut = 10.^(rhoOut./n);
maskOut(maskOut<nRef.^3*minFracOcean) = 0;
maskOut(maskOut>0) = 1;
rhoOut(maskOut==1) = seaRes;
