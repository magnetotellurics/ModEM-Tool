makeSlab
load CascadiaInv
latLims = [42,49];
lonLims = [236,244];
[ny,nx,nz] = size(Res.X);
iY = fix(ny/2);
temp = abs(Res.X(iY,:,1)-lonLims(1));
iX1 = find(min(temp)==temp);
temp = abs(Res.X(iY,:,1)-lonLims(2));
iX2 = find(min(temp)==temp);
temp = abs(Res.Y(:,1,1)-latLims(1));
iY1 = find(min(temp)==temp);
temp = abs(Res.Y(:,1,1)-latLims(2));
iY2 = find(min(temp)==temp);

latsUse = Res.lat(iY1:iY2,iX1:iX2);
lonsUse = Res.lon(iY1:iY2,iX1:iX2);
[nyUse,nxUse] = size(latsUse);

zSlab = zeros(nyUse,nxUse);

for j = 1:nyUse
   for k = 1:nxUse
      zSlab(j,k) = interpSlab(slabs,depths,latsUse(j,k),lonsUse(j,k));
   end
end
zTrans = .25*ones(size(zSlab));
surface(lonsUse,latsUse,zSlab,...
   'facecolor',[.5,.5,.5],'edgecolor','none', ...
      'alphadata',zTrans)

