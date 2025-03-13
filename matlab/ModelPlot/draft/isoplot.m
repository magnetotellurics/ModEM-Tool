%model = readCond_3D('Yellowstone_14km_ErrFl5_14freq_restart3_NLCG_064.rho',2);
%model = readCond_3D('Yellowstone_10km_errfl5T3_200ohmm_smooth_NLCG_175.rho',2);
dataFile = 'Small_USA_5%_3%_run3_NLCG_027_rewrite.dat';
%dataFile = 'Yellowstone_14freq_paper_errfl5T3.dat';
%modelFile = 'Yellowstone_10km_errfl5T3_smooth_NLCG_049.rho';
modelFile = 'Nested_NWUSA_Z_5%_Hz_3%_run3_NLCG_027.rho';
%modelFile = 'Yellowstone_10km_errfl5T3_smooth_NLCG_049.rho';
model = readCond_3D(modelFile,2);
if strcmp(model.paramType,'LOGE')
    model.paramType = 'LOG10';
    model.v = model.v / log(10);
    model.AirCond = model.AirCond / log(10);
elseif strcmp(model.paramType,'LINEAR')
    model.paramType = 'LOG10';
    model.v = log10(model.v);
    model.AirCond = log10(model.AirCond);
end
depth = cumsum(model.grid.dz); zlen = sum(model.grid.dz);
x = cumsum(model.grid.dx); xlen = sum(model.grid.dx);
y = cumsum(model.grid.dy); ylen = sum(model.grid.dy);
[y,x,z] = meshgrid(model.grid.origin(2)+y,model.grid.origin(1)+x,-depth);
origin = [42.016 -116.477 0];
[lat,lon] = xy2latlon(x,y,origin(1),origin(2),'km');
% approximate longitude to mid grid to obtain a regular gridding
nlat = size(lon,1);
nlon = size(lon,2);
midlon = lon(nlat/2,:,:);
for i=1:nlat
    lon(i,:,:) = midlon;
%    lon(i,:,:) = lon(1,:,:);
end
% actually, Yellowstone is 44.4,-110.7 but the projection is screwed up!
[xY,yY] = latlon2xy(44.6,-110.1,origin(1),origin(2),'km');
%% geographic coords large domain
figure
zm = 7; izm = find(depth >= zm, 1, 'first');
zp = 200; izp = find(depth <= zp, 1, 'last');
iz = izm:izp;
ym = -119; iym = 1; %find(lon(nlat/2,:,1) >= ym, 1, 'first');
yp = -106; iyp = nlon; %find(lon(nlat/2,:,1) <= yp, 1, 'last');
iy = iym:iyp;
xm = 40; ixm = 1; %find(lat(:,1,1) >= xm, 1, 'first');
xp = 47; ixp = nlat; %find(lat(:,1,1) <= xp, 1, 'last');
ix = ixm:ixp;
hpatch = patch(isosurface(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),model.v(ix,iy,iz),-2));
isonormals(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),model.v(ix,iy,iz),hpatch)
set(hpatch,'FaceColor','red','EdgeColor','none')
hold on
hpatch2 = patch(isosurface(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),-model.v(ix,iy,iz),3));
isonormals(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),model.v(ix,iy,iz),hpatch2)
set(hpatch2,'FaceColor','blue','EdgeColor','none')
daspect([1,1,40])
view([10,5])
axis tight
box on
camlight left; 
set(gcf,'Renderer','zbuffer'); lighting phong
hold on;
data = readZ_3D('Yellowstone_14freq_paper_errfl5T3.dat','[mV/km]/[nT]');
xdat = data{2}.siteLoc(:,1)./1000;
ydat = data{2}.siteLoc(:,2)./1000;
zdat = data{2}.siteLoc(:,3)./1000;
[latdat,londat] = xy2latlon(xdat,ydat,origin(1),origin(2),'km');
scatter3(londat,latdat,zdat,'filled');
scatter3(-110.3,44.6,0,400,'k','LineWidth',5);
%% geographic coords close up
figure
zm = 7; izm = find(depth >= zm, 1, 'first');
zp = 200; izp = find(depth <= zp, 1, 'last');
iz = izm:izp;
ym = -118; iym = find(lon(nlat/2,:,1) >= ym, 1, 'first');
yp = -106; iyp = nlon; %find(lon(nlat/2,:,1) <= yp, 1, 'last');
iy = iym:iyp;
xm = 40; ixm = 1; %find(lat(:,1,1) >= xm, 1, 'first');
xp = 47; ixp = nlat; %find(lat(:,1,1) <= xp, 1, 'last');
ix = ixm:ixp;
hpatch = patch(isosurface(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),model.v(ix,iy,iz),-1.48));
isonormals(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),model.v(ix,iy,iz),hpatch)
set(hpatch,'FaceColor','red','EdgeColor','none')
hpatch2 = patch(isosurface(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),-model.v(ix,iy,iz),3));
isonormals(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),model.v(ix,iy,iz),hpatch2)
set(hpatch2,'FaceColor','blue','EdgeColor','none')
daspect([1,1,20])
%view([10,20])
view([-30,40])
axis tight
box off
camlight right; 
set(gcf,'Renderer','zbuffer'); lighting phong
hold on;
data = readZ_3D('Yellowstone_14freq_paper_errfl5T3.dat','[mV/km]/[nT]');
xdat = data{2}.siteLoc(:,1)./1000;
ydat = data{2}.siteLoc(:,2)./1000;
zdat = data{2}.siteLoc(:,3)./1000;
[latdat,londat] = xy2latlon(xdat,ydat,origin(1),origin(2),'km');
scatter3(londat,latdat,zdat,'filled');
scatter3(-110.3,44.6,0,400,'k','LineWidth',5);
%% full domain
figure
zm = 7; izm = find(depth >= zm, 1, 'first');
zp = 800; izp = find(depth <= zp, 1, 'last');
iz = izm:izp;
hpatch = patch(isosurface(y(:,:,iz),x(:,:,iz),z(:,:,iz),model.v(:,:,iz),-2));
isonormals(y(:,:,iz),x(:,:,iz),z(:,:,iz),model.v(:,:,iz),hpatch)
set(hpatch,'FaceColor','red','EdgeColor','none')
daspect([ylen,xlen,zlen])
view([10,20])
axis tight
camlight left; 
set(gcf,'Renderer','zbuffer'); lighting phong
hold on;
data = readZ_3D('Yellowstone_14freq_paper_errfl5T3.dat','[mV/km]/[nT]');
xdat = data{2}.siteLoc(:,1)./1000;
ydat = data{2}.siteLoc(:,2)./1000;
zdat = data{2}.siteLoc(:,3)./1000;
scatter3(ydat,xdat,zdat,'filled');
scatter3(yY,xY,0,800,'y','LineWidth',5);
%% close up
figure
zm = 7; izm = find(depth >= zm, 1, 'first');
zp = 200; izp = find(depth <= zp, 1, 'last');
iz = izm:izp;
hpatch = patch(isosurface(y(:,:,iz),x(:,:,iz),z(:,:,iz),model.v(:,:,iz),-1.6)); % 40 ohmm
isonormals(y(:,:,iz),x(:,:,iz),z(:,:,iz),model.v(:,:,iz),hpatch)
set(hpatch,'FaceColor','red','EdgeColor','none')
daspect([4,4,1])
%view([10,20])
view([-30,40])
axis tight
camlight right; 
set(gcf,'Renderer','zbuffer'); lighting phong
hold on;
data = readZ_3D('../WesternUS/Small_USA_5%_3%_run3_NLCG_027_rewrite.dat','[mV/km]/[nT]');
xdat = data{2}.siteLoc(:,1)./1000;
ydat = data{2}.siteLoc(:,2)./1000;
zdat = data{2}.siteLoc(:,3)./1000;
scatter3(ydat,xdat,zdat,'filled');
data = readZ_3D('Yellowstone_14freq_paper_errfl5T3.dat','[mV/km]/[nT]');
xdat = data{2}.siteLoc(:,1)./1000;
ydat = data{2}.siteLoc(:,2)./1000;
zdat = data{2}.siteLoc(:,3)./1000;
scatter3(ydat,xdat,zdat,'filled');
scatter3(yY,xY,0,800,'k','LineWidth',5);