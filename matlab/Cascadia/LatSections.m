load CascadiaInv

figure('Position',[100,100,900,900],'PaperPosition',[1,1,8,8])
%figure('Position',[100,100,900,500],'PaperPosition',[1,1,8,4])

%   script for plotting fixed latitude sections
%    specify list of latitudes sections to plot
%sectionLats = [42 43,44,45,46,47, 48,49];
sectionLats = [45.5 46 46.5 47 47.5];
%  resistivity limits
rhoLims = [0.5,3.5];
%  latitude, longitude, depth limits (deg/deg/km)
latLims = [45,48];
lonLims = [236.5,239];
zLims = [0,100];

[ny,nx,nz] = size(Res.X);
nx2 = fix(nx/2);
yCenter = squeeze(Res.Y(:,nx2,1));

nSec = length(sectionLats);
cmap = jet(64);
cmap = cmap(end:-1:1,:);

for k = 1:nSec
    %   find section nearest to sectionLat
    temp = abs(yCenter-sectionLats(k));
    iY = find(min(temp)==temp);
    temp = abs(Res.X(iY,:,1)-lonLims(1));
    iX1 = find(min(temp)==temp);
    temp = abs(Res.X(iY,:,1)-lonLims(2));
    iX2 = find(min(temp)==temp);
    temp = abs(Res.Z(1,1,:)-zLims(1));
    iZ1 = find(min(temp)==temp);
    temp = abs(Res.Z(1,1,:)-zLims(2));
    iZ2 = find(min(temp)==temp);
    x = squeeze(Res.X(iY,iX1:iX2,iZ1:iZ2));
    y = squeeze(Res.Y(iY,iX1:iX2,iZ1:iZ2));
    z = squeeze(Res.Z(iY,iX1:iX2,iZ1:iZ2));
    rho = squeeze(Res.log10_rho(iY,iX1:iX2,iZ1:iZ2));
    surface(x,y,z,rho)
    caxis(rhoLims);
    shading interp
    %hold on
end
zLims(2) = Res.z(iZ2);
xLims(1) = lonLims(1)-.5;xLims(2) = lonLims(2)+.5
yLims(1) = latLims(1)-.5;yLims(2) = latLims(2)+.5
set(gca,'xlim',xLims,'ylim',yLims,'zlim',zLims,...
    'fontweight','demi','fontsize',14,'Zdir','reverse')

set(gca,'xlim',lonLims,'ylim',latLims,'zlim',zLims,...
    'fontweight','demi','fontsize',14,'Zdir','reverse')

colormap(cmap)
hold on
plotMap
plotSlab
%colorbar
grid on
