load CascadiaInv

%   script for plotting fixed depth sections
%    specify depth sections to plot
%Depths = [20 40 60 80];
Depths = [30  65 100 135  ];
%Depths = [30];
%  resistivity limits
rhoLims = [.5,3.5];
%  latitude, longitude, depth limits (deg/deg/km)
latLims = [42,49];
lonLims = [236,244];
zLims = [0,150];

[ny,nx,nz] = size(Res.X);
z1 = Res.Z(1,1,:);
iY = fix(ny/2);
figure('Position',[100,100,900,900],'PaperPosition',[1,1,8,8])

nDepth = length(Depths);
cmap = jet(64);
cmap = cmap(end:-1:1,:);

temp = abs(Res.X(iY,:,1)-lonLims(1));
iX1 = find(min(temp)==temp);
temp = abs(Res.X(iY,:,1)-lonLims(2));
iX2 = find(min(temp)==temp);
temp = abs(Res.Y(:,1,1)-latLims(1));
iY1 = find(min(temp)==temp);
temp = abs(Res.Y(:,1,1)-latLims(2));
iY2 = find(min(temp)==temp);

for k = 1:nDepth
    %   find section nearest to Depths(k)
    temp = abs(z1-Depths(k));
    iZ = find(min(temp)==temp);
    x = squeeze(Res.X(iY1:iY2,iX1:iX2,iZ));
    y = squeeze(Res.Y(iY1:iY2,iX1:iX2,iZ));
    z = squeeze(Res.Z(iY1:iY2,iX1:iX2,iZ));
    rho = squeeze(Res.log10_rho(iY1:iY2,iX1:iX2,iZ));
    surface(x,y,z,rho)
    caxis(rhoLims);
    shading interp
    %hold on
end
zLims(2) = Res.z(iZ);
xLims(1) = lonLims(1)-.5;xLims(2) = lonLims(2)+.5
yLims(1) = latLims(1)-.5;yLims(2) = latLims(2)+.5
set(gca,'xlim',xLims,'ylim',yLims,'zlim',zLims,...
    'fontweight','demi','fontsize',14,'Zdir','reverse')
colormap(cmap)
hold on
plotMap
plotSlab
hcb = colorbar
%cbPos = get(hcb,'Position')
%cbPos(3) = cbPos(3)/2;
set(hcb,'FontWeight','demi','FontSize',16)%,'Position',cbPos);

