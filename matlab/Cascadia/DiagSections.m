load CascadiaInv

figure('Position',[100,100,900,500],'PaperPosition',[1,1,9,5], ...
'PaperOrientation','Landscape')

%   script for plotting diagonal sections
%    specify lat/lon endpoints
ll1 = [48,236; 46.6, 236; 49 236.6];
ll2 = [42,242; 42, 240; 43.4 242];
% latitude, longitude, depth limits   
LatLims = [42,50];
lonLims = [235,244];
zLims = [00,150];
%  resistivity limit
rhoLims = [.5,3.5];


[nDiag,dum] = size(ll1);
np = 60; nz = 40;

for k = 1:nDiag
   [XI,YI,ZI,C] = interpSlice(Res,ll1(k,:),ll2(k,:),zLims,np,nz);
   surface(YI,XI,ZI,C);
   caxis(rhoLims)
   shading flat
end

xLims(1) = lonLims(1)-.5;xLims(2) = lonLims(2)+.5
yLims(1) = latLims(1)-.5;yLims(2) = latLims(2)+.5
set(gca,'xlim',xLims,'ylim',yLims,'zlim',zLims,...
    'fontweight','demi','fontsize',14,'Zdir','reverse')
grid on

cmap = jet(64);
cmap = cmap(end:-1:1,:);
hold on
plotMap
%plotSlab


colormap(cmap)
hcb = colorbar
set(hcb,'FontWeight','demi','FontSize',16)
