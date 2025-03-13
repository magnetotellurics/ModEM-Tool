ll1 = [47.7 236];
ll2 = [47.7 239];
zLims = [0 80];
figure('Position',[100,100,800,300],'PaperPosition',[1,1,8,3],...
	'PaperOrientation','LandScape')
np = 60; nz = 30;
dz = 1/(nz-1);
dy = 1/(np-1);
y = 236+[0:dy:1]*3;
z = [0:dz:1]*80;
[XI,YI,ZI,C] = interpSlice(Res,ll1,ll2,zLims,np,nz);
pcolor(y,z,C');
shading interp
set(gca,'Ydir','reverse','fontweight','demi','fontsize',14);
colormap(cmap)
caxis([0,3]);
hcb = colorbar
set(hcb,'Fontweight','demi','fontsize',14)
