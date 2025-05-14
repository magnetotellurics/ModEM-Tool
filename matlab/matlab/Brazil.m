%% read files
modelFile = 'Mauricio_inv_final/Modular_NLCG_086.rho';
predFile = 'Mauricio_inv_final/Modular_NLCG_086.dat';
dataFile = 'Mauricio_inv_final/data-imp-gds_3300s1_3p-NG-ok.dat';

clear obj
obj = xymodel.read(modelFile,'WS');
uiplot(obj,10);

%% convert to regular lat/lon
lat0 = -10.8815;
lon0 = -39.3203; 
obj = setOrigin(obj,lat0,lon0);

llobj = llmodel(obj);

%[header,units,isign,origin,info] = readZ_3Dinfo(dataFile);
%[lat0,lon0] = latlontools.latlon0(info{1}.lat,info{1}.lon,info{1}.loc(:,1),...
%                      info{1}.loc(:,2));
%sitelatlon = latlontools.xy2ll(info{1}.loc'/1000,lat0,lon0);
pred = readZ_3D(predFile,'[mV/km]/[nT]');
data = readZ_3D(dataFile,'[mV/km]/[nT]');
sitelatlon = latlontools.xy2ll(data{2}.siteLoc'/1000,lat0,lon0);
sitedepth = data{2}.siteLoc(:,3)';

%% explore RMS
[rms,info] = DataFit(predFile,dataFile,5);

total_rms = sqrt( nansum(nansum(nansum(info.res))) / sum(sum(sum(~isnan(info.res)))) );

sitenum = str2num(info.code);
line{1} = find((sitenum>=0 & sitenum<=12) | sitenum>=61);
line{2} = find(sitenum>=13 & sitenum<=25);
line{3} = [find(sitenum==007); find(sitenum>=26 & sitenum<=30); find(sitenum==23)];
line{4} = find(sitenum>=31 & sitenum<=36);
line{5} = find(sitenum>=37 & sitenum<=47);
line{6} = flipud(find(sitenum>=48 & sitenum<=60));

l1_rms = sqrt( nansum(nansum(nansum(info.res(line{1},:,:)))) / sum(sum(sum(~isnan(info.res(line{1},:,:))))) );
l2_rms = sqrt( nansum(nansum(nansum(info.res(line{2},:,:)))) / sum(sum(sum(~isnan(info.res(line{2},:,:))))) );
l3_rms = sqrt( nansum(nansum(nansum(info.res(line{3},:,:)))) / sum(sum(sum(~isnan(info.res(line{3},:,:))))) );
l4_rms = sqrt( nansum(nansum(nansum(info.res(line{4},:,:)))) / sum(sum(sum(~isnan(info.res(line{4},:,:))))) );
l5_rms = sqrt( nansum(nansum(nansum(info.res(line{5},:,:)))) / sum(sum(sum(~isnan(info.res(line{5},:,:))))) );
l6_rms = sqrt( nansum(nansum(nansum(info.res(line{6},:,:)))) / sum(sum(sum(~isnan(info.res(line{6},:,:))))) );

ZXX_rms = sqrt( nansum(nansum(info.res(:,:,1))) / sum(sum(~isnan(info.res(:,:,1)))) );
ZXY_rms = sqrt( nansum(nansum(info.res(:,:,2))) / sum(sum(~isnan(info.res(:,:,2)))) );
ZYX_rms = sqrt( nansum(nansum(info.res(:,:,3))) / sum(sum(~isnan(info.res(:,:,3)))) );
ZYY_rms = sqrt( nansum(nansum(info.res(:,:,4))) / sum(sum(~isnan(info.res(:,:,4)))) );
Z_rms = sqrt( nansum(nansum(nansum(info.res(:,:,1:4)))) / sum(sum(sum(~isnan(info.res(:,:,1:4))))) );
TX_rms = sqrt( nansum(nansum(info.res(:,:,5))) / sum(sum(~isnan(info.res(:,:,5)))) );
TY_rms = sqrt( nansum(nansum(info.res(:,:,6))) / sum(sum(~isnan(info.res(:,:,6)))) );
T_rms = sqrt( nansum(nansum(nansum(info.res(:,:,5:6)))) / sum(sum(sum(~isnan(info.res(:,:,5:6))))) );

disp( ' Subset  |  RMS'); 
disp( '------------------'); 
disp([' ALL     |  ' num2str(total_rms)]); 
disp([' 000-012 & 061-067 |  ' num2str(l1_rms)]); 
disp([' 013-025 |  ' num2str(l2_rms)]); 
disp([' 026-030 & 023 |  ' num2str(l3_rms)]); 
disp([' 031-036 |  ' num2str(l4_rms)]); 
disp([' 037-047 |  ' num2str(l5_rms)]); 
disp([' 048-060 |  ' num2str(l6_rms)]); 
disp([' ZXX     |  ' num2str(ZXX_rms)]); 
disp([' ZXY     |  ' num2str(ZXY_rms)]); 
disp([' ZYX     |  ' num2str(ZYX_rms)]); 
disp([' ZYY     |  ' num2str(ZYY_rms)]); 
disp([' ALL Z   |  ' num2str(Z_rms)]); 
disp(['  TX     |  ' num2str(TX_rms)]); 
disp(['  TY     |  ' num2str(TY_rms)]); 
disp([' ALL T   |  ' num2str(T_rms)]); 

for i = 1:length(info.code)
    for j = 1:6
    	rms(i,j) = sqrt( nansum(nansum(info.res(i,:,j))) / sum(sum(~isnan(info.res(i,:,j)))) );
    end
    temp = sprintf('%s\t%f %f %f %f %f %f',info.code(i,:),rms(i,1),rms(i,2),rms(i,3),rms(i,4),rms(i,5),rms(i,6));
    disp(temp)
end

lnum = 6;
len = length(line{lnum});
figure('Position',[300,300,1200,400],...
        'PaperPosition',[1,1,12,4]);
h(1)=plot(1:len,rms(line{lnum},1),'linewidth',3,'color','b'); hold on;
h(2)=plot(1:len,rms(line{lnum},2),'linewidth',3,'color','r'); hold on;
h(3)=plot(1:len,rms(line{lnum},3),'linewidth',3,'color','k'); hold on;
h(4)=plot(1:len,rms(line{lnum},4),'linewidth',3,'color','g'); hold on;
h(5)=plot(1:len,rms(line{lnum},5),'linewidth',3,'color','c'); hold on;
h(6)=plot(1:len,rms(line{lnum},6),'linewidth',3,'color','m'); hold off;
set(gca,'xtick',1:len,'xticklabel',info.code(line{lnum},:))
set(gca,'fontweight','demi','fontsize',14)
ylabel('RMS')
legend(h,'Zxx','Zxy','Zyx','Zyy','Tx','Ty','Location','Best');
print('-djpeg95','-r300',['line' num2str(lnum) '_rms_by_site.jpg']);


%% plotting data as colored circles
% 
% ampmax=4.0;
% figure(1);clf
% pcolor(LON,LAT,amp');shading flat;
% colormap('jet');caxis([0, ampmax]);
% map=colormap;[mc,dum]=size(map);
% %
% tcol='k';
% for k=1:length(tgamp)
%  ik=max(1,round(tgamp(k)/ampmax*mc));
%  ik=min(ik,mc);
%  col=map(ik,:);
%  pt=plot(tlon(k),tlat(k),'o','MarkerFaceColor',col,...
%                           'MarkerEdgeColor','k','MarkerSize',12);
% end 

%% make covariances and synthetic models
COV = mask(obj.grid,'layer',[60 1000],9);
uiplot(obj.grid,COV,10);
writeCov_3D('Brazil_freeze_below_60km.cov',COV);

clear sizes layers indices
sizes = 3*2.^(0:4);
layers(:,1) = 2*[0; cumsum(sizes(1:end-1))'];
layers(:,2) = layers(:,1) + 2*sizes';
indices(:,1) = (10:14);
indices(:,2) = (20:24);
obj.grid.xpadding = 16;
obj.grid.ypadding = 16;
COV = mask(obj.grid,'checkerboard',layers,sizes,indices);
uiplot(obj.grid,COV);

halfspace = -log10(200);
modobj = obj;
ii = COV >= 20 & ~mod(COV,2);
modobj.v(ii) = -log10(10000);
ii = COV >= 20 & mod(COV,2);
modobj.v(ii) = -log10(1000);
ii = COV < 20 & ~mod(COV,2);
modobj.v(ii) = -log10(1);
ii = COV < 20 & mod(COV,2);
modobj.v(ii) = -log10(10);
ii = COV == 1;
modobj.v(ii) = halfspace;
uiplot(modobj,10);
write('Brazil_checkerboard_with_distortion.ws',modobj,'WS');
ii = ~mod(COV,10);
modobj.v(ii) = halfspace;
uiplot(modobj,10);
write('Brazil_checkerboard_without_distortion.ws',modobj,'WS');

COV = mask(obj.grid,'layer',[60 1000],2);
ii = COV == 2;
halfspace = mean(mean(obj.v(ii)));
modobj = obj;
modobj.v(ii) = halfspace;
uiplot(modobj,10);
write('Brazil_inverse_averaged_below_60km.ws',modobj,'WS');

prior = xymodel.prior(grid,'uniform',200);
write('Brazil_3km_200ohmm_pr.ws',prior,'WS');

%% plot a horizontal slice for testing
% k = find(llobj.grid.depth>28,1);
% pcolor(llobj.grid.lon,llobj.grid.lat,llobj.v(:,:,k));
% shading flat; colorbar

sections = zSlice(obj,16);
[zcond,xy,ll] = slice(obj,sections{1}.Corners,sections{1}.NM);
figure; pcolor(squeeze(ll(2,:,:)),squeeze(ll(1,:,:)),zcond); shading flat;
hold on; plot(sitelatlon(2,:),sitelatlon(1,:),'k^','MarkerSize',10,'LineWidth',2);
colorbar
caxis([-3,0])

%% isoplot geographic coords large domain
figure
padding = 10;
zm = 0; izm = find(llobj.grid.depth >= zm, 1, 'first');
zp = 200; izp = find(llobj.grid.depth <= zp, 1, 'last');
iz = izm:izp;
ym = -9.5; iym = padding; %find(lon(nlat/2,:,1) >= ym, 1, 'first');
yp = -12.; iyp = llobj.grid.nlon-padding; %find(lon(nlat/2,:,1) <= yp, 1, 'last');
iy = iym:iyp;
xm = -41; ixm = padding; %find(lat(:,1,1) >= xm, 1, 'first');
xp = -38; ixp = llobj.grid.nlat-padding; %find(lat(:,1,1) <= xp, 1, 'last');
ix = ixm:ixp;
[LON,LAT,Z] = meshgrid(llobj.grid.lon,llobj.grid.lat,llobj.grid.depth);
hpatch = patch(isosurface(LON(ix,iy,iz),LAT(ix,iy,iz),-Z(ix,iy,iz),llobj.v(ix,iy,iz),-1.48));
isonormals(LON(ix,iy,iz),LAT(ix,iy,iz),-Z(ix,iy,iz),llobj.v(ix,iy,iz),hpatch)
set(hpatch,'FaceColor','red','EdgeColor','none')
hold on
%hpatch2 = patch(isosurface(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),-model.v(ix,iy,iz),2.47));
%isonormals(lon(ix,iy,iz),lat(ix,iy,iz),z(ix,iy,iz),model.v(ix,iy,iz),hpatch2)
%set(hpatch2,'FaceColor','blue','EdgeColor','none')
%daspect([1,1,40])
view([7,29])
axis tight
box on
camlight left; 
set(gcf,'Renderer','zbuffer'); lighting phong
set(gca,'fontsize',16,'fontweight','demi');
hold on;
scatter3(sitelatlon(2,:),sitelatlon(1,:),sitedepth,'filled');
scatter3(sitelatlon(2,:),sitelatlon(1,:),sitedepth,'filled');
%scatter3(-110.3,44.6,0,400,'k','LineWidth',5);

%% Anomalia norte
clear sections
corners = zeros(3,2,2);

ll1 = [-10.8 -38.4];
ll2 = [-10.05 -39.65];
depths =  [0 60];
corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{1} = struct('Corners',corners,'NM',[100,125]);


ll1 = [-11.5 -38.3];
ll2 = [-9.8 -39.7];
depths =  [0 60];
corners(:,1,1)  = [ll1(1) ll1(2) 10];
corners(:,2,1)  = [ll2(1) ll1(2) 10];
corners(:,1,2)  = [ll1(1) ll2(2) 10];
corners(:,2,2)  = [ll2(1) ll2(2) 10];
sections{2} = struct('Corners',corners,'NM',[100,125]);


%ll1 = [-10.8 -38];
%ll2 = [-10.8 -41];
%depths =  [0 60];
%corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
%corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
%corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
%corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
%sections{2} = struct('Corners',corners,'NM',[100,125]);



%ll1 = [-10.8 -38.5];
%ll2 = [-10.8 -40.0];
%depths =  [0 60];
%corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
%corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
%corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
%corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
%sections{3} = struct('Corners',corners,'NM',[100,125]);


fencePlot(obj,sections)
colorbar
caxis([-3,0])


%% Tomo - thrree depth slices (2, 10 e 30)
clear sections
corners = zeros(3,2,2);

%ll1 = [-10.8 -38.4];
%ll2 = [-10.8 -39.7];
ll1 = [-12.25 -38.3];
ll2 = [-9.5 -39];
depths =  [0 60];
%corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
%corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
%corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
%corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
corners(:,1,1)  = [ll1(1) ll1(2) 2];
corners(:,2,1)  = [ll2(1) ll1(2) 2];
corners(:,1,2)  = [ll1(1) ll2(2) 2];
corners(:,2,2)  = [ll2(1) ll2(2) 2];
sections{1} = struct('Corners',corners,'NM',[100,125]);


ll1 = [-12 -38.3];
ll2 = [-10 -40.35];
depths =  [0 60];
corners(:,1,1)  = [ll1(1) ll1(2) 10];
corners(:,2,1)  = [ll2(1) ll1(2) 10];
corners(:,1,2)  = [ll1(1) ll2(2) 10];
corners(:,2,2)  = [ll2(1) ll2(2) 10];
sections{2} = struct('Corners',corners,'NM',[100,125]);


ll1 = [-12 -38];
ll2 = [-10 -40.35];
depths =  [0 60];
corners(:,1,1)  = [ll1(1) ll1(2) 30];
corners(:,2,1)  = [ll2(1) ll1(2) 30];
corners(:,1,2)  = [ll1(1) ll2(2) 30];
corners(:,2,2)  = [ll2(1) ll2(2) 30];
sections{3} = struct('Corners',corners,'NM',[100,125]);

%ll1 = [-10.8 -38];
%ll2 = [-10.8 -41];
%depths =  [0 60];
%corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
%corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
%corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
%corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
%sections{2} = struct('Corners',corners,'NM',[100,125]);



%ll1 = [-10.8 -38.5];
%ll2 = [-10.8 -40.0];
%depths =  [0 60];
%corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
%corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
%corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
%corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
%sections{3} = struct('Corners',corners,'NM',[100,125]);


fencePlot(obj,sections)
colorbar
caxis([-3,0])


%%  vertical sections
clear sections
corners = zeros(3,2,2);


ll1 = [-10.15 -38.3];
ll2 = [-10.15 -40.35];
depths =  [5 40];
corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{1} = struct('Corners',corners,'NM',[100,125]);

ll1 = [-10.5 -38.3];
ll2 = [-10.5 -40.35];
depths =  [0 40];
corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{4} = struct('Corners',corners,'NM',[100,125]);


ll1 = [-10.8 -38.3];
ll2 = [-10.8 -40.35];
depths =  [0 40];
corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{2} = struct('Corners',corners,'NM',[100,125]);

ll1 = [-10.85 -38.6];
ll2 = [-11.75 -40.2];
depths =  [0 40];
corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{3} = struct('Corners',corners,'NM',[100,125]);

%depth slicE
%ll1 = [-12 -38.3];
%ll2 = [-10 -40.35];
%depths =  [5 40];
%corners(:,1,1)  = [ll1(1) ll1(2) 10];
%corners(:,2,1)  = [ll2(1) ll1(2) 10];
%corners(:,1,2)  = [ll1(1) ll2(2) 10];
%corners(:,2,2)  = [ll2(1) ll2(2) 10];
%sections{4} = struct('Corners',corners,'NM',[100,125]);


fencePlot(obj,sections)
colorbar
caxis([-3,0])
cb=colorbar;
title(cb,'\rho [{\Omega}m]','FontWeight','demi','FontSize',18)
set(cb,'ytick',[-3 -2.5 -2 -1.5 -1 -0.5 0])
set(cb,'yticklabel',[1000 300 100 30 10 3 1],'FontWeight','demi','FontSize',18)
%hold on; plot(info{1}.lon+360,info{1}.lat,'k^','MarkerSize',10,'LineWidth',2);

