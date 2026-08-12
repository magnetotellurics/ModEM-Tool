mydir = '/Users/akelbert/Developer/ModelPlot/earthmodels/';
%dataFile = 'Small_USA_5%_3%_run3_NLCG_027_rewrite.dat';
dataFile = [mydir 'Small_USA_5%_3%_smoothed_data_run3_NLCG_027.dat'];
%dataFile = 'Yellowstone_10km_errfl5T3_smooth_NLCG_049.dat';
[header,units,isign,origin,info] = readZ_3Dinfo(dataFile);
[lat0,lon0] = latgrid.origin(info{1}.lat,info{1}.lon,info{1}.loc(:,1),...
                      info{1}.loc(:,2));

lims.depthmin = 10;
lims.depthmax = 250;
lims.latmin = 41;
lims.latmax = 46.5;
lims.lonmin = -119;
lims.lonmax = -107;
obj.limits = lims;
z = 3:3:400;
%% 
ROTATE = 0;
i=27;
%inversion = 'Yellowstone_10km_errfl5T3_200ohmm_NLCG';
%inversion = 'Yellowstone_10km_errfl5T3_smooth_NLCG';
inversion = 'Nested_NWUSA_Z_5%_Hz_3%_run3_NLCG';
modelFile = ['' inversion '_' sprintf('%03d',i) '.rho'];                 
while exist(modelFile,'file')
    [Cond] = readCond_3D(modelFile,2);
    obj = Cond_xy_ll(Cond,lat0,lon0);
    obj = CONDxyz2ll(obj,'depths',z);
    %[zcond,xy,ll] = slice(obj,sections{k}.Corners,sections{k}.NM);
    %
    %sections = zSlice(obj);
    %saveCONDll(obj,'YellowstoneMT',lims)
    % saveSlice
    depths = [lims.depthmin; obj.z(obj.z > lims.depthmin & obj.z <= lims.depthmax)];
    lons = [lims.lonmin, lims.lonmax];
    lats = [lims.latmin, lims.latmax];
    sections = zSlice(obj,depths,lons,lats);
    for k = 1:length(sections)
        sections{k}.NM = [64,102];
    end
clear sections
corners = zeros(3,2,2);

% ll1 = [41.0 241.0];
% ll2 = [46.5 253.5];
ll1 = [39.0 236.0];
ll2 = [49.5 256.5];
depths =  [100 100];

corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll2(1) ll1(2) depths(1)];
corners(:,1,2)  = [ll1(1) ll2(2) depths(2)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{1} = struct('Corners',corners,'NM',[100,125]);

ll1 = [46.0 246.0];
ll2 = [41.5 252.0];
depths =  [6 90];

corners(:,1,1)  = [ll1(1) ll1(2)+3 depths(1)];
corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
corners(:,1,2)  = [ll2(1) ll2(2)+3 depths(1)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{2} = struct('Corners',corners,'NM',[100,125]);

% ll1 = [41.0 241.0];
% ll2 = [46.5 253.5];
% depths =  [70 70];
% 
% corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
% corners(:,2,1)  = [ll2(1) ll1(2) depths(1)];
% corners(:,1,2)  = [ll1(1) ll2(2) depths(2)];
% corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
% sections{2} = struct('Corners',corners,'NM',[100,125]);

ll1 = [41.0 242.5];
ll2 = [46.5 253.5];
depths =  [6 500];

corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{3} = struct('Corners',corners,'NM',[100,125]);

ll1 = [39.0 236.0];
ll2 = [49.5 256.5];
depths =  [200 200];

corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll2(1) ll1(2) depths(1)];
corners(:,1,2)  = [ll1(1) ll2(2) depths(2)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{4} = struct('Corners',corners,'NM',[100,125]);

ll1 = [39.0 236.0];
ll2 = [49.5 256.5];
depths =  [300 300];

corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll2(1) ll1(2) depths(1)];
corners(:,1,2)  = [ll1(1) ll2(2) depths(2)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{5} = struct('Corners',corners,'NM',[100,125]);

% ll1 = [46.0 244.0];
% ll2 = [41.5 250.0];
% depths =  [6 70];
% 
% corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
% corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
% corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
% corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
% sections{4} = struct('Corners',corners,'NM',[100,125]);

fencePlot(obj,sections)
%colorbar
caxis([-3,0])
%daspect([1,1,2.7])
axis tight
box on
cb=colorbar;
%t=title(['NLCG iteration: ' num2str(i)],'FontWeight','demi','FontSize',24);
%set(t,'position',[-115 45 -80]);
title(cb,'\rho [{\Omega}m]','FontWeight','demi','FontSize',18)
set(cb,'ytick',[-3 -2.5 -2 -1.5 -1 -0.5 0])
set(cb,'yticklabel',[1000 300 100 30 10 3 1],'FontWeight','demi','FontSize',18)
set(cb,'position',get(cb,'position')+[0.05 0 0 0])
view([-40,15])
print('-dpng','-r300',[inversion '_' sprintf('%03d',i) '.png']);
if ROTATE
    for k=1:360
        view([-40+k,15])
        print('-dpng','-r300',[inversion '_' sprintf('%03d',i) '_' sprintf('%03d',k) '.png']);
    end
end
    fclose all; %close all;
    i = i+1;
    modelFile = ['models/' inversion '_' sprintf('%03d',i) '.rho'];
end
%%
obj = Cond_xy_ll(Cond,lat0,lon0);
k = find(obj.llgrid.depth>=300, 1, 'first');
figure; pcolor(obj.llgrid.lon,obj.llgrid.lat,squeeze(obj.CONDll(:,:,k)));
hold on; plot(info{1}.lon,info{1}.lat,'k^','MarkerSize',10,'LineWidth',2);
colorbar
caxis([-3,0])
%%
k = find(z==28);
[zcond,xy,ll] = slice(obj,sections{k}.Corners,sections{k}.NM);
figure; pcolor(squeeze(ll(2,:,:))-360,squeeze(ll(1,:,:)),zcond);
hold on; plot(info{1}.lon,info{1}.lat,'k^','MarkerSize',10,'LineWidth',2);
colorbar
caxis([-3,0])
%%
nslice = length(sections);
for k = 1:nslice
    [cond,xy,ll] = slice(obj,sections{k}.Corners,sections{k}.NM);
    surf(squeeze(ll(2,:,:))-360,squeeze(ll(1,:,:)),squeeze(ll(3,:,:)),cond);
    shading flat
    hold on
end
set(gca,'Fontweight','demi','FontSize',14,'zdir','reverse')
colorbar
caxis([-3,0])
hold on; plot(info{1}.lon,info{1}.lat,'k^','MarkerSize',10,'LineWidth',2);
hold off;

if ~exist(fname,'dir'); mkdir(fname); end
nslice = length(sections);
for k = 1:nslice
    [cond,xy,ll] = slice(obj,sections{k}.Corners,sections{k}.NM);
    str = sprintf('%06.2f',depths(k));
    f = fopen([fname '/' fname '_' str 'km.cond'],'w');
    for i = 1:sections{k}.NM(1)
        for j = 1:sections{k}.NM(2)
            fprintf(f,'%f\t%f\t%g\n',ll(1,i,j),ll(2,i,j)-360,cond(i,j));
        end
    end
    fclose(f);
end

%%
sections = zSlice(obj,[20 50]);
for k = 1:length(sections)
    sections{k}.NM = [100,125];
end
fencePlot(obj,sections);

%% horizonal slice for coords checking & comparison
corners = zeros(3,2,2);
ll1 = [41.0 241.0];
ll2 = [46.5 253.5];
depths =  [20 20];
corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll2(1) ll1(2) depths(1)];
corners(:,1,2)  = [ll1(1) ll2(2) depths(2)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
[zcond,xy,ll] = slice(obj,corners,[100 125]);
figure; pcolor(squeeze(ll(2,:,:))-360,squeeze(ll(1,:,:)),zcond);
hold on; plot(info{1}.lon,info{1}.lat,'k^','MarkerSize',10,'LineWidth',2);

%%   test plotting of a few slices
latitudes = [40,41,42,43,44,45,46,47,48];
depths = [5,100];
sections = latSlice(obj,latitudes,depths);
for k = 1:length(sections)
    sections{k}.NM = [100,125];
end
fencePlot(obj,sections)

%%
clear sections
corners = zeros(3,2,2);

ll1 = [41.0 241.0];
ll2 = [46.5 253.5];
depths =  [20 20];

corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll2(1) ll1(2) depths(1)];
corners(:,1,2)  = [ll1(1) ll2(2) depths(2)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{1} = struct('Corners',corners,'NM',[100,125]);

ll1 = [46.0 246.0];
ll2 = [41.5 252.0];
depths =  [6 90];

corners(:,1,1)  = [ll1(1) ll1(2)+3 depths(1)];
corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
corners(:,1,2)  = [ll2(1) ll2(2)+3 depths(1)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{2} = struct('Corners',corners,'NM',[100,125]);

% ll1 = [41.0 241.0];
% ll2 = [46.5 253.5];
% depths =  [70 70];
% 
% corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
% corners(:,2,1)  = [ll2(1) ll1(2) depths(1)];
% corners(:,1,2)  = [ll1(1) ll2(2) depths(2)];
% corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
% sections{2} = struct('Corners',corners,'NM',[100,125]);

ll1 = [41.0 242.5];
ll2 = [46.5 253.5];
depths =  [6 90];

corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
sections{3} = struct('Corners',corners,'NM',[100,125]);

% ll1 = [46.0 244.0];
% ll2 = [41.5 250.0];
% depths =  [6 70];
% 
% corners(:,1,1)  = [ll1(1) ll1(2) depths(1)];
% corners(:,2,1)  = [ll1(1) ll1(2) depths(2)];
% corners(:,1,2)  = [ll2(1) ll2(2) depths(1)];
% corners(:,2,2)  = [ll2(1) ll2(2) depths(2)];
% sections{4} = struct('Corners',corners,'NM',[100,125]);

fencePlot(obj,sections)
%colorbar
%caxis([-3,0])

cb=colorbar
title(cb,'\rho [{\Omega}m]','FontWeight','demi','FontSize',18)
set(cb,'ytick',[-3 -2.5 -2 -1.5 -1 -0.5 0])
set(cb,'yticklabel',[1000 300 100 30 10 3 1],'FontWeight','demi','FontSize',18)


%hold on; plot(info{1}.lon+360,info{1}.lat,'k^','MarkerSize',10,'LineWidth',2);

