dataFile = 'Yellowstone_14freq_paper_errfl5T3.dat';
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
i=0;
inversion = 'Yellowstone_10km_errfl5T3_200ohmm_smooth_NLCG';
%inversion = 'Yellowstone_10km_200ohmm_NLCG';
%inversion = 'Yellowstone_10km_errfl5T3_smooth_NLCG';
%inversion = 'Yellowstone_10km_final_radial_NLCG';
modelFile = ['' inversion '_' sprintf('%03d',i) '.rho'];                 
while exist(modelFile,'file')
    [Cond] = readCond_3D(modelFile,2);
    if findstr(Cond.paramType,'LINEAR')
        Cond.v = log10(Cond.v);
        Cond.paramType = 'LOG10';
    end
    obj = Cond_xy_ll(Cond,lat0,lon0);
    
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
caxis([-3,0])
daspect([1,1,15])
axis tight
box on
cb=colorbar;
if ~ ROTATE
    t=title(['NLCG iteration: ' num2str(i)],'FontWeight','demi','FontSize',24);
    %set(t,'position',[-115 45 -80]);
end
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
    fclose all; close all;
    i = i+1;
    modelFile = ['' inversion '_' sprintf('%03d',i) '.rho'];
end
