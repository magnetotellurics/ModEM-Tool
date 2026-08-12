%%   This is an example script that computes (anistropic) average
%        conductivity as a function of depth for a given patch
%    To run: need to specify patch file, model file (which must be
%    consistent, of course), data file (to set lat/lon) + patch number

clear
iPatch = 1;
%  mindiff is minimum relative difference between conductive and resistive
%      directions to plot azimuth
mindiff = .25;
%   maximum number of vertical layers to plot
izMax=36;
%  set PRINT = true to print plot file
PRINT = 0;

%   load patch file
load PAPERpat.mat
nPatch = length(patchStruct.names);
result = cell(nPatch,1);
%   load model file
modelFile = 'Nested_NWUSA_Z_5%_Hz_3%_run3_NLCG_027.rho';
%   also need  data file
dataFile = 'Z_Hz_5%_3%_22Per_final_Data_set.dat';
[~,~,~,~,~,info] = readZ_3D(dataFile);
[lat0,lon0] = latgrid.origin(info{1}.lat,info{1}.lon,info{1}.loc(:,1),...
                      info{1}.loc(:,2));
[Cond] = readCond_3D(modelFile,2); 
obj = Cond_xy_ll(Cond,lat0,lon0);
for iPatch = 1:nPatch
%   call depthProfile to get structure containing conductivity vs. depth +
%      various intermediate results
result{iPatch} = depthProfile(obj,patchStruct,iPatch);

%   plot max/min conductivity vs. depth, azimuth of max conductivity 
figure('Position',[100,100,550,500],'PaperPosition',[1,1,5.5,5]);
axes('Position',[ 0.15    0.1100    0.4    0.8150])
semilogx(10.^result{iPatch}.sigmaMax(1:izMax),result{iPatch}.z(1:izMax),'r')
hold on
axis('ij')
set(gca,'FontWeight','demi','FontSize',14,'xlim',10.^[-3.5,-.5],'Xtick',[10^-3,10^-2,10^-1])
plot(10.^result{iPatch}.sigmaMin(1:izMax),result{iPatch}.z(1:izMax),'b')
plot(10.^((result{iPatch}.sigmaMax(1:izMax)+result{iPatch}.sigmaMin(1:izMax))/2),result{iPatch}.z(1:izMax),'k')
fatlines(gca,2)
yl = get(gca,'ylim');
xlabel('Conductivity (S/m)')
ylabel('depth (km)')
title(patchStruct.names(iPatch),'FontSize',16)
%    plot significant azimuths
axes('Position',[ 0.6    0.11    0.3    0.8150])
indPlot = (result{iPatch}.sigmaMax-result{iPatch}.sigmaMin)*log(10) >mindiff;
plot(result{iPatch}.thetaMax(indPlot),result{iPatch}.z(indPlot),'k^','markersize',10,...
    'linewidth',2)
axis('ij')
set(gca,'FontWeight','demi','FontSize',14,'xlim',[-90,90],...
    'Xtick',[-90,-45,0,45,90],'yticklabel',[],'ylim',yl)
xlabel('degree (E = 0)')
title('Conductive Azimuth')
if PRINT
    cfile = [DIR '/CondAvgAnis.eps'];
    eval(['print -depsc ' cfile]);
end
end
