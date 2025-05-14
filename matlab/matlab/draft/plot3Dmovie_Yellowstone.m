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
    
filt = [modeldir '/*.ws; *.cpr; *.rho; *.prm'];
[filename, pathname] = uigetfile(filt, 'Resistivity Model File');

%  read in chosen file
cfile = [pathname filename];

Cond = readCond_3D(cfile,2);

% convert to log10
if strcmp(Cond.paramType,'LOGE')
    Cond.paramType = 'LOG10';
    Cond.v = Cond.v / log(10);
    Cond.AirCond = Cond.AirCond / log(10);
elseif strcmp(Cond.paramType,'LINEAR')
    Cond.paramType = 'LOG10';
    Cond.v = log10(Cond.v);
    Cond.AirCond = log10(Cond.AirCond);
end

filt = './*.dat; *.res; *.imp';
[filename, pathname] = uigetfile(filt, 'Load Data File');
Zfile = [pathname filename];

[AllData] = readZ_3D(Zfile,'[mV/km]/[nT]');

plot3Dview

print('-djpeg','-r300', cfile(1:end-3));
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

