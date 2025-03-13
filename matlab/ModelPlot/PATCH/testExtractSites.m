%%   This is an example script ... loads data file, patch file, and then
%    extracts all sites that fall within the patch, and makes a little plot
%    of phases (in the present version both observed and predicted)
%    for off-diagonal components.  Uses TTrFunZ (impedance)
%    transfer function class; could easily modify this to do other things
clear

PRINT = false;
iPatch = 1;
load EBRpat.mat
dataFile = 'Z_Hz_5%_3%_22Per_final_Data_set.dat';
predFile = 'Small_USA_5%_3%_run3_NLCG_027_rewrite.dat';
[allData,header,units,isign,origin,info] = readZ_3D(dataFile);
[allDataPred,header,units,isign,origin,infoPred] = readZ_3D(predFile);

% first find which sites are in the area covered by the patch
%   this uses function  whichPatch_ll
ISTA = [];
lon = 360+info{1}.lon;
lat = info{1}.lat;
nSites = 0;
latsUse = [];
lonsUse = latsUse;
for l = 1:length(lat)
    [k,name] = whichPatch_ll(lat(l),lon(l),patchStruct);
    if k == iPatch
        ISTA = [ISTA l];
        nSites = nSites+1;
        sites{nSites} = info{1}.code(l,:);
        latsUse = [latsUse lat(l)];
        lonsUse = [lonsUse lon(l)];
    end
end
ZTF = cell(nSites,1);
ZTFpred =  ZTF;
DIR = ['patch' num2str(iPatch)];
if ~exist(DIR,'dir')
    system(['mkdir ' DIR ])
end
nPeriodMin = 10;
%%
for l = 1:nSites
    ZTF{l} = TTrFunZ;
    ZTFpred{l} = TTrFunZ;
    ModEMdata2TTrFunZ(ZTF{l},info{1},ISTA(l));
    ModEMdata2TTrFunZ(ZTFpred{l},infoPred{1},ISTA(l));
    nPeriod = length(ZTF{l}.T);
    
    if nPeriod >=10
        ZTF{l}.ap_res;
        ZTFpred{l}.ap_res;
        %ZTF{l}.phasePlot;
        
        
        %ZTF{l}.phasePlot(ZTFpred{l})
        if PRINT
            cfile = [ DIR '/'  ZTF{l}.Header.Sites{1} '_phase.jpg'];
            eval(['print -djpeg90 ' cfile])
        end
    end
end

%% put a bunch of little rho plots on one page
siteSpacing = 70;
[N,M,rects] = arrayMultiPlotSet(latsUse,lonsUse,siteSpacing);
screenRect = [100,100 N*300 M*200];
paperRect = [1,1,N*2,M*1.5];
hfig = figure('Position',screenRect,'PaperPosition',paperRect);
for k = 1:nSites
    if length(ZTF{k}.T) >=10
        hax = phaseSubPlot(ZTF{k},hfig,rects(k,:),ZTFpred{k});
        set(get(hax,'Xlabel'),'string',[]);
        set(get(hax,'Ylabel'),'string',[]);
    end
end

%% plot all phases for the whole array
nSites = length(lat);
ZTF = cell(nSites,1);
ZTFpred = ZTF;
for l = 1:nSites
    ZTF{l} = TTrFunZ;
    ZTFpred{l} = TTRfunZ;
    ModEMdata2TTrFunZ(ZTF{l},info{1},l);
    ModEMdata2TTrFunZ(ZTFpred{l},infoPred{1},l);
    ZTF{l}.ap_res;
    ZTFpred{l}.ap_res;
    nPeriod = length(ZTF{l}.T);
end

%%
siteSpacing = 70;
[N,M,rects] = arrayMultiPlotSet(lat,lon,siteSpacing);
screenRect = [100,100 N*300 M*200];
paperRect = [1,1,N*2,M*1.5];
hfig = figure('Position',screenRect,'PaperPosition',paperRect);
for k = 1:nSites
    if length(ZTF{k}.T) >=10
        hax = phaseSubPlot(ZTF{k},hfig,rects(k,:),ZTFpred{k});
        set(get(hax,'Xlabel'),'string',[]);
        set(get(hax,'Ylabel'),'string',[]);
    end
end
