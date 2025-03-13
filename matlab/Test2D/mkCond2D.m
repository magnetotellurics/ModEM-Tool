function [model,data,prior] = mkCond2D(conf,modelFile,dataFile)

% Usage: [model,data,prior] = mkCond2D(conf,modelFile,dataFile)
%
% Inputs are file names and configuration structure created by
% Config2D. Outputs are matlab structures; prior is the optional
% prior model structure.
% If dataFile is not specified, assume that no data are required.

if nargin == 1
    WRITE_FILES = 0;
else
    WRITE_FILES = 1;
end

% read the conf structure
LOGCOND = conf.LOGCOND;
COMPUTE_DATA = conf.COMPUTE_DATA;
ADD_NOISE = conf.ADD_NOISE;
MODE = conf.mode;
FracError = conf.data_error;
NODE_SITE = conf.NODE_SITE;
PLOT_MODEL = conf.PLOT_MODEL;
PLOT_DATA = conf.PLOT_DATA;

if nargin < 3 && COMPUTE_DATA
    display('Unable to compute data: no data file specified')
    COMPUTE_DATA = 0;
end

%  strip all extensions off file names
if WRITE_FILES
    [modelPath, modelBaseName] = fileparts(modelFile);
    if COMPUTE_DATA
        [dataPath, dataBaseName, ext] = fileparts(dataFile);
        dataTemplate = fullfile(dataPath,['Template_' dataBaseName ext]);
    end
else
    modelPath = pwd;
    modelBaseName = 'model';
end

%    compute horizontal grid from input parameters in conf
DyCenter = conf.dyBase*ones(1,conf.Ny);
DyPad = conf.dyBase*(conf.yPadFac.^[1:conf.NyPad]);
Dy = [DyPad(end:-1:1) DyCenter DyPad];
ny = length(Dy);

%    compute vertical grid from input parameters in conf
Dzb  = conf.dzBase*(conf.zFac.^[0:conf.Nz-1]);
nzb = length(Dzb);

%    NOTE:   Air layers are hard-coded in ModEM2D ... leave this alone!
nza = 10;
Dza(1) = 10.0;
for k = 2:nza
   Dza(k) = Dza(k-1)*3;
end
conf.NzAir = nza;
%Dza = conf.dzBase*(conf.zAirFac.^[0:conf.NzAir-1]);
%nza = length(Dza);
Dz = [Dza(end:-1:1) Dzb];
nz = nza+nzb;

%  create grid 
gridDef = struct('Nz',nz,'Nza',nza,'Ny',ny,...
		'Dy',Dy,'Dz',Dz);
grid = gridDef;

%   set up background
airCond = 1e-10;
Cond = ones(ny,nzb)*conf.bg_cond;

%   insert blocks; now using a loop, allowing any number of blocks

for k = 1:conf.nBlocks
    Iy1 = (conf.y0(k):conf.y1(k));  
    Iz1 = (conf.z0(k):conf.z1(k));
    Cond(Iy1,Iz1) = conf.anomaly(k);
end

%  write out conductivity parameter
earthCond = struct('paramType','LINEAR',...
	'v',Cond,'AirCond',airCond,'grid',grid);
logEarthCond = earthCond;
logEarthCond.v = log(earthCond.v);
logEarthCond.AirCond = log(earthCond.AirCond);
logEarthCond.paramType = 'LOGE';
if LOGCOND
    model = logEarthCond;
    prior = model;
    prior.v(:,:) = log(conf.bg_cond);
else
    model = earthCond;
    prior = model;
    prior.v(:,:) = conf.bg_cond;
end
if WRITE_FILES
    [status] = writeCond_2D(modelFile,model,grid);
end



y = cumsum(Dy);
z = cumsum(Dz);


%  write out data file (might be a template only ... 
%          with no actual data values)
%  data will be stored (in matlab) in a cell array: one cell
%   for each period; each cell contains a structure containing
%   period, data sites, data, and error bars
%  For the test nSites are uniformly spaced centered over anomaly
if NODE_SITE
   nNodes = 8;
   yNodes = y(20+(1:nNodes)*7);
   siteLoc = zeros(1,2*nNodes-1);
   siteLoc(1:2:end) = yNodes;
   siteLoc(2:2:end-1) = (yNodes(1:end-1)+yNodes(2:end))/2;
   siteLoc = siteLoc';
   nSites = length(siteLoc);
   siteChar = sprintf('%03d\n',1:nSites);
   siteChar = regexp(siteChar(1:end-1),'\n','split');
else
   nSites = conf.nSites;
   y1 = y(conf.FirstSite); y2 = y(conf.LastSite);   
   ys = (y2-y1)/(nSites-1);
   siteLoc = (y1:ys:y2)';
   siteChar = sprintf('%03d\n',1:nSites);
   siteChar = regexp(siteChar(1:end-1),'\n','split');
end

if PLOT_MODEL
    mPlot = logEarthCond;
    nZplot = 25;
    nYskip = conf.NyPad;
    cax = [-3.5,-.5];  %   +log10(2);  Model 2
    ctitle = ['Synthetic Model : ' modelBaseName]; 
    OPTIONS = struct('nYskip',nYskip,'nZplot',nZplot,...
        'siteLoc',siteLoc,'cax',cax,'title',ctitle);
    plotCond(mPlot,grid,OPTIONS);
    CondPlotFile = [modelBaseName '.jpg'];
    eval(['print -djpeg90 ' fullfile(modelPath,CondPlotFile)]);
end

siteLoc = [siteLoc zeros(nSites,1)];
Z = zeros(nSites,1);
Zerr = ones(nSites,1);
T1 = conf.T1;
T2 = conf.T2;
nPer = conf.nPer;
if nPer > 1
    dt = (log10(T2/T1))/(nPer-1);
    periods = 10.^[log10(T1):dt:log10(T2)];
else
    periods = T1;
end

if MODE == 'JT'
   %  joint TE/TM (at all sites)
   %   first all TE
   allData = cell(2*nPer,1);
   for j=1:nPer
      allData{j} = struct('T',periods(j),...
	'Mode','TE',...
	'siteLoc',siteLoc,...
	'siteChar',char(siteChar),...
	'Z',Z,'Zerr',Zerr,...
    'Cmplx',1);
   end
   %   then all TM
   for j=1:nPer
      allData{j+nPer} = struct('T',periods(j),...
	'Mode','TM',...
	'siteLoc',siteLoc,...
	'siteChar',char(siteChar),...
	'Z',Z,'Zerr',Zerr,...
    'Cmplx',1);
   end
else
   % only one of TE and TM ...
   allData = cell(nPer,1);
   for j=1:nPer
      allData{j} = struct('T',periods(j),...
	'Mode',MODE,...
	'siteLoc',siteLoc,...
	'siteChar',char(siteChar),...
	'Z',Z,'Zerr',Zerr,...
    'Cmplx',1);
   end
end

% define the output data structure
data = allData;
   
if COMPUTE_DATA
   %  compute predicted data
   %[d] = fwdPred(logEarthCond,allData);
   [status] = writeZ_2D(dataTemplate,allData);
   [status] = Test2D('FORWARD',modelFile,dataTemplate,dataFile)
   [data] = readZ_2D(dataFile);
   if ADD_NOISE
      display(['adding noise ... FracError = ' num2str(FracError)])
      data = addNoise_2D(data,FracError);
   end
   % update output data structure
   [status] = writeZ_2D(dataFile,data);
elseif WRITE_FILES
   [status] = writeZ_2D(dataFile,data);
end

if COMPUTE_DATA && PLOT_DATA
    OPTIONS.rho_cax = [0 3.5];
    OPTIONS.phi_cax = [0 90];
    OPTIONS.title = ...
        [modelBaseName ' : Relative Error = ' num2str(100*FracError) '%'];
    clear h
    h = plotPseud(data,OPTIONS);
    if length(h) == 2
        figure(h(1))
        PseudPlotFile = ['TE_PSEUD_' dataBaseName '.jpg'];
        eval(['print -djpeg90 ' fullfile(dataPath,PseudPlotFile)]);
        figure(h(2))
        PseudPlotFile = ['TM_PSEUD_' dataBaseName '.jpg'];
        eval(['print -djpeg90 ' fullfile(dataPath,PseudPlotFile)]);
    else
        figure(h)
        PseudPlotFile = [MODE '_PSEUD_' dataBaseName '.jpg'];
        eval(['print -djpeg90 ' fullfile(dataPath,PseudPlotFile)]);
    end
end
