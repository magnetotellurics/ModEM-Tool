function [model,data,prior] = mkCond3D(conf,modelFile,dataFile,EsolnFile)

% Usage: [model,data,prior] = mkCond3D(conf,modelFile,dataFile,EsolnFile)
%
% Inputs are file names and configuration structure created by
% Config3D. Outputs are matlab structures; prior is the optional
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
OFF_DIAG = conf.OFF_DIAG;
FracError = conf.data_error;
%NODE_SITE = conf.NODE_SITE;
%PLOT_MODEL = conf.PLOT_MODEL;
%PLOT_DATA = conf.PLOT_DATA;

if nargin < 3 && COMPUTE_DATA
    'Unable to compute data: no data file specified'
    COMPUTE_DATA = 0;
end

%  strip all extensions off file names
if WRITE_FILES
    [modelPath, modelBaseName, ext] = fileparts(modelFile);
    if COMPUTE_DATA
        [dataPath, dataBaseName, ext] = fileparts(dataFile);
        dataTemplate = fullfile(dataPath,['Template_' dataBaseName ext]);
    end
else
    modelPath = pwd;
    modelBaseName = 'model';
end

%   here the basic grid geometry is defined
Dx = [20000 10000 5000 3000 2000 2000 2500 4500 4500 4500 ...
        4000 4500 4500 4500 2500 2000 2000 3000 5000 10000 20000];
Dy = [ 20000 10000 5000 2500 1250 1250 1250 1250 2500 5000 5000 ...
        2500 1250 1250 1250 1250 2500 5000 5000 2500 1250 1250 ...
        1250 1250 2500 5000 10000 20000];
Dz = [ 500 1000 1500 3000 4000 5000 5000 10000 10000 20000 40000];

nx = length(Dx);
ny = length(Dy);
nzEarth = length(Dz);
nzAir = 5;
nz = nzEarth+nzAir;

%   here the conductivity model is defined: background with
%    rectangular anomalies of the form Ix# x IY# x IZ#
% from larger blocks to smaller, up to down
% example starts
% Ix1 = [6:16]; Iy1 = [7:14];  Iz1 = [1:5];
% Ix2 = [6:16]; Iy2 = [15:22];  Iz2 = [1:5];
% example ends

Ix1 = (conf.x0(1):conf.x1(1)); Iy1 = (conf.y0(1):conf.y1(1));  Iz1 = (conf.z0(1):conf.z1(1));
Ix2 = (conf.x0(2):conf.x1(2)); Iy2 = (conf.y0(2):conf.y1(2));  Iz2 = (conf.z0(2):conf.z1(2));
Ix3 = (conf.x0(3):conf.x1(3)); Iy3 = (conf.y0(3):conf.y1(3));  Iz3 = (conf.z0(3):conf.z1(3));
Ix4 = (conf.x0(4):conf.x1(4)); Iy4 = (conf.y0(4):conf.y1(4));  Iz4 = (conf.z0(4):conf.z1(4));
Ix5 = (conf.x0(5):conf.x1(5)); Iy5 = (conf.y0(5):conf.y1(5));  Iz5 = (conf.z0(5):conf.z1(5));
Ix6 = (conf.x0(6):conf.x1(6)); Iy6 = (conf.y0(6):conf.y1(6));  Iz6 = (conf.z0(6):conf.z1(6));

backGroundCond = conf.bg_cond;
airCond = 1e-10;
Cond = backGroundCond*ones(nx,ny,nzEarth);
Cond(Ix1,Iy1,Iz1) = conf.anomaly(1);
Cond(Ix2,Iy2,Iz2) = conf.anomaly(2);
Cond(Ix3,Iy3,Iz3) = conf.anomaly(3);
Cond(Ix4,Iy4,Iz4) = conf.anomaly(4);
Cond(Ix5,Iy5,Iz5) = conf.anomaly(5);
Cond(Ix6,Iy6,Iz6) = conf.anomaly(6);

%  create grid structure
grid = struct('Nx',nx,'Ny',ny,'NzEarth',nzEarth,'NzAir',nzAir,...
		'dx',Dx,'dy',Dy,'dz',Dz);

%  write out model/grid file in RM format ... 
%origin = [0,0,0];
%rotation = 0;
%rho = 1./Cond;
%writeRM(modelFile,Dx,Dy,Dz,nzAir,rho,origin,rotation);

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
    prior.v(:,:,:) = log(backGroundCond);
else
    model = earthCond;
    prior = model;
    prior.v(:,:,:) = backGroundCond;
end
if WRITE_FILES
    [status] = writeCond_3D(modelFile,model,2);
end

% if PLOT_MODEL
%     mPlot = logEarthCond;
%     nZplot = 18;
%     nYskip = 7;
%     cax = [-3.5,-.5];  %   +log10(2);  Model 2
%     ctitle = ['Relative Error = ' num2str(100*FracError) '%'];
%     OPTIONS = struct('nYskip',nYskip,'nZplot',nZplot,...
%         'cax',cax,'title',ctitle);
%     plotCond(mPlot,grid,OPTIONS);
%     CondPlotFile = [modelBaseName '.jpg'];
%     eval(['print -djpeg90 ' fullfile(modelPath,CondPlotFile)]);
% end

%  set up data vector template
x = cumsum(Dx);
y = cumsum(Dy);
z = cumsum(Dz);
x = [0 x];
y = [0 y];
z = [0 z];
xctr = x(1:end-1) + diff(x)/2;
yctr = y(1:end-1) + diff(y)/2;

%   DATA LOCATIONS:
if conf.nSites == 6
    %  For the initial test 6 sites are centered over anomaly
    Cells = [11 6 ; 11 9 ; 11 12; 11 17; 11 20 ;11 23];
    siteLoc = [(x(Cells(:,1))+x(Cells(:,1)+1))/2; ...
        (y(Cells(:,2))+y(Cells(:,2)+1))/2; ...
        zeros(size(x(Cells(:,1))))];
    siteChar = num2str([1:size(Cells,1)]');
    % first dimension is always the number of sites
    siteLoc = siteLoc';
elseif conf.nSites == 3
    %  .... or mini data set with only 3 sites ...
    Cells = [11 6 ; 11 12; 11 20 ];
    siteLoc = [(x(Cells(:,1))+x(Cells(:,1)+1))/2; ...
        (y(Cells(:,2))+y(Cells(:,2)+1))/2; ...
        zeros(size(x(Cells(:,1))))];
    siteChar = num2str([1:size(Cells,1)]');
    % first dimension is always the number of sites
    siteLoc = siteLoc';
elseif conf.nSites == 30
    %  ... or a large synthetic test with 30 sites
    Cells = [8 5; 8 8 ; 8 10 ; 8 11; 8 13; 8 16; 8 18; 8 19; 8 21; 8 24; ...
        11 5; 11 8 ; 11 10 ; 11 11; 11 13; 11 16; 11 18; 11 19 ;11 21; 11 24; ...
        14 5; 14 8 ; 14 10 ; 14 11; 14 13; 14 16; 14 18; 14 19 ;14 21; 14 24; ...
        ];
    siteLoc = [(x(Cells(:,1))+x(Cells(:,1)+1))/2; ...
        (y(Cells(:,2))+y(Cells(:,2)+1))/2; ...
        zeros(size(x(Cells(:,1))))];
    siteChar = num2str([1:size(Cells,1)]');
    siteChar(siteChar==' ') = '0';
    % first dimension is always the number of sites
    siteLoc = siteLoc';
elseif conf.nSites == 50
    %  ... or a large synthetic test with 30 sites
    Cells = [4 4; 4 7 ; 4 10 ; 4 11; 4 13; 4 15; 4 18; 4 19; 4 21; 4 24; ...
        8 4; 8 7 ; 8 10 ; 8 11; 8 13; 8 15; 8 18; 8 19; 8 21; 8 24; ...
        11 4; 11 7 ; 11 10 ; 11 11; 11 13; 11 15; 11 18; 11 19; 11 21; 11 24; ...
        14 4; 14 7 ; 14 10 ; 14 11; 14 13; 14 15; 14 18; 14 19; 14 21; 14 24; ...
        18 4; 18 7 ; 18 10 ; 18 11; 18 13; 18 15; 18 18; 18 19; 18 21; 18 24; ...
        ];
    siteLoc = [(x(Cells(:,1))+x(Cells(:,1)+1))/2; ...
        (y(Cells(:,2))+y(Cells(:,2)+1))/2; ...
        zeros(size(x(Cells(:,1))))];
    siteChar = num2str([1:size(Cells,1)]');
    siteChar(siteChar==' ') = '0';
    % first dimension is always the number of sites
    siteLoc = siteLoc';
else
    % create a data point at all cell centers except padding
    xpad = 5;
    ypad = 5;
    xdiff = 1; % create sites every xdiff cells (latitude)
    ydiff = 1; % create sites every ydiff cells (longitude)
    siteLoc = [];
    siteChar = '';
    for i = xpad+1:xdiff:grid.Nx-xpad
        for j = ypad+1:ydiff:grid.Ny-ypad
            siteLoc = [siteLoc; xctr(i) yctr(j) 0.0];
            siteChar = [siteChar; sprintf('%03d',i) '-' sprintf('%03d',j)];
        end
    end
end
nSites = size(siteLoc,1);


%   PERIODS:
T1 = conf.T1; %10;
T2 = conf.T2; %100;
nPer = conf.nPer; %2;
if nPer > 1
    dt = (log10(T2/T1))/(nPer-1);
    periods = 10.^[log10(T1):dt:log10(T2)];
else
    periods = T1;
end

if OFF_DIAG
    %  set up for off-diagonal impedances only
    nImp = 2;
    compChar = ['Re(Zxy)';'Im(Zxy)';'Re(Zyx)';'Im(Zyx)'];
else
    nImp = 6;
    compChar = ['Re(Zxx)';'Im(Zxx)';'Re(Zxy)';'Im(Zxy)'; ...
        'Re(Zyx)';'Im(Zyx)';'Re(Zyy)';'Im(Zyy)'; ...
        'Re(Tx) ';'Im(Tx) ';'Re(Ty) ';'Im(Ty) '];
end
Z = zeros(nSites,nImp)+1i*zeros(nSites,nImp);
Zerr = ones(nSites,nImp);

allData = cell(nPer,1);
for j=1:nPer
    allData{j} = struct('T',periods(j),...
        'nComp',nImp*2,...
        'compChar',compChar,...
        'siteLoc',siteLoc,...
        'siteChar',siteChar,...
        'Z',Z,'Zerr',Zerr,...
        'Cmplx',1);
end
%   end data vector template setup

data = allData;

if COMPUTE_DATA
   %  compute predicted data
   %[d] = fwdPred(logEarthCond,allData);
   [status] = writeZ_3D(dataTemplate,allData);
   if nargin > 3
    [status] = Test3D('FORWARD',modelFile,dataTemplate,dataFile,EsolnFile)
   else
    [status] = Test3D('FORWARD',modelFile,dataTemplate,dataFile)
   end
   [data] = readZ_3D(dataFile);
   if ADD_NOISE
      ['adding noise ... FracError = ' num2str(FracError)]  
      data = addNoise_3D(data,FracError);
   end
end
if WRITE_FILES
   [status] = writeZ_3D(dataFile,data);
end    

% if COMPUTE_DATA && PLOT_DATA
%     OPTIONS.rho_cax = [0 3.5];
%     OPTIONS.phi_cax = [0 90];
%     clear h
%     h = plotPseud(data,OPTIONS);
%     if length(h) == 2
%         figure(h(1))
%         PseudPlotFile = ['JT_TE_PSEUD_' dataBaseName '.jpg'];
%         eval(['print -djpeg90 ' fullfile(dataPath,PseudPlotFile)]);
%         figure(h(2))
%         PseudPlotFile = ['JT_TM_PSEUD_' dataBaseName '.jpg'];
%         eval(['print -djpeg90 ' fullfile(dataPath,PseudPlotFile)]);
%     else
%         figure(h)
%         PseudPlotFile = [MODE '_PSEUD_' dataBaseName '.jpg'];
%         eval(['print -djpeg90 ' fullfile(dataPath,PseudPlotFile)]);
%     end
% end