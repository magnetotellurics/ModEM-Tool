%  Script invTest is used to test NLCG inversion in fortran 90
%  NOTE: Test2D must be in your path or current directory
%% settings and file directory
testDir = 'SENSTEST/';
config = [testDir 'config.mat'];
%% create conf structure
if exist(config,'file')
    load(config)
else
    Config2D
end
MODE = conf.mode;
LOGCOND = 1;
PLOT_MODEL = 1;

%% File names
m0_File = [testDir 'm0.rho']; % model (true)
pr_File = [testDir 'pr.rho']; % model (prior)
mi_File = [testDir 'mi.rho']; % model (inverse)
de_File = [testDir 'de_' MODE '.dat']; % data (with errors) 
di_File = [testDir 'di_' MODE '.dat']; % data (for the inverse model) 
if exist(testDir,'dir')==0
    mkdir(testDir);
end

[m0_Path, m0_BaseName] = fileparts(m0_File); 
[mi_Path, mi_BaseName] = fileparts(mi_File); 
[de_Path, de_BaseName] = fileparts(de_File); 
[di_Path, di_BaseName] = fileparts(di_File); 

%% create the background model and data files and structures
conf.COMPUTE_DATA = 1;
conf.ADD_NOISE = 1;
conf.PLOT_MODEL = 0;
conf.PLOT_DATA = 0;
[m0,de,pr] = mkCond2D(conf,m0_File,de_File)

%% save the (manually, if required) edited configuration
if conf.SAVE_CONF
    save(config,'conf')
end

%% create a random perturbation in the background model
writeCond_2D(pr_File,pr);

%% Run forward solver
Test2D('INVERSE_NLCG',pr_File,de_File)

%% Save inverse model and responses
system(['cp  Modular_NLCG_042.rho ' mi_File])
system(['cp  Modular_NLCG_042.dat ' di_File])

%% Read output file
mi = readCond_2D(mi_File);
di = readZ_2D(di_File);

%% Plot inverse model
if PLOT_MODEL
    nZplot = 25;
    nYskip = 7;
    cax = [-3.5,-.5];  %   +log10(2);  Model 2

    ctitle = 'Inverse model';
    OPTIONS = struct('nYskip',nYskip,'nZplot',nZplot,...
        'siteLoc',de{1}.siteLoc(:,1),'cax',cax,'title',ctitle);
    plotCond(mi,mi.grid,OPTIONS);
    CondPlotFile = [mi_BaseName '.jpg'];
    eval(['print -djpeg90 ' fullfile(mi_Path,CondPlotFile)]);

    ctitle = 'True model';
    OPTIONS = struct('nYskip',nYskip,'nZplot',nZplot,...
        'siteLoc',de{1}.siteLoc(:,1),'cax',cax,'title',ctitle);
    plotCond(m0,m0.grid,OPTIONS);
    CondPlotFile = [m0_BaseName '.jpg'];
    eval(['print -djpeg90 ' fullfile(m0_Path,CondPlotFile)]);

    ctitle = 'Inverse model';
end

%% Plot data comparison
plotBands = [1 4 7 10];
for k = 1:length(plotBands)
   T = de{plotBands(k)}.T;
   figure('Position',[100,100,400,700])
   axes('Position',[.1,.55,.8,.35])
   plot(real(de{plotBands(k)}.Z),'ks')
   hold on
   plot(real(di{plotBands(k)}.Z),'r-')
   title(['Period = ', num2str(T)])
   axes('Position',[.1,.1,.8,.35])
   plot(imag(de{plotBands(k)}.Z),'ks')
   hold on
   plot(imag(di{plotBands(k)}.Z),'r-')
end
