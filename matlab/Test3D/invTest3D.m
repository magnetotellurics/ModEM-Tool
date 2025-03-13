%  Script invTest is used to test NLCG inversion in fortran 90
%  NOTE: Test3D must be in your path or current directory
%% settings and file directory
testDir = 'SENSTEST/';
config = [testDir 'config.mat'];
%% create conf structure
if exist(config,'file')
    load(config)
else
    Config3D
end
LOGCOND = 1;

%% File names
m0_File = [testDir 'm0.ws']; % model (true)
pr_File = [testDir 'pr.ws']; % model (prior)
mi_File = [testDir 'mi.ws']; % model (inverse)
de_File = [testDir 'de.dat']; % data (with errors) 
di_File = [testDir 'di.dat']; % data (for the inverse model) 
if exist(testDir,'dir')==0
    mkdir(testDir);
end

[m0_Path, m0_BaseName] = fileparts(m0_File); 
[mi_Path, mi_BaseName] = fileparts(mi_File); 
[de_Path, de_BaseName] = fileparts(de_File); 
[di_Path, di_BaseName] = fileparts(di_File); 

%% Create the background model and data files and structures
conf.COMPUTE_DATA = 1;
conf.ADD_NOISE = 1;
conf.PLOT_MODEL = 0;
conf.PLOT_DATA = 0;
[m0,de,pr] = mkCond3D(conf,m0_File,de_File)

%% Write input files
writeCond_3D(pr_File,pr,2);
writeZ_3D(de_File,de);

%% Run forward solver
Test3D('INVERSE_NLCG',pr_File,de_File)
