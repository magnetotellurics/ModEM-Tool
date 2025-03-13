Config3D

%% File names
conf.COMPUTE_DATA = 0;
[m0,dT] = mkCond3D(conf)
m0_File = 'm0.cpr';
dT_File = 'Template_d0_3D.imp';
d0_File = 'd0_3D.imp';

%% Write input files
writeCond_3D(m0_File,m0);
writeZ_3D(dT_File,dT);

%% Run forward solver
Test3D('FORWARD',m0_File,dT_File,d0_File)

%% Read output file
d0 = readZ_3D(d0_File)
