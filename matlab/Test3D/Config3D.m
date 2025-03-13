% Usage: mkCond3D(conf) to create the 3-D model and data
%        sensTest3D to run the full sensitivity test

% choose to use or not use log cond for computations
conf.LOGCOND = 1;

% choose to generate the model or run the forward solver
conf.COMPUTE_DATA = 0;
conf.ADD_NOISE = 0;

% choose the nature of data and data errors
conf.OFF_DIAG = 0;
conf.data_error = 0.03;

% set the values of background and anomalous conductivity
bg = .01;
a1 = .1;
a2 = .03;
a3 = .001;
a4 = .0003;

% choose a subset of anomalies
conf.bg_cond = bg;
conf.anomaly = [a2 a3 bg bg bg bg]; % 6 values

% (x0:x1) max range from 1 to 21
conf.x0      = [ 6  6  6  6  6  6];
conf.x1      = [16 16 16 16 16 16];

% (y0:y1) max range from 1 to 28
conf.y0      = [ 7 15  7 15  7 15];
conf.y1      = [14 22 14 22 14 22];

% (z0:z1) max range from 1 to 11
conf.z0      = [ 2  2  6  6 10 10];
conf.z1      = [ 6  6  8  8 11 11];

% for data, choose the periods and number of sites etc
conf.NODE_SITE = 0;
conf.nSites = 'all'; % 3, 6, 30 or 'all'
conf.nPer = 4;
conf.T1 = 1;
conf.T2 = 1000;

% choose whether to plot things
conf.PLOT_MODEL = 1;
conf.PLOT_DATA = 1;

% choose whether to save this configuration
conf.SAVE_CONF = 1;

% use ascii or binary files
conf.BINARY = 0;

% architecture dependency for binary files
conf.ForPC = 1;