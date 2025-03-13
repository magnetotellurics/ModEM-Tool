% Usage: mkCond2D(conf) to create the 2-D model and data
%        sensTest2D to run the full sensitivity test

% choose to use or not use log cond for computations
conf.LOGCOND = 1;

% choose to generate the model or run the forward solver
conf.COMPUTE_DATA = 0;
conf.ADD_NOISE = 0;

% choose the nature of data and data errors
conf.mode = 'JT'; % 'TE', 'TM', 'JT'
conf.data_error = 0.03;

% set the values of background and anomalous conductivity
bg = .01;
a1 = .1;
a2 = .03;
a3 = .001;
a4 = .0003;
%bg = 0.01;
%a1 = 0.5;

% choose a subset of anomalies
conf.nBlocks = 9;
conf.bg_cond = bg;
%conf.anomaly = [bg a1 a3 a2 bg bg a3 bg bg]; % 3 blocks
%conf.anomaly = [bg bg bg a1 a1 a3 a3 bg bg]; % 2 blocks
%conf.anomaly = [a3 a1 a3 a1 a3 a1 a3 a1 a3]; % rubik's cube
conf.anomaly = [a1 a3 a1 a3 a1 a3 a1 a3 a1]; % rubik's cube inverted

% (y0:y1) max range from 1 to 100
conf.y0      = [16 16 16 46 46 46  76  76  76]; % rubik's cube
conf.y1      = [40 40 40 70 70 70 100 100 100]; % rubik's cube
%conf.y0      = [1  11 21 31 41 51 61 71 81]; % 2 blocks
%conf.y1      = [10 20 30 40 50 60 70 80 90]; % 2 blocks
%conf.y0      = [5  16 26 36 46 56 66 76 86]; % 3 blocks
%conf.y1      = [15 35 35 55 55 65 75 95 95]; % 3 blocks
%conf.y0      = [1  16 26 36 46 56 66 76 86]; % random
%conf.y1      = [50 35 35 55 55 65 75 95 95]; % random

% (z0:z1) max range from 1 to 31
conf.z0      = [1  11 18  1 11 18  1 11 18]; % rubik's cube
conf.z1      = [8  16 23  8 16 23  8 16 23]; % rubik's cube
%conf.z0      = [1  10  1 15  1 10 10  1 15]; % 3 blocks
%conf.z1      = [15 18  5 18 10 18 15 10 18]; % 3 blocks
%conf.z0      = [1  10 21 15 21 10 10 21 15]; % random
%conf.z1      = [9  18 25 18 30 18 15 30 18]; % random

% for data, choose the periods and number of sites etc
conf.NODE_SITE = 0;
conf.nSites = 30;
conf.FirstSite = 10;
conf.LastSite = 106;
conf.nPer = 12;
conf.T1 = 0.3; %0.3;
conf.T2 = 3000; %3000;

% horizontal grid parameters
conf.Ny = 106;
conf.NyPad = 5;
conf.dyBase = 1500;
conf.yPadFac = 1.5;

% vertical grid parameters
conf.Nz = 40;
conf.dzBase = 500; 
conf.zFac = 1.1;

% choose whether to plot things
conf.PLOT_MODEL = 1;
conf.PLOT_DATA = 1;

% choose whether to save this configuration
conf.SAVE_CONF = 1;

% use ascii or binary files
conf.BINARY = 0;

% architecture dependency for binary files
conf.ForPC = 1;