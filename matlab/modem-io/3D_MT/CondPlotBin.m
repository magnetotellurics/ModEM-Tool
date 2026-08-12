%  prompt user for model file name
 modeldir = '.';
%  note: if we use suffixes to name grid files (e.g. *.grd)
%   we can change following line to
%   filt = [modeldir '/*.grd'];
filt = [modeldir '/*.cpr'];
[filename, pathname] =  ....
uigetfile(filt, 'Binary Conductivity Model Parameter File');

%  read in chosen file
cfile = [pathname filename];
condParam = readCond_3D(cfile)
grid = struct('dx',condParam.dx,'dy',condParam.dy,...
	'dz',condParam.dz,'nZair',0,'origin',[0,0,0],...
        'theta',0);

%[grid,rho] = rdModelRM(cfile);
Cond = exp(condParam.v);

options.slice = 'Z';
options.Np = 1;
[Nx,Ny,Nz] = size(Cond);
options.iXlim(1) = 1;
options.iXlim(2) = Nx+1;
options.iYlim(1) = 1;
options.iYlim(2) = Ny+1;
options.iZlim(1) = 1;
options.iZlim(2) = Nz+1;
CondPlotSet(Cond,grid,options);
