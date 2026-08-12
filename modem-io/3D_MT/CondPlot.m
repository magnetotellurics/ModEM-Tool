%  prompt user for model file name
 modeldir = '.';
%  note: if we use suffixes to name grid files (e.g. *.grd)
%   we can change following line to
%   filt = [modeldir '/*.grd'];
filt = [modeldir '/*.rm; *.cpr'];
[filename, pathname] = uigetfile(filt, 'Resistivity Model File');

%  read in chosen file
cfile = [pathname filename];

% [grid,rho] = rdModelRM(cfile);
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

options.slice = 'Z';
options.Np = 1;
[Nx,Ny,Nz] = size(Cond.v);
options.iXlim(1) = 1;
options.iXlim(2) = Nx+1;
options.iYlim(1) = 1;
options.iYlim(2) = Ny+1;
options.iZlim(1) = 1;
options.iZlim(2) = Nz+1;
options.cblabel = 'log_{10} \sigma';
CondPlotSet(Cond.v,Cond.grid,options);
