%  prompt user for model file name
 modeldir = '.';
%  note: if we use suffixes to name grid files (e.g. *.grd)
%   we can change following line to
%   filt = [modeldir '/*.grd'];
filt = [modeldir '/*.ws; *.cpr; *.rho; *.prm'];
[filename, pathname] = uigetfile(filt, 'Choose Resistivity Model File');
[filename_pr, pathname_pr] = uigetfile(filt, 'Now Choose the Prior Model File');

%  read in chosen files
cfile = [pathname filename];
cfile_pr = [pathname_pr filename_pr];

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

Prior = readCond_3D(cfile_pr,2);

% convert to log10
if strcmp(Prior.paramType,'LOGE')
    Prior.paramType = 'LOG10';
    Prior.v = Prior.v / log(10);
    Prior.AirCond = Prior.AirCond / log(10);
elseif strcmp(Prior.paramType,'LINEAR')
    Prior.paramType = 'LOG10';
    Prior.v = log10(Prior.v);
    Prior.AirCond = log10(Prior.AirCond);
end

% take the difference!
Cond.v = Cond.v - Prior.v;

% now plot
options.slice = 'Z';
options.Np = 1;
[Nx,Ny,Nz] = size(Cond.v);
options.iXlim(1) = 1;
options.iXlim(2) = Nx+1;
options.iYlim(1) = 1;
options.iYlim(2) = Ny+1;
options.iZlim(1) = 1;
options.iZlim(2) = Nz+1;
options.cblabel = 'log_{10} (\sigma/\sigma_0)';
options.cmap = 'redblue';
CondPlotSet(Cond.v,Cond.grid,options);