%  read grid file at start
ReadGrid = 0;
%  prompt user for model file name
fdir = '.';
%  note: if we use suffixes to name grid files (e.g. *.grd)
%   we can change following line to
%   filt = [fdir '/*.grd'];
filt = [fdir '/*.soln'];
if(ReadGrid)
   %   First get grid file  (need grid until new headers work)
   [filename, pathname] = uigetfile(filt, 'Grid File');
   cfile = [pathname filename];
   [ModelGrid,rho] = rdModelRM(cfile);
   dzAir = 10*(3.^[1:ModelGrid.NzAir]);
   ModelGrid.dz = [dzAir(end:-1:1)' ; ModelGrid.dz];
   S = EdgeCond(ModelGrid,rho);
else
   S = struct('x',[],'y',[],'z',[],'CondRead',0);
end

%  then get solution file
[filename, pathname] = uigetfile(filt, 'Solution File');
cfile = [pathname filename];
fNum = 1; mNum=1;
[E,T,Grid,Modes] = rdExyz(cfile, fNum, mNum);

options.T = T;
options.iPer = fNum;
options.nPer = length(T);
options.slice = 'Z';
options.Np = Grid.NzAir+1;
options.mode = mNum;
options.Comp = 'Y';
options.SolnFile = cfile;
options.Modes = Modes;

[Nx,Ny,Nz] = size(E.y);
options.iXlim(1) = 1;
options.iXlim(2) = Nx;
options.iYlim(1) = 1;
options.iYlim(2) = Ny;
options.iZlim(1) = Grid.NzAir+1;
options.iZlim(2) = Nz;
EplotSet(E,S,Grid,options);
