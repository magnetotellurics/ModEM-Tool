function [cond] = readCond_3D(cfile,ForPC)
%
% Usage:  [cond] = readCond_3D(cfile)
%         [cond] = readCond_3D(cfile,ForPC)
%
%  Reads in Cond3D object, returns as structure
%
%
% Set optional input ForPC = 1 for binary
% files compatible with PC = 0 for sun workstation
% Default (i.e., ForPC absent) assumes native

if nargin == 1
   fid = fopen(cfile,'r','n');
else
   if ForPC
      fid = fopen(cfile,'r','l');
   else
      fid = fopen(cfile,'r','b');
   end
end

ll = fread(fid,1,'long');
nm = fread(fid,4,'long');
nx = nm(1)
ny = nm(2);
nzearth = nm(3);
nzair = nm(4);
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
dx = fread(fid,[nx,1],'float64');
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
dy = fread(fid,[ny,1],'float64');
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
dz = fread(fid,[nzearth,1],'float64');
ll = fread(fid,1,'long');
grid = struct('Nx',nx,'Ny',ny,'NzEarth',nzearth,'NzAir',nzair,...
    'dx',dx,'dy',dy','dz',dz');

ll = fread(fid,1,'long');
s = fread(fid,80,'char');
s = deblank(char(s'));
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
AirCond = fread(fid,1,'float64');
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
v = fread(fid,[nx*ny*nzearth,1],'float64');
v = reshape(v,[nx,ny,nzearth]);
cond = struct('v',v,'paramType',s,'AirCond',AirCond,'grid',grid);
fclose(fid);
