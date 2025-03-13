function [cond,grid] = readCond_2D(cfile,ForPC)
%
% Usage:  [cond,grid] = readCond_2D(cfile)
%         [cond,grid] = readCond_2D(cfile,ForPC)
%
%  Reads in Cond2D object, returns as structure
%  Optional output argument grid also part of cond
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
Ny = fread(fid,1,'long');
Nz = fread(fid,1,'long');
Nza = fread(fid,1,'long');
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
Dy = fread(fid,Ny,'float64');
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
Dz = fread(fid,Nz,'float64');
ll = fread(fid,1,'long');
grid = struct('Ny',Ny,'Nz',Nz,'Nza',Nza,'Dy',Dy','Dz',Dz');

ll = fread(fid,1,'long');
s = fread(fid,80,'char');
s = deblank(char(s'));
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
nm = fread(fid,2,'long');
ny = nm(1);
nz = nm(2);
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
v = fread(fid,[ny,nz],'float64');
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
AirCond = fread(fid,1,'float64')
cond = struct('v',v,'paramType',s,'AirCond',AirCond,'grid',grid);
fclose(fid);
