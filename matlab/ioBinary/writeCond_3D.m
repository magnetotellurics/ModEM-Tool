function [status] = writeCond_3D(cfile,cond,ForPC)
%
% Usage:  [status] = writeCond_3D(cfile,cond)
%         [status] = writeCond_3D(cfile,cond,FOrPC)
%
%  writes ModelParam object, provided as structure
%  status is total number of bytes written
%
% Set optional input ForPC = 1 for binary
% files compatible with PC = 0 for sun workstation
% Default (i.e., ForPC absent) assumes native

if nargin == 2
   fid = fopen(cfile,'w','n');
else
   if ForPC
      fid = fopen(cfile,'w','l');
   else
      fid = fopen(cfile,'w','b');
   end
end

grid = cond.grid;
status = 0;

%   write array size record
Nxyz = [grid.Nx grid.Ny grid.NzEarth grid.NzAir];
ll = 16;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,Nxyz,'long');
status = status+4*count;
count = fwrite(fid,ll,'long');
status = status+4*count;

%   write dx record
ll = length(grid.dx)*8;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,grid.dx,'float64');
count = fwrite(fid,ll,'long');
status = status+4*count;

%   write dy record
ll = length(grid.dy)*8;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,grid.dy,'float64');
count = fwrite(fid,ll,'long');
status = status+4*count;

%   write dz record
ll = length(grid.dz)*8;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,grid.dz,'float64');
count = fwrite(fid,ll,'long');
status = status+4*count;

%   write paramType record
ll = 80;
count = fwrite(fid,ll,'long');
status = status+4*count;
s = blanks(80);
s(1:length(cond.paramType)) = cond.paramType;
count = fwrite(fid,s,'char');
status = status+count;
count = fwrite(fid,ll,'long');
status = status+4*count;

%  write AirCond record
ll = 8;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,cond.AirCond,'float64');
status = status+8*count;
count = fwrite(fid,ll,'long');
status = status+4*count;

%  write conductivity array record
ll = Nxyz(1)*Nxyz(2)*Nxyz(3)*8;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,cond.v,'float64');
status = status+8*count;
count = fwrite(fid,ll,'long');
status= status+4*count;

fclose(fid);
