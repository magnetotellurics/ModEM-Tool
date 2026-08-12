function [status] = writeCond_2D(cfile,cond,grid,ForPC)
%
% Usage:  [status] = writeCond_2D(cfile,cond,grid)
%         [status] = writeCond_2D(cfile,cond,grid,ForPC)  
%
% Writes CondParam object, provided as structure
% status is total number of bytes written
%
% If optional input argument grid is present, use it;
% otherwise, try to extract it from cond
%
% Set optional input ForPC = 1 for binary
% files compatible with PC = 0 for sun workstation
% Default (i.e., ForPC absent) assumes native

if nargin < 3
    if isfield(cond,'grid')
        grid = cond.grid;
    else
        error('Please pass grid as the 3rd input argument to writeCond2D')
    end
end

if nargin < 4
   fid = fopen(cfile,'w','n');
else
   if ForPC
      fid = fopen(cfile,'w','l');
   else
      fid = fopen(cfile,'w','b');
   end
end

status = 0;
ll = 12;
count = fwrite(fid,ll,'long');
status = status+4*count;
Iword = [grid.Ny,grid.Nz,grid.Nza];
count = fwrite(fid,Iword,'long');
status = status+4*count;
count = fwrite(fid,ll,'long');
status = status+ 4*count;
ll = 8*grid.Ny;
count = fwrite(fid,ll,'long');
status = status+ 4*count;
count = fwrite(fid,grid.Dy,'float64');
status = status+ 8*count;
count = fwrite(fid,ll,'long');
status = status+ 4*count;
ll = 8*grid.Nz;
count = fwrite(fid,ll,'long');
status = status+ 4*count;
count = fwrite(fid,grid.Dz,'float64');
status = status+ 8*count;
count = fwrite(fid,ll,'long');
status = status+ 4*count;

ll = 80;
count = fwrite(fid,ll,'long');
status = status+4*count;
s = blanks(80);
s(1:length(cond.paramType)) = cond.paramType;
count = fwrite(fid,s,'char');
status = status+count;
count = fwrite(fid,ll,'long');
status = status+4*count;
Nyz = size(cond.v);
ll = 8;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,Nyz,'long');
status = status+4*count;
count = fwrite(fid,ll,'long');
status = status+4*count;
ll = Nyz(1)*Nyz(2)*8;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,cond.v,'float64');
status = status+8*count;
count = fwrite(fid,ll,'long');
status= status+4*count;
ll = 8;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,cond.AirCond,'float64');
status = status+8*count;
count = fwrite(fid,ll,'long');
status = status+4*count;

fclose(fid);
