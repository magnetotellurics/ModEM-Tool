function [CONDS,header] = readCondMTX_2D(cfile,ForPC)
%
% Usage:  [CONDS,header] = readCondMTX_2D(cfile)
%         [CONDS,header] = readCondMTX_2D(cfile,ForPC)
%
%  Reads in sequence of Cond2D objects, return as cell
%   array of structures ... one Cond2D object for each
%   real observation
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

fread(fid,1,'long');
header = fread(fid,80,'char');
fread(fid,1,'long');
fread(fid,1,'long');
nSigma = fread(fid,1,'long');
fread(fid,1,'long');

CONDS = cell(nSigma,1);
for k = 1:nSigma
   fread(fid,1,'long');
   s = fread(fid,80,'char');
   s = deblank(char(s'));
   fread(fid,1,'long');
   fread(fid,1,'long');
   nm = fread(fid,2,'long');
   ny = nm(1);
   nz = nm(2);
   fread(fid,1,'long');
   fread(fid,1,'long');
   v = fread(fid,[ny,nz],'float64');
   fread(fid,1,'long');
   fread(fid,1,'long');
   AirCond = fread(fid,1,'float64');
   fread(fid,1,'long');
   CONDS{k} = struct('v',v,'paramType',s,'AirCond',AirCond);
end
fclose(fid);
