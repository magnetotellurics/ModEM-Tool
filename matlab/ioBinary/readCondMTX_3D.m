function [CONDS,header] = readCondMTX_3D(cfile,ForPC)

% Usage:  [CONDS,header] = readCondMTX_3D(cfile)
%         [CONDS,header] = readCondMTX_3D(cfile,ForPC)
%
%  Reads in sequence of 3D ModelParam objects, returns as cell
%   array of structures ... one ModelParam object for each
%   real observation, when used to read a sensitivity matrix
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
fread(fid,80,'char');
fread(fid,2,'long');
nSigma = fread(fid,1,'long');
CONDS  = cell(nSigma,1);
fread(fid,1,'long');

m = MT3DmodelParam();
CONDS = cell(nSigma,1);
for k = 1:nSigma
   fread(fid,1,'long');
   s = fread(fid,80,'char');
   paramType = deblank(char(s'));
   fread(fid,2,'long');
   nm = fread(fid,3,'long');
   nx = nm(1);
   ny = nm(2);
   nz = nm(3);
   fread(fid,2,'long');
   dx = fread(fid,[nx,1],'float64');
   fread(fid,2,'long');
   dy = fread(fid,[ny,1],'float64');
   fread(fid,2,'long');
   dz = fread(fid,[nz,1],'float64');
   fread(fid,2,'long');
   AirCond = fread(fid,1,'float64');
   fread(fid,2,'long');
   v = fread(fid,[nx*ny*nz,1],'float64');
   v = reshape(v,[nx,ny,nz]);
   if k==1
       Grid = struct('dx',dx,'dy',dy,'dz',dz,'nx',nx,'ny',ny,...
           'nzEarth',nz,'NzAir',0,'rotation',0,'origin',[0,0,0]);
   end
   CONDS{k} = m.setModelParam(v,paramType,AirCond,Grid);
   fread(fid,1,'long');
end
fclose(fid);
