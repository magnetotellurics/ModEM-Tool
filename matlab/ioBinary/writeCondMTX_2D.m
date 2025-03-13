function [status] = writeCondMTX_2D(cfile,CONDS,header,ForPC)
%
% Usage: [status] = writeCondMTX_2D(cfile,CONDS,header)
%        [status] = writeCondMTX_2D(cfile,CONDS,header,ForPC)
%
%  Write sequence of Cond2D objects stored in cell 
%   array of structures 
%
% Set optional input ForPC = 1 for binary
% files compatible with PC = 0 for sun workstation
% Default (i.e., ForPC absent) assumes native


if nargin == 3
   fid = fopen(cfile,'w','n');
else
   if ForPC
      fid = fopen(cfile,'w','l');
   else
      fid = fopen(cfile,'w','b');
   end
end

ll = 80;
status = 0;
count = fwrite(fid,ll,'long');
status = status+4*count;
s = blanks(80);
s(1:length(header)) = header;
count = fwrite(fid,s,'char');
status = status+count;
count = fwrite(fid,ll,'long');
status = status+4*count;
nSigma = length(CONDS);
ll = 4;
count = fwrite(fid,ll,'long');
status = status+count;
count = fwrite(fid,nSigma,'long')
status = status+count;
count = fwrite(fid,ll,'long');
status = status+count;

for k = 1:nSigma
   ll = 80;
   s = blanks(80);
   s(1:length(CONDS{k}.paramType)) = CONDS{k}.paramType;
   count = fwrite(fid,ll,'long');
   status = status+4*count;
   count = fwrite(fid,s,'char');
   status = status+count;
   count = fwrite(fid,ll,'long');
   status = status+4*count;

   Nyz = size(CONDS{k}.v);
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
   count = fwrite(fid,CONDS{k}.v,'float64');
   status = status+8*count;
   count = fwrite(fid,ll,'long');
   status= status+4*count;
   ll = 8;
   count = fwrite(fid,ll,'long');
   status = status+4*count;
   count = fwrite(fid,CONDS{k}.AirCond,'float64');
   status = status+8*count;
   count = fwrite(fid,ll,'long');
   status = status+4*count;
end
fclose(fid);
