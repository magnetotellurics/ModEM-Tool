function [status] = writeZ_3D(cfile,allData,ForPC)
%
%  Usage:  [status] = writeZ_3D(cfile,allData);
%          [status] = writeZ_3D(cfile,allData,ForPC);
%
%   write contents of cell array allData to file
%   cfile, for 3D MT data.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   Slightly modified from writeZ for TE+TM 2D MT
%
%   ForPC is an optional argument, true (1) to produce
%   files which can be read on PC (linux/windows)

if nargin == 2
   fid = fopen(cfile,'w','n');
else
   if ForPC
      fid = fopen(cfile,'w','l');
   else
      fid = fopen(cfile,'w','b');
   end
end

nTx = length(allData);
status = 0;

%  record # 1: number of transmitters 
ll = 4;
count = fwrite(fid,ll,'long');
status = status + 4*count; 
count = fwrite(fid,nTx,'long');
status = status + 4*count; 
count = fwrite(fid,ll,'long');
status = status + 4*count; 

for j = 1:nTx
%   for each transmitter (period/mode) ....

   [nSites,nZ] = size(allData{j}.Z);
   nComp = nZ*2;
   %  record # 1 in transmitter block: period, # of components, number of sites 
   ll =  16;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   count = fwrite(fid,allData{j}.T,'float64');
   status = status + 8*count; 
   count = fwrite(fid,nComp,'long');
   status = status + 4*count; 
   count = fwrite(fid,nSites,'long');
   status = status + 4*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

   %  record # 2 in transmitter block: site locations
   ll = 3*8*nSites;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   count = fwrite(fid,allData{j}.siteLoc','float64');
   status = status + 8*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

   %  record # 3 in transmitter block: data
   ll = nComp*8*nSites;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   Zri = reshape(allData{j}.Z,[1,nZ*nSites]);
   Zri = [real(Zri);imag(Zri)];
   count = fwrite(fid,Zri,'float64');
   status = status + 8*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

   %  record # 4 in transmitter block: error bars
   %  NOTE: write out same errors for real and imag parts!
   ll = nComp*8*nSites;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   %  NOTE: write out same errors for real and imag parts!
   Zri = reshape(allData{j}.Zerr,[1,nZ*nSites]);
   Zri = [real(Zri);real(Zri)];
   count = fwrite(fid,Zri,'float64');
   status = status + 8*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

end
fclose(fid);
