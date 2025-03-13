function [status] = writeZ_2D(cfile,allData,ForPC)
%
%  Usage:  [status] = writeZ_2D(cfile,allData);
%          [status] = writeZ_2D(cfile,allData,ForPC);
%
%   write contents of cell array allData to file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%
%   Slightly modified from writeTEZ to allow for TE+TM
%
%   Set optional input ForPC = 1 for binary
%   files compatible with PC = 0 for sun workstation

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

   %  record # 1 in transmitter block: period, mode, number of sites 
   [nSites] = size(allData{j}.siteLoc,1);
   ll =  14;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   count = fwrite(fid,allData{j}.T,'float64');
   status = status + 8*count; 
   count = fwrite(fid,allData{j}.Mode,'char');
   status = status + count; 
   count = fwrite(fid,nSites,'long');
   status = status + 4*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

   %  record # 2 in transmitter block: site locations
   ll = 2*8*nSites;
   count = fwrite(fid,ll,'long');
   status = status + 4*count;
   siteLoc = allData{j}.siteLoc';
   count = fwrite(fid,siteLoc,'float64');
   status = status + 8*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

   %  record # 3 in transmitter block: data
   ll = 2*8*nSites;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   Zri = [real(allData{j}.Z) imag(allData{j}.Z)].';
   count = fwrite(fid,Zri,'float64');
   status = status + 8*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

   %  record # 4 in transmitter block: error bars
   %  NOTE: write out same errors for real and imag parts!
   ll = 2*8*nSites;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   %  NOTE: write out same errors for real and imag parts!
   Zri = [real(allData{j}.Zerr) real(allData{j}.Zerr)].';
   count = fwrite(fid,Zri,'float64');
   status = status + 8*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

end
fclose(fid)
