function [status] = writeTEZ(cfile,allData)
%  Usage:  [status] = writeTEZ(cfile,allData);
%   write contents of cell array allData to file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)

nPer = length(allData);
fid = fopen(cfile,'w');
status = 0;

%  record # 1: number of periods 
ll = 4;
count = fwrite(fid,ll,'long');
status = status + 4*count; 
count = fwrite(fid,nPer,'long');
status = status + 4*count; 
count = fwrite(fid,ll,'long');
status = status + 4*count; 

for j = 1:nPer
%   for each period ....

   %  record # 1 in period block: period, number of sites 
   [dum,nSites] = size(allData{j}.siteLoc);
   ll =  12;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   count = fwrite(fid,allData{j}.T,'float64');
   status = status + 8*count; 
   count = fwrite(fid,nSites,'long');
   status = status + 4*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

   %  record # 2 in period block: site locations
   ll = 2*8*nSites;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   count = fwrite(fid,allData{j}.siteLoc,'float64');
   status = status + 8*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

   %  record # 3 in period block: data
   ll = 2*8*nSites;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   Zri = [real(allData{j}.Z);imag(allData{j}.Z)];
   count = fwrite(fid,Zri,'float64');
   status = status + 8*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

   %  record # 4 in period block: error bars
   %  NOTE: write out same errors for real and imag parts!
   ll = 2*8*nSites;
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 
   %  NOTE: write out same errors for real and imag parts!
   Zri = [real(allData{j}.Zerr);real(allData{j}.Zerr)];
   count = fwrite(fid,Zri,'float64');
   status = status + 8*count; 
   count = fwrite(fid,ll,'long');
   status = status + 4*count; 

end
fclose(fid)
