function [allData] = readTEZ(cfile)
%  Usage: [allData] = readTEZ(cfile)
%   reads data from impedance file, returns 
%    as cell array allData
%   There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)

fid = fopen(cfile,'r');
%  record # 1: number of periods 
ll = fread(fid,1,'long');
nPer = fread(fid,1,'long');
ll = fread(fid,1,'long');

for j = 1:nPer
%   for each period ....

   %  record # 1 in period block: period, number of sites 
   ll = fread(fid,1,'long');
   T = fread(fid,1,'float64');
   nSites = fread(fid,1,'long');
   ll = fread(fid,1,'long');

   %  record # 2 in period block: site locations
   ll = fread(fid,1,'long');
   siteLoc = fread(fid,[2,nSites],'float64');
   ll = fread(fid,1,'long');

   %  record # 3 in period block: data
   ll = fread(fid,1,'long');
   Zri = fread(fid,[2,nSites],'float64');
   ll = fread(fid,1,'long');
   Z = Zri(1,:)+i*Zri(2,:);

   %  record # 4 in period block: error bars
   ll = fread(fid,1,'long');
   Zerr = fread(fid,[2,nSites],'float64');
   ll = fread(fid,1,'long');
   Zerr = [Zerr(1,:)];

   allData{j} = struct('T',T,'siteLoc',siteLoc,'Z',Z,'Zerr',Zerr);
end
fclose(fid);
