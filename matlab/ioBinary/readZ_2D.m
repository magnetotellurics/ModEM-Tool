function [allData] = readZ_2D(cfile,ForPC)
%
%  Usage: [allData] = readZ_2D(cfile)
%         [allData] = readZ_2D(cfile,ForPC)
%
%   reads data from impedance file, returns 
%    as cell array allData
%   There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   slightly modified from readTEZ to allow for TE+TM 
%
%   ForPC is an optional argument, true (1) to produce
%   files which can be read on PC (linux/windows)

if nargin == 1
   fid = fopen(cfile,'r','n');
else
   if ForPC
      fid = fopen(cfile,'r','l');
   else
      fid = fopen(cfile,'r','b');
   end
end

%  record # 1: number of transmitters 
ll = fread(fid,1,'long');
nTx = fread(fid,1,'long');
ll = fread(fid,1,'long');

for j = 1:nTx
%   for each period/mode ....

   %  record # 1 in transmitter block: period, mode, number of sites 
   ll = fread(fid,1,'long');
   T = fread(fid,1,'float64');
   MODE = char(fread(fid,2,'char'));
   MODE = MODE';
   nSites = fread(fid,1,'long');
   ll = fread(fid,1,'long');

   %  record # 2 in transmitter block: site locations
   ll = fread(fid,1,'long');
   siteLoc = fread(fid,[2,nSites],'float64');
   ll = fread(fid,1,'long');

   %  record # 3 in transmitter block: data
   ll = fread(fid,1,'long');
   Zri = fread(fid,[2,nSites],'float64');
   ll = fread(fid,1,'long');
   Z = Zri(1,:)+i*Zri(2,:);

   %  record # 4 in transmitter block: error bars
   ll = fread(fid,1,'long');
   Zerr = fread(fid,[2,nSites],'float64');
   ll = fread(fid,1,'long');
   Zerr = [Zerr(1,:)];

   allData{j} = struct('T',T,'Cmplx',1,...
	'Mode',MODE,'siteLoc',siteLoc','Z',Z.','Zerr',Zerr.');
end
fclose(fid);
