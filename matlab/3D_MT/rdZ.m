function [Z,loc,T] = rdZ(Zfile)

%  loads all impedances, and if available,
%   vertical field TFs from impedance file Zfile
%    (new format)
% USAGE:  [Z,loc,T] = rdZ(Zfile);
% Zfile = input file name
%  Z(2,2,nSites,nPer) = array of impedances 
% (plus vertical field TFs if in file ... then Z(3,2,nSites,nPer))
%  loc(3,nSites) = array of site locations

%  open file for reading
fid = fopen(Zfile,'r');

% Overall file header: nPer, nSites, ifBzTF
lrec = fread(fid,1,'int32');
if lrec == 32
   %  read version number
   fVer = fread(fid,20,'char');
end
N = fread(fid,3,'int32');
lrec = fread(fid,1,'int32');
nPer = N(1);
nSites = N(2);
ifBzTF = N(3);
if ifBzTF
   NTF = 3;
else
   NTF = 2;
end
Z = zeros(NTF,2,nSites,nPer)+i*zeros(NTF,2,nSites,nPer);
T = zeros(nPer,1);

% site locations
lrec = fread(fid,1,'int32');
loc = fread(fid,[3,nSites],'float64');
lrec = fread(fid,1,'int32');

for iPer = 1:nPer
   %  frequency band header
   lrec = fread(fid,1,'int32');
   omega = fread(fid,1,'float64');
   T(iPer) = 2*pi/omega;
   iFreq = fread(fid,1,'int32');
   lrec = fread(fid,1,'int32');
   % impedances for one band
   lrec = fread(fid,1,'int32');
   temp = fread(fid,[4*NTF,nSites],'float64');
   temp = reshape(temp,[2,NTF,2,nSites]);
   temp(:,1:2,:,:) = temp(:,1:2,:,:)/1000;
   Z(:,:,:,iPer) = (temp(1,:,:,:)+i*temp(2,:,:,:));
   lrec = fread(fid,1,'int32');
end 
