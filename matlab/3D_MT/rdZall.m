function [Z,loc,T] = rdZall(Zfile);
%  Loads all impedances, and if available, vertical field TFs
%   from impedance file Zfile
%  Usage:   [Z,loc,T] = rdZall(Zfile);

[Ztemp,loc,T1,nPer] = rdZmodel(Zfile, 1);
%   divide by 1000 to get into standard units (mV/km*nT)
Ztemp(1:2,:,:) = Ztemp(1:2,:,:)/1000;;
[nTF,npol,nSites] = size(Ztemp);
Z = zeros(nTF,npol,nSites,nPer);
Z = Z+i*Z;
Z(:,:,:,1) = Ztemp;
T = zeros(nPer);
T(1) = T1;
for k = 2:nPer
   [Ztemp,loc,T1,nPer] = rdZmodel(Zfile, k);
   %   divide by 1000 to get into standard units (mV/km*nT)
   Ztemp(1:2,:,:) = Ztemp(1:2,:,:)/1000;
   Z(:,:,:,k) = Ztemp;
   T(k) = T1;
end
