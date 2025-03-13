function [rho,phi,T,y,MODE,rho2,phi2,T2,y2] = mkPseud(Z)
%   Usage:  [rho,phi] = mkPseud(Z);
%   make rho/phi pseudosection for one mode from
%   general 2D impedance data structure
%   add capability to sort out modes later!
%    to start assume just one mode + same # of sites for all
%    periods
nT = length(Z);
nSites = length(Z{1}.Z);
nTE = 0;
nTM = 0;
for k = 1:nT
   if strcmp(Z{k}.Mode,'TE')
      nTE = nTE+1;
   end 
   if strcmp(Z{k}.Mode,'TM')
      nTM = nTM+1;
   end
end
if nTE > 0
   if nTM > 0
      MODE = 'JT';
   else
      MODE = 'TE';
   end
else
   if nTM > 0
      MODE = 'TM';
   else
      error('Error: modes not specified correctly in mkPseud')
   end
end

kTE  = 0;
kTM = 0;
if strcmp(MODE,'JT')
   rho = zeros(nTE,nSites);
   phi = zeros(nTE,nSites);
   T = zeros(nTE,1);
   rho2 = zeros(nTM,nSites);
   phi2 = zeros(nTM,nSites);
   T2 = zeros(nTM,1);
   for k = 1:nT
      if strcmp(Z{k}.Mode,'TE')
         kTE = kTE + 1;
         T(kTE) = Z{k}.T;
         rho(kTE,:) = Z{k}.T*(abs(Z{k}.Z/1e+3).^2)/5;
         phi(kTE,:) = atan2(-imag(Z{k}.Z),real(Z{k}.Z))*180/pi;
         y = Z{k}.siteLoc(:,1);
      else
         kTM = kTM+1;
         T2(kTM) = Z{k}.T;
         rho2(kTM,:) = Z{k}.T*(abs(Z{k}.Z/1e+3).^2)/5;
         %    switch to first quadrant ... assuming TM mode is in third
         %    quadrant!
         phi2(kTM,:) = atan2(imag(Z{k}.Z),-real(Z{k}.Z))*180/pi;
         y2 = Z{k}.siteLoc(:,1);
      end
   end
   rho = rho(1:kTE,:);
   phi = phi(1:kTE,:);
else
   %   just one mode
   y = Z{1}.siteLoc(:,1);
   nT = max(nTE,nTM);
   rho = zeros(nT,nSites);
   phi = zeros(nT,nSites);
   T = zeros(nTE);
   for k = 1:nT
      T(k) = Z{k}.T;
      rho(k,:) = Z{k}.T*(abs(Z{k}.Z/1e+3).^2)/5;
      if strcmp(MODE,'TE')
          phi(k,:) = atan2(-imag(Z{k}.Z),real(Z{k}.Z))*180/pi;
      else
          phi(k,:) = atan2(imag(Z{k}.Z),-real(Z{k}.Z))*180/pi;
      end
   end
   rho2 = [];
   phi2 = [];
   y2 = [];
   T2 = [];
end
