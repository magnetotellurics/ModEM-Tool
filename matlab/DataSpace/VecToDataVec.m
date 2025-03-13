function [dOut] = VecToDataVec(dIn,dTemplate);
%  Usage: [dOut] = VecToDataVec(dIn,dTemplate);
%  converts a standard real vector to a complex impedance 
%   data vector object
%  Inputs:  dIn is a real vector
%           dTemplate is a template data vector structure
nTx = length(dTemplate);
dOut = dTemplate;
ii = 1;
for k = 1:nTx
   nSites = length(dTemplate{k}.siteLoc);
   for j = 1:nSites
      if dTemplate{k}.Cmplx
         dOut{k}.Z(j) = dIn(ii)+i*dIn(ii+1);
         ii = ii + 2;
      else
         dOut{k}.Z(j) = dIn(ii);
         ii = ii + 1;
      end
   end
end
