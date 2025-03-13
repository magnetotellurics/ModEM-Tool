function [dOut] = CdInv_d(dIn,Ntimes);
% Usage: [dOut] = CdInv_d(dIn,Ntimes);
%      Ntimes is 1 to divide by standard deviation
%                2 to divide by variance
%             This is optional, default is 1
%   NOW ALSO modifies error stadandard deviation
%     consistent with rescaling of data values

if nargin == 1
   Ntimes = 1
end
dOut = dIn;
nTx = length(dIn);
for k = 1:nTx
   nSite = length(dIn{k}.siteLoc);
   if Ntimes == 1
      dOut{k}.Z = dIn{k}.Z./dIn{k}.Zerr;
      dOut{k}.Zerr = ones(size(dIn{k}.Zerr));
   else
      dOut{k}.Z = dIn{k}.Z./(dIn{k}.Zerr).^2;
      dOut{k}.Zerr = 1./dIn{k}.Zerr;
   end
end
