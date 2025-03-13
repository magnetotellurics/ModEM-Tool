function [mOut] = JT_times_d(dIn,mIn,NORMALIZE)
%  Usage : [mOut] = JT_times_d(dIn,mIn,NORMALIZE);
%    NORMALIZE is opitonal; default is to not normalize
%     by Cd^{-1} before multiplying by J^T
% multiply by transpose of sensitivity matrix (stored as 
%     cell array of structures)
%   mIn is used as a template for the output
%   Sensitivity J is global
global J

if nargin < 3
   NORMALIZE = 0;
end

nTx = length(dIn);
kk = 0;
mOut = mIn;
mOut.v = zeros(size(mOut.v));
for k = 1:nTx
    nSite = length(dIn{k}.siteLoc);
    for j = 1:nSite
       kk = kk + 1;
       if NORMALIZE
          mOut.v = mOut.v + (real(dIn{k}.Z(j))./d{k}.Zerr)*J{kk}.v;
       else
          mOut.v = mOut.v + real(dIn{k}.Z(j))*J{kk}.v;
       end
       kk = kk + 1;
       if NORMALIZE
          mOut.v = mOut.v + (imag(dIn{k}.Z(j))./d{k}.Zerr)*J{kk}.v;
       else
          mOut.v = mOut.v + imag(dIn{k}.Z(j))*J{kk}.v;
       end
    end
end
