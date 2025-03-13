function [mOut] = InitHalfSpace(mIn,sig)
%   Usage: [mOut] = InitHalfSpace(mIn,sig);
%   Initialize a model parameter mOut as half space
%     mIn = template for model parameter
%     sig = conductivity (NOT log conductivity, 
%		or resistivity) for half space
%
%  set prior, and starting conductivity (uniform half space)
mOut = mIn;
if strcmp(mIn.paramType,'LOGE')  
   mOut.v = log(sig)*ones(size(mIn.v));
else
   mOut.v = sig*ones(size(mIn.v));
end
