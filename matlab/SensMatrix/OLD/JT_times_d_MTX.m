function [g] = JT_times_d_MTX(dIn,mIn,NORMALIZE);
%  variant on JT_times_d ... does not sum over
%   frequencies, and also returns "imaginary part"
%   output is a cell array of 2*nTx model parameters
%    mIn is a template for the output model parameters
%   Sensitivity J is global (already set)

global J

if nargin < 3
  NORMALIZE = 0;
end

nTx = length(dIn);
g = cell(2,nTx);
kk = 0;
for k = 1:nTx
   g{1,k} = mIn; g{1,k}.v = 0;
   g{2,k} = g{1,k}; 
   nSite = length(dIn{k}.siteLoc);
   for j = 1:nSite
      if NORMALIZE
         temp = (dIn{k}.Z(j)/(dIn{k}.Zerr(j))^2);
      else
         temp = dIn{k}.Z(j);
      end
      kk = kk + 1;
      g{1,k}.v = g{1,k}.v + real(temp)*J{kk}.v;
      g{2,k}.v = g{2,k}.v + imag(temp)*J{kk}.v;
      kk = kk + 1;
      g{1,k}.v = g{1,k}.v + imag(temp)*J{kk}.v;
      g{2,k}.v = g{2,k}.v - real(temp)*J{kk}.v;
   end
end
