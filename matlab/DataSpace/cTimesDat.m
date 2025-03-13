function [dOut] = cTimesDat(alpha,dIn)
%   a very tentative implementation of multiplication
%   by scalars for data vector objects.  Ultimately
%   this might be extended to allow also for multiplication
%   by matrices.  Also now the scalar has to come first,
%   but this could be generalized ... and no error checking so far!
%
%   Usage: [dOut] = cTimesDat(alpha,dIn)

dOut = dIn;
nTx = length(dIn);
for k = 1:nTx
   dOut{k}.Z = alpha*dOut{k}.Z;
end
