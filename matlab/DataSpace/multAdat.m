function [dOut] = multAdat(A,d);
%   Matrix vector multiply  required for DCG:
%       [Cd^{-1/2] J Cm J^T Cd^{-1/2} + mu I ] dIn = dOut
%
% Usage:  [dOut] = multA_d(A,d);
%   A = structure with fields Cm, m, mu

dOut = dIn;
nTx = length(dOut);
%   divide by data error standard dev
for k = 1:nTx
   dOut{k}.Z = dIn{k}.Z./dIn{k}.Zerr;
end
m_temp = JmultT(A.m,dOut);
m_temp = CovMult(m_temp,A.Cm,2);
[dOut] = J_times_m(A.m,dIn);
%   divide by data error standard dev, add mu*dIn.
for k = 1:nTx
   dOut{k}.Z = dOut{k}.Z./dOut{k}.Zerr;
   dOut{k}.Z = dOut{k}.z + A.mu*dIn{k}.Z;
end
