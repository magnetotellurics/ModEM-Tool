function [Jtilde] = NormalizeJ(J,CmHalf);
%   Normalizes Jacobian structure, outputing
%   J_tilde = C_m^{1/2} J C_d^{-1/2} and d_tilde = C_d{-1/2} d
%
%   Usage:  [Jtilde] = NormalizeJ(J,CmHalf);

if nargin == 1
   CM_MULT = 0;
else
   CM_MULT = 1;
end

CdInv = InvErrCov(J{1});
Jtilde= cell(3,1);
Jtilde{3} = J{3};
Jtilde{1} = CdInv_d(J{1});
M = length(J{2});
Jtilde{2} = cell(1,M);
for k = 1:M
   if CM_MULT
      m = CovMult(J{2}{k},CmHalf,1);
   else
      m = J{2}{k};
   end
   m.v = m.v*CdInv(k);
   Jtilde{2}{k} = m;
end
