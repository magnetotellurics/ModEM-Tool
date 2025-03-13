function [dOut] = J_times_m_MTX(J,mIn)
%  Usage : [dOut] = J_times_m_MTX(J,mIn);
% multiply cell array of model paramters (one for
%   each frequency) by sensitivity matrix 
%   (stored as cell array of structures)

nTx = length(J{1});
dOut = J{1};
kk = 0;
for k = 1:nTx
   nSite = length(dOut{k}.siteLoc);
   for j = 1:nSite
      kk = kk + 1;
      Zr = sum(sum(J{2}{kk}.v .* mIn{k}.v));
      if dOut{k}.Cmplx
         kk = kk + 1;
         Zi = sum(sum(J{2}{kk}.v .* mIn{k}.v));
         dOut{k}.Z(j)  = Zr+i*Zi;
      else
         dOut{k}.Z(j)  = Zr;
      end
   end
end
