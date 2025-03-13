function [dOut] = J_times_m(J,m)
%  Usage : [dOut] = J_times_m(J,m);
%    multiply m by sensitivity matrix J;

dOut = J{1};
nTx = length(dOut);
kk = 0;
for k = 1:nTx
    nSite = length(dOut{k}.siteLoc);
    for j = 1:nSite
       kk = kk + 1;
       Zr = sum(sum(J{2}{kk}.v .* m.v));
       if dOut{k}.Cmplx
          kk = kk + 1;
          Zi = sum(sum(J{2}{kk}.v .* m.v));
          dOut{k}.Z(j)  = Zr+i*Zi;
       else
          dOut{k}.Z(j)  = Zr;
       end
    end
end
