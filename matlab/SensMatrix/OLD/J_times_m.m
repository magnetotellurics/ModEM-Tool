function [dOut] = J_times_m(m,d)
%  Usage : [J] = J_times_m(m,d);
%    multiply m by sensitivity matrix; d is only used as a template
%    for initializing output vector (stored as cell array of structures)
%    sensitivity matrix is assumed consistent with d, and is global

global J

dOut = d;
nTx = length(dOut);
kk = 0;
for k = 1:nTx
    nSite = length(d{k}.siteLoc);
    for j = 1:nSite
       kk = kk + 1;
       Zr = sum(sum(J{kk}.v .* m.v));
       kk = kk + 1;
       Zi = sum(sum(J{kk}.v .* m.v));
       dOut{k}.Z(j)  = Zr+i*Zi;
    end
end
