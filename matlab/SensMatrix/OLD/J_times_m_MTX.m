function [dOut] = J_times_m_MTX(mIn,d)
%  Usage : [dOut] = J_times_m_MTX(mIn,d);
% multiply cell array of model paramters (one for
%   each frequency) by sensitivity matrix 
%   (stored as cell array of structures)
%   d is template for data vector
%   sensitivity matrix J is global

global J

nTx = length(d);
dOut = d;
kk = 0;
for k = 1:nTx
    nSite = length(d{k}.siteLoc);
    for j = 1:nSite
       kk = kk + 1;
       Zr = sum(sum(J{kk}.v .* mIn{k}.v));
       kk = kk + 1;
       Zi = sum(sum(J{kk}.v .* mIn{k}.v));
       dOut{k}.Z(j)  = Zr+i*Zi;
    end
end
