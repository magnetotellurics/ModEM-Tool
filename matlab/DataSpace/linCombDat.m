function [d] = linCombDat(c1,d1,c2,d2);
%  Usage : [d] = linCombDat(d1,c1,d2,c2);
%   computes d = c1*d1+c2*d2 for two data vector object
%    so far error bars  are just copied from d1;

d = d1;
nTx = length(d);
for k = 1:nTx
    d{k}.Cmplx = d1{k}.Cmplx | d2{k}.Cmplx;
    d{k}.Z = c1*d1{k}.Z + c2*d2{k}.Z;
end
