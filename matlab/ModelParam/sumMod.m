function [m] = sumMod(m1,m2)
%  Usage: [m] = sumMod(m1,m2);
%   computes sum of model parameters m1 and m2

m = m1;
m.v = m1.v+m2.v;
