function [m] = linCombMod1(c1,m1,c2,m2)
%
%   computes linear combination of model parameters
%                 m = c1*m1+c2*m2
%   This version is for a single model parameter
%
%  Usage: [m] = linCombMod(c1,m1,c2,m2);

m = m1;
m.v = c1*m1.v + c2*m2.v;
