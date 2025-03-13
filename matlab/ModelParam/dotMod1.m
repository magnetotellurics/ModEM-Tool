function [ip] = dotMod1(m1,m2);
%  Variant on dotMod for single model parameters m1, m2
%
%  Usage:  [ip] = dotMod(m1,m2);

ip = sum(sum(m1.v .* m2.v));
