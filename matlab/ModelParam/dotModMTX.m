function [ip] = dotModMTX(m1,m2);
% model inner products for cell arrays of model parameters m1, m2
% (does not use model covariance!) NOT FOR SINGLE MODEL PARAMTERS
% 
% Usage:  [ip] = dotMod(m1,m2);
%         ip is output real vector of model space inner products,
%         one for each pair of model parameters in cell arrays m1, m2

nMod = length(m1);
ip = zeros(nMod,1);
for k = 1:nMod
   ip(k) = sum(sum(m1{k}.v .* m2{k}.v));
end
