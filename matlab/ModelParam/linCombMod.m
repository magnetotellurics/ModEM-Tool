function [m] = linCombMod(c1,m1,c2,m2)
%  Computes linear combination of cell arrays of model parameter
%  structures
%  NOT FOR A SINGLE MODEL PARAMETER
%  
%  Usage: [m] = linCombMod(c1,m1,c2,m2);
%     m1 and m2 are input cell arrays of model parameters; c1 and c2
%     are (in general complex) coefficients (vectors)

m = m1;
nMod = length(m);
for k = 1:nMod
   m{k}.v = c1(k)*m1{k}.v + c2(k)*m2{k}.v;
end
