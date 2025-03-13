function [d] = minusDat(d1,d2)
%  Implements subtraction for data vector objects
%   so far error bars  are just copied from d1;
%
%  Usage : [d] = minusDat(d1,d2);

nTx1 = length(d1);
nTx2 = length(d2);
if nTx1 ~= nTx2
   fprintf(1,'%s\n','Error: data vector objects not compatible')
   return
end

d = d1;
for k = 1:nTx1
    d{k}.Cmplx = d1{k}.Cmplx | d2{k}.Cmplx;
    d{k}.Z = d1{k}.Z - d2{k}.Z;
end
