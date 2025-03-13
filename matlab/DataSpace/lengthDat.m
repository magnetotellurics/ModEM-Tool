function [n] = lengthDat(d);
%  Total length of a data vector object
%
%  Usage [n] = lengthDat(d);

nTx = length(d);
n = 0;
for k = 1:nTx
   if d{k}.Cmplx
      n = n + 2*length(d{k}.siteLoc);
   else
      n = n + length(d{k}.siteLoc);
   end
end
