function ip = ipDat(d1,d2)
%   inner product of two (complex impedance)
%    data objects, using error standard deviation 
%    defined from d1 
%  Usage :  ip = ipDat(d1,d2);

nTx = length(d1);
ip= 0.;
for k = 1:nTx
    x = d1{k}.Z./d1{k}.Zerr;
    y = d2{k}.Z./d1{k}.Zerr;
    ip = ip + x'*y;
end