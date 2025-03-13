function ip = ipDatFreq(d1,d2)
%   inner product of two (complex impedance)
%    data objects, separately for each frequency
%  using error standard deviation  defined from d1 
%  Usage :  ip = ipDatFreq(d1,d2);
%   returns vector of inner products, one for each frequency

nTx = length(d1);
ip= zeros(nTx,1);
for k = 1:nTx
    x = d1{k}.Z./d1{k}.Zerr;
    y = d2{k}.Z./d1{k}.Zerr;
    ip(k) = x'*y;
end
