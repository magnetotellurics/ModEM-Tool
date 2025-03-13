function ip = ipDatNoCov(d1,d2)
%   inner product of two (complex impedance)
%    data objects, ignoring data error covariance.
%  Usage :  ip = ipDat(d1,d2);

nTx = length(d1);
ip= 0.;
for k = 1:nTx
    x = d1{k}.Z;
    y = d2{k}.Z;
    ip = ip + d1{k}.Z'*d2{k}.Z;
end
