function ip = dotDat(d1,d2)
%   real inner product of two (complex impedance)
%    data objects.  Note that vectors are stored as
%    complex objects, but the inner product corresponds
%    to the real inner product of the corresponding real
%    vectors of real and imaginary parts.
%  Usage :  ip = dotDat(d1,d2);

nTx = length(d1);
ip= 0.;
for k = 1:nTx
    ip = ip + real(d1{k}.Z'*d2{k}.Z);
end