function [dOut,dErr] = dataVecToVec(dIn);
%  Makes a standard real vector out of a complex impedance
%  data vector object; if called with two arguments
%  Also returns data standard error as a vector
%
%  Usage: [dOut,dErr] = dataVecToVec(dIn);

nTx = length(dIn);
dOut = [];
for k = 1:nTx
   if dIn{k}.Cmplx
      Zr = real(dIn{k}.Z);
      Zi = imag(dIn{k}.Z);
      Z = zeros(2*length(Zr),1);
      Z(1:2:end) = Zr;
      Z(2:2:end) = Zi;
   else
      Z = real(dIn{k}.Z);
   end
   dOut = [dOut; Z];
end
if nargout ==2
   dErr = [];
   for k = 1:nTx
      nsta = length(dIn{k}.Zerr);
      if dIn{k}.Cmplx
         Zerr = zeros(1,2*nsta);
         Zerr(1:2:end-1) = dIn{k}.Zerr;
         Zerr(2:2:end) = dIn{k}.Zerr;
      else
         Zerr = dIn{k}.Zerr;
      end
      dErr = [dErr; Zerr'];
   end
end
