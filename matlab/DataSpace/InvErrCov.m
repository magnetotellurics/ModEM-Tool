function CdInv = InvErrCov(dIn);
% Makes a standard real vector containing the inverse of the data
% error standard deviation from and input impedance data vector object
% 
% Usage:  CdInv = InvErrCov(dIn);

nTx = length(dIn);
CdInv = [];
for k = 1:nTx
   Zr = real(dIn{k}.Zerr);
   if dIn{k}.Cmplx
      % errors for real and imaginary part are assumed equal
      Cd = zeros(length(Zr)*2,1);
      Cd(1:2:end) = Zr;
      Cd(2:2:end) = Zr;
   else
      Cd = Zr;
   end
   CdInv = [CdInv; Cd];
end
CdInv = 1./CdInv;
