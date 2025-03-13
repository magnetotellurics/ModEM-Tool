function B = TEmag(E,Dz,omega)
%  Usage:  B = TEmag(E)
%    compute magnetic induction on edges, given electric 
%    fields on nodes 
%  (still need to check sign ...)
B = E(:,2:end)-E(:,1:end-1);
B = B*diag(1./Dz);
B = B/(i*omega);
