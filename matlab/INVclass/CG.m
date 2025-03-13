function [x,CGiter] = CG(A,b,CGiter)
%
%   Solves the symetric system Ax = b by conjugate gradients
%
%   Usage: [x,CGiter] = CG(A,b,CGiter);
%
%  A is a "representer matrix" object
%  b is a data vector object

%  CGiter is an iterControl object

x = b;
x = zeroVec(x);
r = b;
bNorm = sqrt(dot(b,b));
iter = 0;
DONE = 0;
while ~DONE
   iter = iter + 1;
   if iter == 1
      rho = dot(r,r);
      p = r;
   else
      rho = dot(r,r);
      beta = rho/rhoOld;
      p = r + beta*p;
   end
   q = A*p;
   alpha = rho/dot(p,q);
   x = x + alpha*p;
   r = r - alpha*q;
   rhoOld = rho;
   CGiter.niter = iter;
   CGiter.relErr(iter) = sqrt(rho)/bNorm;
   CGiter.UserData{iter} = struct('x',x,'p',p,'q',q);
   DONE = checkConvergence(CGiter);  
end

