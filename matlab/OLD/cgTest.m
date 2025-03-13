dh=resid;
NORMALIZE = 1;
g = JTfreq(dh,NORMALIZE);
g = g(1,:);
g0 = g;
h = g;
P = zeros(size(g{1}.v));

%   test symmetry of multA
x = h;
y = h;
for k = 1:length(h)
   x{k}.v = randn(size(x{k}.v));
   y{k}.v = randn(size(x{k}.v));
end
Ax = multA(x);
Ay = multA(y);
yAx = ipMod(y,Ax);
xAy = ipMod(x,Ay);

%   cgIter
g_n = g;
h_n = h;

Ah_n = multA(h_n);
lambda_n = ipMod(g_n,g_n)./ipMod(g_n,Ah_n);
lambda_n = -lambda_n;
g_n1 = linCombMod(ones(size(lambda_n)),g_n,lambda_n,Ah_n);
gamma_n = ipMod(g_n1,Ah_n)./ipMod(h_n,Ah_n);
gamma_n = -gamma_n;
h_n1 = linCombMod(ones(size(gamma_n)),g_n1,gamma_n,h_n);
