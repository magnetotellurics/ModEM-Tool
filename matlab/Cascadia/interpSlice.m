function [XI,YI,ZI,C] = interpSlice(Res,ll1,ll2,zLims,np,nz)

%  ll1, ll2 are (lat,lon) for two ends of profile
%  zLims gives depth range for section

XI = zeros(np,nz);
YI = XI;
ZI = YI;
XI = [0:1/(np-1):1]'*ones(1,nz)*(ll2(1)-ll1(1))+ll1(1);
YI = [0:1/(np-1):1]'*ones(1,nz)*(ll2(2)-ll1(2))+ll1(2);
ZI = ones(np,1)*[0:1/(nz-1):1]*(zLims(2)-zLims(1))+zLims(1);

[nx,ny] = size(Res.lon);
nz = length(Res.z);
X = Res.X;
Y = Res.Y;
Z = Res.Z;
C = interp3(X,Y,Z,Res.log10_rho,YI,XI,ZI);
