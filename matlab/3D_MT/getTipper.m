function [Tx,Ty,X,Y] = getTipper(P)
%  Gets Tx/TY and construct tipper components for plan view

Tx = -squeeze(P.Z(3,1,P.i1:P.i2,P.j1:P.j2,P.Np)).';
Ty = -squeeze(P.Z(3,2,P.i1:P.i2,P.j1:P.j2,P.Np)).';
X =  P.X(P.i1:P.i2);
Y = P.Y(P.j1:P.j2);
nX = length(X); nY = length(Y)
X = ones(nY,1)*X; 
Y = Y'*ones(1,nX); 
