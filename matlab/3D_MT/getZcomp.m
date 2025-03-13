function [C,xLab,yLab,X,Y,AxMode,lev] = getZcomp(P)
%  Gets the appropriate slice/component of impedance
%  array Z; all information is stored in data structure
%  P, which is attached as "userdata" to the figure.
% Usage [C,xLab,yLab,X,Y,AxMode,lev] = getZcomp(P);

if(P.Slice == 'T')
  C = squeeze(P.Z(P.Comp,P.Mode,P.i1:P.i2,P.j1:P.j2,P.Np));
  X =  P.Y(P.j1:P.j2);
  Y = P.X(P.i1:P.i2);
  lev = P.T;
  xLab = 'Y';
  yLab = 'X';
  AxMode = 'xy';
elseif (P.Slice == 'X')
  C = squeeze(P.Z(P.Comp,P.Mode,P.Np,P.j1:P.j2,P.k1:P.k2)).';
  X =  P.Y(P.j1:P.j2);
  Y = log10(P.T(P.k1:P.k2));
  lev = P.X;
  xLab = 'Y';
  yLab = 'log_{10}T';
  AxMode = 'ij';
elseif(P.Slice == 'Y')
  C = squeeze(P.Z(P.Comp,P.Mode,P.i1:P.i2,P.Np,P.k1:P.k2)).';
  X =  P.X(P.i1:P.i2);
  Y = log10(P.T(P.k1:P.k2));
  lev = P.Y;
  xLab = 'X';
  yLab = 'log_{10}T';
  AxMode = 'ij';
else
  %  write error msg, return
end
