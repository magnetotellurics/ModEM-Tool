 function [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P)
 % Usage [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);

%  now set up component to plot
if(P.Slice == 'Z')
  if(P.Comp == 'X')
     C = squeeze(P.E.x(P.i1:P.i2,P.j1:P.j2,P.Np));
     X =  P.yEdge(P.j1:P.j2);
     Y = P.xCenter(P.i1:P.i2);
     lev = P.zEdge;
  elseif(P.Comp == 'Y')
     C = squeeze(P.E.y(P.i1:P.i2,P.j1:P.j2,P.Np));
     X =  P.yCenter(P.j1:P.j2);
     Y = P.xEdge(P.i1:P.i2);
     lev = P.zEdge;
  else
     C = squeeze(P.E.z(P.i1:P.i2,P.j1:P.j2,P.Np));
     X =  P.yEdge(P.j1:P.j2);
     Y = P.xEdge(P.i1:P.i2);
     lev = P.zCenter;
  end
  xLab = 'Y';
  yLab = 'X';
  AxMode = 'xy';
elseif (P.Slice == 'X')
  if(P.Comp == 'X')
     C = squeeze(P.E.x(P.Np,P.j1:P.j2,P.k1:P.k2)).';
     X =  P.yEdge(P.j1:P.j2);
     Y = P.zEdge(P.k1:P.k2);
     lev = P.xCenter;
  elseif(P.Comp == 'Y')
     C = squeeze(P.E.y(P.Np,P.j1:P.j2,P.k1:P.k2)).';
     X =  P.yCenter(P.j1:P.j2);
     Y = P.zEdge(P.k1:P.k2);
     lev = P.xEdge;
  else
     C = squeeze(P.E.z(P.Np,P.j1:P.j2,P.k1:P.k2)).';
     X =  P.yEdge(P.j1:P.j2);
     Y = P.zCenter(P.k1:P.k2);
     lev = P.xEdge;
  end
  xLab = 'Y';
  yLab = 'Z';
  AxMode = 'ij';
elseif(P.Slice == 'Y')
  if(P.Comp == 'X')
     C = squeeze(P.E.x(P.i1:P.i2,P.Np,P.k1:P.k2)).';
     X =  P.xCenter(P.i1:P.i2);
     Y = P.zEdge(P.k1:P.k2);
     lev = P.yEdge;
  elseif(P.Comp == 'Y')
     C = squeeze(P.E.y(P.i1:P.i2,P.Np,P.k1:P.k2)).';
     X =  P.xEdge(P.i1:P.i2);
     Y = P.zEdge(P.k1:P.k2);
     lev = P.yCenter;
  else
     C = squeeze(P.E.z(P.i1:P.i2,P.Np,P.k1:P.k2)).';
     X =  P.xEdge(P.i1:P.i2);
     Y = P.zCenter(P.k1:P.k2);
     lev = P.yEdge;
  end
  xLab = 'X';
  yLab = 'Z';
  AxMode = 'ij';
else
  %  write error msg, return
end
