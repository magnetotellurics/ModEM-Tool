 function [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P)
 % Usage [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
 %  gets a slice of a specified component of mag field H
 %   assumed defined on cell faces

%  now set up component to plot
if(P.Slice == 'Z')
  if(P.Comp == 'X')
     C = squeeze(P.H.x(P.i1:P.i2,P.j1:P.j2,P.Np));
     X =  P.yCenter(P.j1:P.j2);
     Y = P.xEdge(P.i1:P.i2);
     lev = P.zCenter;
  elseif(P.Comp == 'Y')
     C = squeeze(P.H.y(P.i1:P.i2,P.j1:P.j2,P.Np));
     X =  P.yEdge(P.j1:P.j2);
     Y = P.xCenter(P.i1:P.i2);
     lev = P.zCenter;
  else
     C = squeeze(P.H.z(P.i1:P.i2,P.j1:P.j2,P.Np));
     X =  P.yCenter(P.j1:P.j2);
     Y = P.xCenter(P.i1:P.i2);
     lev = P.zEdge;
  end
  xLab = 'Y';
  yLab = 'X';
  AxMode = 'xy';
elseif (P.Slice == 'X')
  if(P.Comp == 'X')
     C = squeeze(P.H.x(P.Np,P.j1:P.j2,P.k1:P.k2)).';
     X =  P.yCenter(P.j1:P.j2);
     Y = P.zCenter(P.k1:P.k2);
     lev = P.xEdge;
  elseif(P.Comp == 'Y')
     C = squeeze(P.H.y(P.Np,P.j1:P.j2,P.k1:P.k2)).';
     X =  P.yEdge(P.j1:P.j2);
     Y = P.zCenter(P.k1:P.k2);
     lev = P.xCenter;
  else
     C = squeeze(P.H.z(P.Np,P.j1:P.j2,P.k1:P.k2)).';
     X =  P.yCenter(P.j1:P.j2);
     Y = P.zEdge(P.k1:P.k2);
     lev = P.xCenter;
  end
  xLab = 'Y';
  yLab = 'Z';
  AxMode = 'ij';
elseif(P.Slice == 'Y')
  if(P.Comp == 'X')
     C = squeeze(P.H.x(P.i1:P.i2,P.Np,P.k1:P.k2)).';
     X =  P.xEdge(P.i1:P.i2);
     Y = P.zCenter(P.k1:P.k2);
     lev = P.yCenter;
  elseif(P.Comp == 'Y')
     C = squeeze(P.H.y(P.i1:P.i2,P.Np,P.k1:P.k2)).';
     X =  P.xCenter(P.i1:P.i2);
     Y = P.zCenter(P.k1:P.k2);
     lev = P.yEdge;
  else
     C = squeeze(P.H.z(P.i1:P.i2,P.Np,P.k1:P.k2)).';
     X =  P.xCenter(P.i1:P.i2);
     Y = P.zCenter(P.k1:P.k2);
     lev = P.yEdge;
  end
  xLab = 'X';
  yLab = 'Z';
  AxMode = 'ij';
else
  %  write error msg, return
end
