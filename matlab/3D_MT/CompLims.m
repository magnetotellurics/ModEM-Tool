function [Nx1,Ny1,Nz1,Np] = CompLims(P)
%  Usage: [Nx1,Ny1,Nz1,Np] = CompLims(P)

if(P.PlotH)
  if(P.Comp == 'X')
     Nx1 = length(P.xEdge);
     Ny1 = length(P.yCenter);
     Nz1 = length(P.zCenter);
  elseif(P.Comp == 'Y')
     Nx1 = length(P.xCenter);
     Ny1 = length(P.yEdge);
     Nz1 = length(P.zCenter);
  else
     Nx1 = length(P.xCenter);
     Ny1 = length(P.yCenter);
     Nz1 = length(P.zEdge);
  end
else
  if(P.Comp == 'X')
     Nx1 = length(P.xCenter);
     Ny1 = length(P.yEdge);
     Nz1 = length(P.zEdge);
  elseif(P.Comp == 'Y')
     Nx1 = length(P.xEdge);
     Ny1 = length(P.yCenter);
     Nz1 = length(P.zEdge);
  else
     Nx1 = length(P.xEdge);
     Ny1 = length(P.yEdge);
     Nz1 = length(P.zCenter);
  end
end

  if(P.Slice == 'X')
     Np = min(P.Np,Nx1);
  elseif(P.Slice == 'Y')
     Np = min(P.Np,Ny1)
  else
     Np = min(P.Np,Nz1)
  end
