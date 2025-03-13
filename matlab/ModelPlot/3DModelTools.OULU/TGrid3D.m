classdef TGrid3D < handle
% Egbert,Smirnov,Cherevatova 2013
% Railo 2014
%
% base class for storing 2D MT data; follows structure used in
% original non-object-oriented "dataSpace" routines; idea is to
% quickly develop a data-space object that can later be replaced by
% something more logical

properties
  Nx, Ny, Nz    %	number of grid cells in x,y,z directions
  Nza           %	number of air layers
  Dx, Dy, Dz    %	cell dimensions: x,y,z-direction (vector, not scalar!)
  origin        %	origin at UTM coordinates (xUTM, yUTM)
end

methods
%*******************************************************************
function obj = TGrid3D(Dx,Dy,Dz,Nza,varargin)
  %   class constructor ... simple
  if nargin >= 4
    obj.Nx = length(Dx);
    obj.Ny = length(Dy);
    obj.Nz = length(Dz);
    obj.Nza = Nza;
    obj.Dx = Dx;
    obj.Dy = Dy;
    obj.Dz = Dz;
    if nargin == 5
        obj.origin = varargin{1};
    end
  end
end

function addOrigin(obj,x,y,zone)
    obj.origin.x = x;
    obj.origin.y = y;
    obj.origin.zone = zone;
end
 %*******************************************************************
%*******************************************************************
function [nx,ny,nz] = setLimits(obj,type)
  %   set arraay limits for components of types
  %                     cell, node,
  %                     xedge, yede, zedge
  %                     xface, yface, zface
   %   should probably add range checks
  switch lower(type)
    case 'cell'
        nx = obj.Nx;
        ny = obj.Ny;
        nz = obj.Nz;
    case 'node'
        nx = obj.Nx+1;
        ny = obj.Ny+1;
        nz = obj.Nz+1;
    case 'xedge'
        nx = obj.Nx;
        ny = obj.Ny+1;
        nz = obj.Nz+1;
    case 'yedge'
        nx = obj.Nx+1;
        ny = obj.Ny;
        nz = obj.Nz+1;
    case 'zedge'
        nx = obj.Nx+1;
        ny = obj.Ny+1;
        nz = obj.Nz;
    case 'xface'
        nx = obj.Nx+1;
        ny = obj.Ny;
        nz = obj.Nz;
    case 'yface'
        nx = obj.Nx;
        ny = obj.Ny+1;
        nz = obj.Nz;
    case 'zface'
        nx = obj.Nx;
        ny = obj.Ny;
        nz = obj.Nz+1;
  end
end
%*******************************************************************
function [I,J,K] = gridIndex(obj,index,type)
  %  cells, nodes, edges, faces are enumerated as column vector
  %   elements in the standard way, consistent with 
  %    matlab reshape(X,Nx,Ny,Nz), where X is an array of
  %  dimension X(Nx,Ny,Nz).   Given array of integer indices of 
  %   vector elements this function, computes the
  %  corresponding i,j,k for an array of given type
  %   allowabel types : cell, node,
  %                     xedge, yede, zedge
  %                     xface, yface, zface
  
  [nx,ny,nz] = obj.setLimits(type);
  I = mod(index,nx);
  I(I==0) = nx;
  J  = mod(ceil(index/nx),ny);
  J(J==0) = ny;
  K = ceil(index/(nx*ny));
end
%*******************************************************************
function [index] = vectorIndex(obj,I,J,K,type)
  %  cells, nodes, edges, faces are enumerated as column vector
  %   elements in the standard way, consistent with 
  %    matlab reshape(X,Nx,Ny,Nz), where X is an array of
  %  dimension X(Nx,Ny,Nz).   Given array of integer indices of 
  %   vector elements this function, computes the
  %  corresponding i,j,k for an array of given type
  %   allowabel types : cell, node,
  %                     xedge, yede, zedge
  %                     xface, yface, zface
  
  [nx,ny,nz] = obj.setLimits(type);
  index = (K-1)*nx*ny+(J-1)*nx+I;
end %vectorIndex

end     % methods
end    % classdef