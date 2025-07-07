function [cond] = readCond_3D(cfile,format)
% readCond_3D - Read condiuctivity file of either WS3D (ModEM Default) or Mackie3D format
% Usage:  [cond] = readCond_3D(cfile, format)
%
% Read from either a Mackie3D file or a Weerachai Siripunvaraporn's "0"  (WS3D) format.
% 
% Input Arguments
%  cfile - File name to read from
%   string value
%  format - Model format either: 1 = Mackie format, or 2 = WS3D format
%   integer value 
%
% Output Argument
%  cond - Struct that describes the file:
%   struct value
%    cond.paramType - Resistivity domain either LINEAR, LOGE, or LOG10
%      string value
%    cond.v - Conductivity values
%      double array
%    cond.AirCond - Value in file of air conductivity (or default: cond.paramType(1e-10))
%      double
%    cond.grid - Struct that describes the grid
%      struct value
%        grid.dx, grid.dy, grid.dz -  Grid distances in x, y and z directions in units (defualt kilometers)
%          double array
%        grid.nx, grid.ny - Number of grid cells in x and y direction
%           integer value
%        grid.NzEarth - Number of layers in the z direction (without air)
%           integer value
%        grid.NzAir - Number of air layers
%           integer value
%        grid.origin - Origin of the model
%           double array [3x1]
%        grid.rotation - Orientation/rotation of the model
%           double value
%        grid.units - Model file units - default 'km'
%           string value
%
% See also writeCond_3D, read_mackie3d_model, read_WS3d_model, write_WS3d_model, write_mackie3d_model

if nargin < 2
    format = 1;
end

if format == 1
    [dx,dy,dz,rho,nzAir,type,origin,rotation] = read_mackie3d_model(cfile);
elseif format == 2
    [dx,dy,dz,rho,nzAir,type,origin,rotation] = read_WS3d_model(cfile);
else
    error('Unknown format');
end

grid.dx = dx/1000;
grid.dy = dy/1000;
grid.dz = dz/1000;
grid.Nx = length(dx);
grid.Ny = length(dy);
grid.NzEarth = length(dz);
grid.NzAir = nzAir;
grid.origin = origin/1000;
grid.rotation = rotation;
grid.units = 'km';

if findstr(type,'LOGE')
    cond.paramType = 'LOGE';
    cond.v = - rho;
    cond.AirCond = log(1e-10);
elseif findstr(type,'LOG10')
    cond.paramType = 'LOG10';
    cond.v = - rho;
    cond.AirCond = log10(1e-10);
else
    cond.paramType = 'LINEAR';
    cond.v = 1./rho;
    cond.AirCond = 1e-10;
end
cond.grid = grid;
