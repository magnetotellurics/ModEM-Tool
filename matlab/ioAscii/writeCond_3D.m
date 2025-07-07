function [status] = writeCond_3D(cfile,cond,format)
% writeCond_3D - Create a conductivity file in either WS3D or Mackie3D format
% Usage:  [status] = writeCond_3D(cfile, cond, format)
%
% This function is best used with readCond_3D.
%
% Input Arguments
%   cfile - File name to write too
%     string value
%   cond - Struct that contains information to write grid and conductivities
%     struct value
%       cond.paramType - Resistivity domain either LINEAR, LOGE, or LOG10
%         string value
%       cond.v - Conductivity values
%         double array
%       cond.AirCond - Value in file of air conductivity (or default: cond.paramType(1e-10))
%         double
%       cond.grid - Struct that describes the grid
%         struct value
%           grid.dx, grid.dy, grid.dz -  Grid distances in x, y and z directions in units (defualt kilometers)
%             double array
%           grid.nx, grid.ny - Number of grid cells in x and y direction
%             integer value
%           grid.NzEarth - Number of layers in the z direction (without air)
%             integer value
%           grid.NzAir - Number of air layers
%             integer value
%           grid.origin - Origin of the modeljk
%             double array [3x1]
%           grid.rotation - Orientation/rotation of the model
%             double value
%           grid.units - Model file units - default 'km'
%             string value
% format - Model format either: 1 = Mackie's format (default) or 2 = WS3D format 
%
% Output Arguments
%   status - Number of bytes written - Status returned by MatLab write routines
%    integer value
%
% See also readCond_3D, read_mackie3d_model, read_WS3d_model, write_WS3d_model, write_mackie3d_model

if nargin < 3
    format = 1;
end

nzAir = cond.grid.NzAir;
dx = cond.grid.dx;
dy = cond.grid.dy;
dz = cond.grid.dz;

if isfield(cond.grid,'origin')
    origin = cond.grid.origin;
else
    origin = [0 0 0];
end

if isfield(cond.grid,'units')
    if strcmp(cond.grid.units,'km')
        % convert everything to meters!
        dx = 1000*dx;
        dy = 1000*dy;
        dz = 1000*dz;
        origin = 1000*origin;
    end
end        

if isfield(cond.grid,'rotation')
    rotation = cond.grid.rotation;
else
    rotation = 0;
end

if findstr(cond.paramType,'LOGE')
    type = 'LOGE';
    rho = - cond.v;
elseif findstr(cond.paramType,'LOG10')
    type = 'LOG10';
    rho = - cond.v;
else
    type = 'LINEAR';
    rho = 1./(cond.v);
end

if format == 1
    status = write_mackie3d_model(cfile,dx,dy,dz,rho,nzAir,type,origin,rotation);
elseif format == 2
    status = write_WS3d_model(cfile,dx,dy,dz,rho,nzAir,type,origin,rotation);
else
    error('Unknown format');
end
