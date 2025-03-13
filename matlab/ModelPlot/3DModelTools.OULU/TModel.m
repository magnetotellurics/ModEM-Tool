classdef TModel < handle
% Egbert, Smirnov, Cherevatova 2013
%
%

properties
  Lat, Long, Z   % number of grid cells in x,y,z directions
  Nx, Ny, Nz     % cell dimensions: x,y,z-direction
  sigma          % 3D array of sigma coductivities
end

methods
%*******************************************************************
function obj = TModel(Filename)
  %   class constructor ... simple
end
%*******************************************************************
%*******************************************************************

function PlotPlainView(obj, Nz)
  %  plot plain view of sigma for Nz depth
  if Nz < obj.Nz+1
    figure;
    contourf(1:obj.Nx, 1:obj.Ny,  squeeze(obj.sigma(:,:,Nz))');
  else
    disp('Specify another depth layer to plot. Nz is out of range');  
  end
  
end

% Plot functions



end     % methods
end    % classdef