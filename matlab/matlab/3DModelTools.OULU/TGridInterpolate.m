classdef TGridInterpolate < handle
% Egbert,Smirnov,Cherevatova 2013   
%
% base class for storing 2D MT data; follows structure used in
% original non-object-oriented "dataSpace" routines; idea is to
% quickly develop a data-space object that can later be replaced by
% something more logical

properties
  Models  % the array of TModel objects 
  Grid    % TGrid3D object  or T3DfwdFDE
          %  T3DfwdFDE contains sigmaCell proerties and a few functions to
          %  work with it
end

methods
%*******************************************************************
function obj = TGridInterpolate(Models, Grid);
  %   class constructor ... simple
  obj.Models = Models;
  obj.Grid = Grid;
end
%*******************************************************************
%*******************************************************************
function Interpolate(obj)
  
  

end


end     % methods
end    % classdef