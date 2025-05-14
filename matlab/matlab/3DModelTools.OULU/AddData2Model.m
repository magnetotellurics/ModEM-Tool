function [SMap,VarMap,grid] = AddData2Model(grid,prior,data,k)
% Interpolation of conductances with Gaussian distributions
% based on a prior and data
%
% INPUTS
% grid      - 3d grid we want to use
% prior     - 3d prior model object
% data      - conductance data object (1d, 2d or 3d)
% k         - smoothing factor
%
% OUTPUTS
% SMap      - estimated conductances with the given grid
% VarMap    - estimated variances of conductances with the given grid 
% grid      - the given grid

