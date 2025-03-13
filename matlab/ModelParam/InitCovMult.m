function [CmHalf] = InitCovMult(grid)
%  Usage: [CmHalf] = InitCovMult(grid);
%  Returns parameters for Cm^{1/2};
rho = .2;
zMax = 10000;
nOuter = 4;
nYinner = 4;
nZinner = 2;
CmHalf = struct('rho',rho,'zMax',zMax,'nOuter',nOuter,...
        'nYinner',nYinner,'nZinner',nZinner,'grid',grid);
