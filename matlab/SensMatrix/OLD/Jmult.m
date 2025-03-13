function [dOut] = Jmult(delSig,sigma0,dIn,eAll)
%  Calcualtes product of sensitivity matrix and conductivity
%   parameter perturbation delSigma, by calling fortran program
%    CondDataMaps;  This emulates subroutine Jmult in module SensMatrix
%
%  Usage: [dOut] = Jmult(delSig,sigma0,dIn,eAll);
%  
%  Inputs:      delSigma = conductivity parameter (structure)
%         	sigma0 = background conductivity parameter (structure)
%		dIn = data vector (cell array of structures; template)
%               eAll is optional (not yet implemented)
%  Output:	dOut = data vector (cell array of structures)
%  Need to create scratch directory, write grid file (Input.grd)
%   before calling this function

%   write out conductivity parameter
cfile = 'scratch/Input1.cpr';
writeCond2D(cfile,sigma0);
%   write out data vector (template)
cfile = 'scratch/Input.imp';
writeZ(cfile,dIn);
%   write out conductivity parameter perturbation
cfile = 'scratch/Input2.cpr';
writeCond2D(cfile,delSig);
%  execute program with -G option to multiply by sensitivity matrix
[s] = unix('CondDataMaps -G');
%  read in resulting data vector
cfile = 'scratch/Out.imp';
[dOut] = readZ(cfile);
