function [dOut] = fwdPred(sigma0,dIn)
%  Computes predicted data for conductivity sigma0
%
%  Usage: [dOut] = fwdPred(sigma0,dIn);
%  
%  Inputs:	sigma0 = conductivity parameter (structure)
%		dIn = data vector (cell array of structures)
%  Output:	dOut = predicted data vector (cell array of structures)
%  Need to create scratch directory, write grid file (scratch/Input.grd)
%   before calling this function

%   write out conductivity parameter
cfile = 'scratch/Input1.cpr';
writeCond2D(cfile,sigma0);
%   write out data vector (template)
cfile = 'scratch/Input.imp';
writeZ(cfile,dIn);
%  execute program with -S option to compute sensitivity matrix
[s] = unix('CondDataMaps -F');
%  read in predicted data
cfile = 'scratch/Out.imp';
[dOut] = readZ(cfile);
