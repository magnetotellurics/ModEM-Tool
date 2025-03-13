function [dOut] = Jmult_MTX(delSigma,sigma,dIn)
%  Usage : [dOut] = Jmult_MTX(delSigma,sigma,dIn);
% multiply cell array of model paramters (one for
%   each frequency) by sensitivity matrix, using call to 
%   fortran program
%  Inputs:      delSigma = conductivity parameters (cell array
%                          of structures)
%               sigma0 = background conductivity parameter (structure)
%               dIn = data vector (cell array of structures; template)
%               eAll is optional (not yet implemented)
%  Output:      dOut = data vector (cell array of structures)


% write out background conductivity parameter
cfile = 'scratch/Input1.cpr';
writeCond_2D(cfile,sigma);
%   write out data vector (template)
cfile = 'scratch/Input.imp';
writeZ_2D(cfile,dIn);
%  compute Jm for each 
cfile = 'scratch/Input2.sns';
header = '';
status = writeCondMTX_2D(cfile,delSigma,header);
%  execute program with -N option to multiply by
%        sensitivity matrix frequency by frequency
[S] = unix('CondDataMaps -N')
%  read in resulting sensitivities
cfile = 'scratch/Out.sns';
%  NOTE: as currently configured only real parts of J^T d are
%    computed by CondDataMaps
[S,header] = readCondMTX_2D(cfile);
%  read in resulting data vector
cfile = 'scratch/Out.imp';
[dOut] = readZ_2D(cfile);
