function [S,header] = JmultT_MTX(sigma0,d)
%  Calculates product of transpose of sensitivity matrix with
%  data vector ... this version does not sum output conductivity
%   parameter over frequencies (and returns real/imag parts!)
%
%  Usage: [S,header] = JmultT_MTX(sigma0,d);
%  
%  Inputs:	sigma0 = conductivity parameter (structure)
%		d = data vector (cell array of structures)
%  Output:	S = cell array of conductivcity structures
%		header = header written to file by fortran program
%  Need to create scratch directory, write grid file (Input.grd)
%   before calling this function

%   write out conductivity parameter
cfile = 'scratch/Input1.cpr';
writeCond2D(cfile,sigma0);
%   write out data vector
cfile = 'scratch/Input.imp';
writeZ(cfile,d);
%  execute program with -M option to multiply by transpose of 
%        sensitivity matrix without summing over frequencies
[S] = unix('CondDataMaps -M')
%  read in resulting sensitivities
cfile = 'scratch/Out.sns';
%  NOTE: as currently configured only real parts of J^T d are
%    computed by CondDataMaps
[S,header] = readCondMTX(cfile);
