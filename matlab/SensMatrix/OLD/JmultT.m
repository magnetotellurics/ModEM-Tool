function [m] = JmultT(sigma0,d)
%  Calculates product of transpose of sensitivity matrix with
%  data vector, returning a single model paramter
%
%  Usage: [m] = JmultT(sigma0,d);
%  
%  Inputs:	sigma0 = conductivity parameter (structure)
%		d = data vector (cell array of structures)
%  Output:	m = conductivcity parameter structure
%  Need to create scratch directory, write grid file (Input.grd)
%   before calling this function

%   write out conductivity parameter
cfile = 'scratch/Input1.cpr';
writeCond2D(cfile,sigma0);
%   write out data vector (template)
cfile = 'scratch/Input.imp';
writeZ(cfile,d);
%  execute program with -T option to multiply by transpose of 
%        sensitivity matrix without summing over frequencies
s = unix('CondDataMaps -T')
%  read in resulting sensitivities
cfile = 'scratch/Out.cpr';
[m] = readCond2D(cfile);
