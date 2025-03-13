function [J,header] = makeSens(sigma0,d)
%  Makes sensitivity matrix for 2D MT TE problem
%
%  Usage: [J,header] = makeSens(sigma0,d);
%  
%  Inputs:	sigma0 = conductivity parameter (structure)
%		d = data vector (cell array of structures)
%  Output:	J = sensitivity matrix (cell array of
%		      conductivity parameter structures)
%  Need to create scratch directory, write grid file (Input.grd)
%   before calling this function

%   write out conductivity parameter
cfile = 'scratch/Input1.cpr';
writeCond2D(cfile,sigma0);
%   write out data vector (template)
cfile = 'scratch/Input.imp';
writeZ(cfile,d);
%  execute program with -S option to compute sensitivity matrix
[s] = unix('CondDataMaps -S');
%  read in sensitivity matrix
cfile = 'scratch/Out.sns';
[J,header] = readCondMTX(cfile);
