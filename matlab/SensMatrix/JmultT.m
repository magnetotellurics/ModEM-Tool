function [JTxd] = JmultT(m0,d)
%  Calculates product of transpose of sensitivity matrix with
%  data vector, returning a single model paramter
%
%  Usage: [m] = JmultT(m0,d);
%  
%  Inputs:	sigma0 = conductivity parameter (structure)
%		d = data vector (cell array of structures)
%  Output:	m = conductivcity parameter structure
%  Need to create scratch directory, write grid file (Input.grd)
%   before calling this function

%   write out conductivity parameter
m0_File = 'scratch/m0.cpr';
writeCond_2D(m0_File,m0);
%   write out data vector (template)
d_File = 'scratch/d.imp';
writeZ_2D(d_File,d);
%  file name for output sensitivity (model parameter)
JTxd_File = 'scratch/JTxd.cpr';

Test2D('MULT_BY_J_T',m0_File,d_File,JTxd_File)
[JTxd] = readCond_2D(JTxd_File);