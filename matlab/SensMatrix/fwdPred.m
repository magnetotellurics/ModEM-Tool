function [d] = fwdPred(m,d0)
%  Computes predicted data for conductivity sigma0
%   Writes sigma0, dIn to files to be read by Fortran program Mod2DMT
%   calls program with option to  compute predicted data, reads result
%    and returns result in dOut 
%
%  Usage: [d] = fwdPred(m,d0);
%  
%  Inputs:	m = conductivity parameter (structure)
%		d0 = data vector (cell array of structures)
%  Output:	d = predicted data vector (cell array of structures)

%    make scratch directory if it doesn't exist ...
if exist('scratch','dir')==0
    mkdir('scratch');
end

%   write out conductivity parameter in scratch directory
m_File = 'scratch/Input.cpr';
writeCond_2D(m_File,m);

%   write out data vector template in scratch directory
d0_File = 'scratch/Input.imp';
writeZ_2D(d0_File,d0);
%   file name for computed data vector, to be output by Mod2DMT
d_File = 'scratch/Out.imp';

%   run Mod2DMT
Test2D('FORWARD',m_File,d0_File,d_File)

%  read in predicted data
d = readZ_2D(d_File);
