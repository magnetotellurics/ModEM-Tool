function [J,header] = makeSens(sigma0,d)
%  Makes sensitivity matrix for 2D MT TE problem
%   Writes sigma0, d to files to be read by Fortran program Mod2DMT
%   calls program with option to  compute sensitivity matrix, reads result
%    and returns result in J 
%
%  Usage: [J,header] = makeSens(sigma0,d);
%  
%  Inputs:	sigma0 = conductivity parameter (structure)
%		d = data vector (cell array of structures)
%  Output:	J = sensitivity matrix (cell array of
%		      conductivity parameter structures)
%             This is modified to output three cells:
%              first contains the data vector (so that we can keep
%                track of what sensitivity goes with what!)
%             second contains the actual model parameters which
%               define the sensitivity (one per observation)
%             third contains the grid

%    make scratch directory if it doesn't exist ...
if exist('scratch','dir')==0
    mkdir('scratch');
end

%   write out conductivity parameter in scratch directory
m_File = 'scratch/Input.cpr';
writeCond_2D(m_File,sigma0);

%   write out data vector template in scratch directory
d0_File = 'scratch/Input.imp';
writeZ_2D(d0_File,d);
%   file name for computed sensitivity matrix, to be output by Mod2DMT
J_File = 'scratch/Out.sns';

%   run Mod2DMT
Test2D('COMPUTE_J',m_File,d0_File,J_File)

%  read in sensitivity matrix
J = cell(2,1);
J{1} = d;
[J{2},header] = readCondMTX_2D(J_File);
J{3} = sigma0.grid; 
