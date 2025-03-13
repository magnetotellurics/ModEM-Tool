function [Jxdm] = Jmult(dm,m0,d0)
%  Calcualtes product of sensitivity matrix and conductivity
%   parameter perturbation delSigma, by calling fortran program
%    CondDataMaps;  This emulates subroutine Jmult in module SensMatrix
%
%  Usage: [Jxdm] = Jmult(dm,m0,d0);
%  
%  Inputs:      dm = conductivity parameter (structure)
%         	    m0 = background conductivity parameter (structure)
%		        d0 = data vector (cell array of structures; template)
%
%  Output:	Jxdm = data vector (cell array of structures)
%
%  Need to create scratch directory, write grid file (Input.grd)
%   before calling this function
% [MULT_BY_J]
%  -M  rFile_Model rFile_dModel rFile_Data wFile_Data
%   Multiplies a model by J to create a data vector

if exist('scratch','dir')==0
    mkdir('scratch');
end
%   write scratch file for input conductivity parameter
m0_File = 'scratch/m0.cpr';
writeCond_2D(m0_File,m0);
%   write scratch file for input data vector (template)
d0_File = 'scratch/d0.imp';
writeZ_2D(d0_File,d0);
%   write scratch file for input conductivity parameter perturbation
dm_File = 'scratch/dm.cpr';
writeCond_2D(dm_File,dm);
%   file for output data, to be read after running program
Jxdm_File = 'scratch/Jxdm.imp';

%  execute program with -G option to multiply by sensitivity matrix
Test2D('MULT_BY_J',m0_File,dm_File,d0_File,Jxdm_File)
[Jxdm] = readZ_2D(Jxdm_File);

%  read in resulting data vector
[Jxdm] = readZ_2D(Jxdm_File);
