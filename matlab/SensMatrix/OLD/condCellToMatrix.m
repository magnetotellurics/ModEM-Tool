function [M] = condCellToMatrix(J);
%   convert cell array of conductivity parameters to a standard
%    matlab matrix: column k of the output matrix corresponds
%    to the conductivity parameter stored in cell k
%  Usage: [M] = condCellToMatrix(J);

nSens = prod(size(J));
nParam = prod(size(J{1}.v));
M = zeros(nSens,nParam);
for k  = 1:nSens
   M(k,:) = reshape(J{k}.v,1,nParam);
end
