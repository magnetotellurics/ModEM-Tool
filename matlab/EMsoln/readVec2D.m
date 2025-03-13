function [vec] = readVec2D(cfile)
% Usage:  [vec] = readVec2D(cfile)
%  Reads in vec2D object, returns as structure
fid = fopen(cfile,'r');
ll = fread(fid,1,'long');
s = fread(fid,80,'char');
s = deblank(char(s'));
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
nm = fread(fid,2,'long');
ny = nm(1);
nz = nm(2);
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
temp = fread(fid,[2*ny,nz],'float64');
v = temp(1:2:end,:)+ i *temp(2:2:end,:); 
vec = struct('v',v,'gridType',s);
fclose(fid)
