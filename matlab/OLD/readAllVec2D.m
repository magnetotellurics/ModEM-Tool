function [VEC] = readAllVec2D(cfile)
% Usage:  [VEC] = readAllVec2D(cfile)
%  Reads in a sequence of vec2D objects,
%    returns as a cell array of structure
fid = fopen(cfile,'r');
ll = fread(fid,1,'long');
nVec = fread(fid,1,'long');
ll = fread(fid,1,'long');

for k = 1:nVec
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
   ll = fread(fid,1,'long');
   v = temp(1:2:end,:)+ i *temp(2:2:end,:); 
   VEC{k} = struct('v',v,'gridType',s);
end
fclose(fid);
