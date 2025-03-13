function [status] = writeVec2D(cfile,vec)
% Usage:  [status] = writeVec2D(cfile,vec)
%  writes vec2D object, provided as structure
%   now status is total number of bytes written
fid = fopen(cfile,'w');
ll = 80;
status = 0;
count = fwrite(fid,ll,'long');
status = status+4*count;
s = blanks(80);
s(1:length(vec.gridType)) = vec.gridType;
count = fwrite(fid,s,'char');
status = status+count;
count = fwrite(fid,ll,'long');
status = status+4*count;
Nyz = size(vec.v);
ll = 8;
count = fwrite(fid,ll,'long');
status = status+4*count;
count = fwrite(fid,Nyz,'long');
status = status+4*count;
count = fwrite(fid,ll,'long');
status = status+4*count;
ll = Nyz(1)*Nyz(2)*16;
count = fwrite(fid,ll,'long');
status = status+4*count;
temp = zeros(2*Nyz(1),Nyz(2));
temp(1:2:end,:) = real(vec.v);
temp(2:2:end,:) = imag(vec.v);
count = fwrite(fid,temp,'float64');
status = status+8*count;
count = fwrite(fid,ll,'long');
status = status+4*count;
fclose(fid);
