function [gridDef] = readGrid2D(cfile)
% Usage:  [gridDef] = readGrid2D(cfile)
fid = fopen(cfile,'r');
ll = fread(fid,1,'long');
Ny = fread(fid,1,'long');
Nz = fread(fid,1,'long');
Nza = fread(fid,1,'long');
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
Dy = fread(fid,Ny,'float64');
ll = fread(fid,1,'long');
ll = fread(fid,1,'long');
Dz = fread(fid,Nz,'float64');
ll = fread(fid,1,'long');
fclose(fid);
gridDef = struct('Ny',Ny,'Nz',Nz,'Nza',Nza,'Dy',Dy','Dz',Dz');
