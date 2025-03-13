function [Dy,Dz,Cond,periods,nza] = rdGridCond(cfile)
%Usage:   [Dy,Dz,Cond,periods,nza] = rdGridCond(cfile);
fid = fopen(cfile,'r');
ll = fread(fid,1,'long');
Iword = fread(fid,4,'long');
nz = Iword(1);
nza = Iword(2);
ny = Iword(3);
nPer = Iword(4);
ll = fread(fid,2,'long');
Dy = fread(fid,ny,'float64');
ll = fread(fid,2,'long');
Dz = fread(fid,nz,'float64');
ll = fread(fid,2,'long');
Cond = fread(fid,[ny,nz],'float64');
ll = fread(fid,2,'long');
periods = fread(fid,nPer,'float64');
fclose(fid);
