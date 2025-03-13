function [grid,rho] = rdModelRM(cfile)
% usage :  [grid,rho] = rdModelRM(cfile);
% reads 3D model input file in R. Mackie format
%  HAS NOT BEEN DEBUGGED FOR RESISTIVITY CODES!

fid = fopen(cfile,'r');
cline = fgets(fid);
I2 = length(cline);
[Nx,COUNT,ERRMSG,NEXTINDEX]  = sscanf(cline,'%d',1);
I1 = NEXTINDEX;
[Ny,COUNT,ERRMSG,NEXTINDEX]  = sscanf(cline(I1:I2),'%d',1);
I1 = I1 + NEXTINDEX-1;
[Nz,COUNT,ERRMSG,NEXTINDEX]  = sscanf(cline(I1:I2),'%d',1);
I1 = I1 + NEXTINDEX-1;
[NzAir,COUNT,ERRMSG,NEXTINDEX]  = sscanf(cline(I1:I2),'%d',1);
I1 = I1 + NEXTINDEX-1;
VALUE = -1;
for k = I1:I2
   if(upper(cline(k)) == 'V')
       VALUE = 1;
       break
   elseif(upper(cline(k)) == 'M')
       VALUE = 0;
       break
   end
end
if(VALUE == -1)
  fprintf(1,'%s\n','Error in rdModelRM: unknown option for VALUE/CODE')
  return
end
dx = fscanf(fid,'%f',Nx);
dy = fscanf(fid,'%f',Ny);
dz = fscanf(fid,'%f',Nz);
Cond = zeros(Nx,Ny,Nz); 
if(VALUE)
   for k = 1:Nz
      iz = fscanf(fid,'%d',1);
      OneLayer = fscanf(fid,'%f',[Nx,Ny]);
      rho(:,:,iz) = OneLayer;
   end
else
   rhoInd = rho;
   for k = 1:Nz
      iz = fscanf(fid,'%d',1);
      OneLayer = fscanf(fid,'%d',[Nx,Ny]);
      rhoInd(:,:,iz) = OneLayer;
   end
%   FOLLOWING IS WRONG: Where do we get # of codes?
   [rhoMap, nrho] = fscanf(fid,'%f',inf);
   for k = 1:nrho
      rho(rhoInd == k-1) = rhoMap(k);
   end
end
%  read rest of file (what is this???)
%   add to grid structure as appropriate
cline = fgets(fid);
cline = fgets(fid);
cline = fgets(fid);
cline = fgets(fid);
[Origin,count] = fscanf(fid,'%f',3);
if count < 3
   nfill = 3-count;
   Origin = [Origin ; zeros(nfill,1)];
end
[theta,count] = fscanf(fid,'%f',1);
if count == 0
   theta = 0.0
end
grid = struct('dx',dx,'dy',dy,'dz',dz,'NzAir',NzAir,...
    'origin',Origin,'rotation',theta);
