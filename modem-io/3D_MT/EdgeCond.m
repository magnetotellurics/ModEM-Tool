function S = EdgeCond(grid,rho);
%   Usage:  S = EdgeCond(grid,rho);
AirCond = 1e-10;
[Nx,Ny,NzE] = size(rho);
NzA = grid.NzAir;
Nz = NzE+NzA;
t = zeros(Ny+1,NzE+1);
s = t; 
A = grid.dy*grid.dz(NzA+1:end)';
S.x = ones([Nx,Ny+1,Nz+1])*AirCond;
for ix = 1:Nx
  temp = A./squeeze(rho(ix,:,:));
  s(1:Ny,1:NzE) = temp;
  t(1:Ny,1:NzE) = A;
  s(2:Ny+1,1:NzE) = s(2:Ny+1,1:NzE)+temp;
  t(2:Ny+1,1:NzE) = t(2:Ny+1,1:NzE)+A;
  s(1:Ny,2:NzE+1) = s(1:Ny,2:NzE+1)+temp;
  t(1:Ny,2:NzE+1) = t(1:Ny,2:NzE+1)+A;
  s(2:Ny+1,2:NzE+1) = s(2:Ny+1,2:NzE+1)+temp;
  t(2:Ny+1,2:NzE+1) = t(2:Ny+1,2:NzE+1)+A;
  S.x(ix,:,NzA+1:Nz+1) = s./t;
end

t = zeros(Nx+1,NzE+1);
s = t+i*t;
A = grid.dx*grid.dz(NzA+1:end)';
S.y = ones([Nx+1,Ny,Nz+1])*AirCond;
for iy = 1:Ny
  temp = A./squeeze(rho(:,iy,:));
  s(1:Nx,1:NzE) = temp;
  t(1:Nx,1:NzE) = A;
  s(2:Nx+1,1:NzE) = s(2:Nx+1,1:NzE)+temp;
  t(2:Nx+1,1:NzE) = t(2:Nx+1,1:NzE)+A;
  s(1:Nx,2:NzE+1) = s(1:Nx,2:NzE+1)+temp;
  t(1:Nx,2:NzE+1) = t(1:Nx,2:NzE+1)+A;
  s(2:Nx+1,2:NzE+1) = s(2:Nx+1,2:NzE+1)+temp;
  t(2:Nx+1,2:NzE+1) = t(2:Nx+1,2:NzE+1)+A;
  S.y(:,iy,NzA+1:Nz+1) = s./t;
end

t = zeros(Nx+1,Ny+1);
s = t+i*t;
A = grid.dx*grid.dy';
S.z = ones([Nx+1,Ny+1,Nz])*AirCond;
for iz = 1:NzE
  temp = A./squeeze(rho(:,:,iz));
  s(1:Nx,1:Ny) = temp;
  t(1:Nx,1:Ny) = A;
  s(2:Nx+1,1:Ny) = s(2:Nx+1,1:Ny)+temp;
  t(2:Nx+1,1:Ny) = t(2:Nx+1,1:Ny)+A;
  s(1:Nx,2:Ny+1) = s(1:Nx,2:Ny+1)+temp;
  t(1:Nx,2:Ny+1) = t(1:Nx,2:Ny+1)+A;
  s(2:Nx+1,2:Ny+1) = s(2:Nx+1,2:Ny+1)+temp;
  t(2:Nx+1,2:Ny+1) = t(2:Nx+1,2:Ny+1)+A;
  S.z(:,:,NzA+iz) = s./t;
end
S.CondRead = 1;
