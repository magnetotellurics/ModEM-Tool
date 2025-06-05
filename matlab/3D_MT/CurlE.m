function [H] = CurlE(E,Grid);
% Usage:  [H] = CurlE(E,Grid);

Nx = length(Grid.dx);
Ny = length(Grid.dy);
Nz = length(Grid.dz);
%  Hx
diff1 = zeros(Nx+1,Ny,Nz)+i*zeros(Nx+1,Ny,Nz);
diff2 = diff1;
for iy = 1:Ny
   diff1(:,iy,:) = (E.z(:,iy+1,:)-E.z(:,iy,:))/Grid.dy(iy);
end
for iz = 1:Nz
   diff2(:,:,iz) = (E.y(:,:,iz+1)-E.y(:,:,iz))/Grid.dz(iz);
end
H.x = diff1-diff2;
%  Hy
diff1 = zeros(Nx,Ny+1,Nz)+i*zeros(Nx,Ny+1,Nz);
diff2 = diff1;
for iz = 1:Nz
   diff1(:,:,iz) = (E.x(:,:,iz+1)-E.x(:,:,iz))/Grid.dz(iz);
end
for ix = 1:Nx
   diff2(ix,:,:) = (E.z(ix+1,:,:)-E.z(ix,:,:))/Grid.dx(ix);
end
H.y = diff1-diff2;
%  Hz
diff1 = zeros(Nx,Ny,Nz+1)+i*zeros(Nx,Ny,Nz+1);
diff2 = diff1;
for ix = 1:Nx
   diff1(ix,:,:) = (E.y(ix+1,:,:)-E.y(ix,:,:))/Grid.dx(ix);
end
for iy = 1:Ny
   diff2(:,iy,:) = (E.x(:,iy+1,:)-E.x(:,iy,:))/Grid.dy(iy);
end
H.z = diff1-diff2;
