% function [] = genslice(MODFILE,direction)
% genslice - function to generate slices from a 3D model file.
%
%  NOTES: 
%  
%
%  See also XXX.
%  --------------------
%  Bo Yang, 2014.
%  Institute of Geophysics and Geomatics,
%  China University of Geosciences, Wuhan, China.
%  Comments, bug reports and questions, please send to:
%  yangbo.cug@163.com.
%  Copyright 2014-2018, Bo Yang, IGG, CUG.
%  $Revision: 1.0 $ $Date: 2014/02/26 $
%  Last changed: 2014/02/26 14:31:49.

%  Revision log:
%  2014/02/26 : Version 1.0 released.
MODFILE = 'Block2-iter41.rho';
direction = 'Y';
%
%  read in the model file.
%
cond = readCond_3D(MODFILE,2);
nx = cond.grid.Nx;
ny = cond.grid.Ny;
nz = cond.grid.NzEarth;
dx = cond.grid.dx;
dy = cond.grid.dy;
dz = cond.grid.dz;
x0 = cond.grid.origin(1);
y0 = cond.grid.origin(2);
z0 = cond.grid.origin(3);
xn = x0 + cumsum([0;dx]);
yn = y0 + cumsum([0;dy]);
zn = z0 + cumsum([0;dz]);
sige = zeros(nx+1,ny+1,nz+1);
sige(1:nx,1:ny,1:nz) = cond.v;
v  = sige;
%
%  generate the slices one by one.
%
switch lower(direction)
case {'x','yz'}
	for k = 1:nx
		sig2d = squeeze(v(k,:,:)).';
		pcolor(yn,zn,sig2d);
		xlabel('y');
		ylabel('z');
		title(['kx = ',num2str(k),'  X = ',num2str(xn(k))]);
		colorbar;
		caxis([-7 0]);
		set(gca,'YDir','reverse');
		% add other reference info.
		% save gcf
		FIGNAME = ['figs/','X',sprintf('%3.3d',k)];
		saveas(gcf,FIGNAME,'fig');
		close(gcf);
	end
case {'y','xz'}
	for k = 1:ny
		sig2d = squeeze(v(:,k,:)).';
		pcolor(xn,zn,sig2d);
		xlabel('x');
		ylabel('z');
		title(['ky = ',num2str(k),'  Y = ',num2str(yn(k))]);
		colorbar;
		caxis([-7 0]);
		set(gca,'YDir','reverse');
		% add other reference info.
		% save gcf
		FIGNAME = ['figs/','Y',sprintf('%3.3d',k)];
		saveas(gcf,FIGNAME,'fig');
		close(gcf);
	end
case {'z','xy'}
	for k = 1:nz
		sig2d = squeeze(v(:,:,k));
		pcolor(yn,xn,sig2d);
		xlabel('y');
		ylabel('x');
		title(['kz = ',num2str(k),'  Z = ',num2str(zn(k))]);
		colorbar;
		caxis([-7 0]);
		% add other reference info.
		% save gcf
		FIGNAME = ['figs/','Z',sprintf('%3.3d',k)];
		saveas(gcf,FIGNAME,'fig');
		close(gcf);
	end
otherwise
	error('Unsupported direction!');
end
%
% End of the function.
%
