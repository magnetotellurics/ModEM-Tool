dx = 10 * ones(40, 1);
dy = 10 * ones(60, 1);
dz = xygrid.logz(0, 100, 36, 1.2);

grid = xygrid(dx, dy, dz);

npad = 4;
grid.pad('NSEW', npad);


ind = grid.mask('cylinder', [5 15], [0 100], 100, 5);

sigma = xymodel(grid);
sigma.v = -2 * ones(grid.nx, grid.ny, grid.nzEarth);

sigma = sigma.log10();

sigma.v(ind == 5) = 1;

sigma.uiplot();

sigma.plot('depth',10); delta=1e-6; caxis([-4-delta 0+delta]);
print('-djpeg','-r300',['example_depth_' num2str(10) 'km.jpg']);

sigma.write('example.model.rho', 'WS');

%
% Creating Projection and Data file
%
utmstruct = defaultm('tranmerc');
utmstruct.origin = [45 -90 0];
utmstruct.scalefactor = 0.9996;
utmstruct.geoid = almanac('earth', 'wgs84', 'meters');
utmstruct.falseeasting=500000;
mstruct = defaultm(utmstruct);

per = 10.^(0.6:0.6:4);
full_imp = mtdata('Full_Impedance', per, 'xy');
full_vert = mtdata('Full_Vertical_Components', per, 'xy');

full_imp = full_imp.gridNodes(grid, per, 'all', mstruct);
full_vert = full_vert.gridNodes(grid, per, 'all', mstruct);

full_imp = zero(full_imp);
full_vert = zero(full_vert);

full_imp.header = 'Full_Impedance, written out by ModEM-Examples Tutorial';
full_vert.header = 'Full_Vertical_Components, written out by ModEM-Examples Tutorial';

full_imp.write('data.full_impedance.6per.dat');
full_vert.write('data.full_vert_comp.6per.dat');