nx = 20;
ny = 20;

dx = 1 * ones(nx, 1);
dy = 1 * ones(ny, 1);
dz = xygrid.logz(0, 100, 36, 1.2);
grid = xygrid(dx, dy, dz);

npad = 4;
grid = grid.pad('NSEW', npad);

mstruct = defaultm('tranmerc')
mstruct.origin = [0.0, 0.0, 0.0];
%mstruct.scalefactor = 0.9996;
mstruct.scalefactor = 1.0;
mstruct.geoid = almanac('earth', 'wgs84', 'meters');
mstruct = defaultm(mstruct);

data_type = "Full_Impedance"; is_complex = 0;
%data_type = "Off_Diagonal_Impedance"; is_complex = 0;
%data_type = "Full_Vertical_Components"; is_complex = 0;
%data_type = "Full_Interstation_TF"; is_complex = 0;
%data_type = "Off_Diagonal_Rho_Phase"; is_complex = 1;
%data_type = "Phase_Tensor"; is_complex = 1;

disp("Testing: " + data_type)

per = 10.^(0.6:0.6:4);
dat = mtdata(data_type, per, 'xy');
dat = dat.gridNodes(grid, per, 'all', mstruct);
dat = zero(dat);
dat.header = 'Miles Test';

expected_nsites = nx * ny + (nx+ny) + 1;

assert(dat.nSites == expected_nsites, "FAILED!");

status = dat.write('test.dat', 'list', is_complex);
%writeApres_3d('test.apres3d.dat', dat.d, "header")

mt_read = mtdata.read("test.dat", "list", dat.units, data_type);
