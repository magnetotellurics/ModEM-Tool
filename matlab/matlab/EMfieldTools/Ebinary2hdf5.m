% ESOLN: binary to hdf5; hdf5 to binary

nPer = 1;
%fname_hdf5 = 'OUTPUT/2023-08-17-hdf5-SP2-serial/CylinderModel_1per.h5';
fname_hdf5 = 'CylinderModel_1per.h5';
fname_binary = 'OUTPUT/2023-08-17-stable-prerelease-serial/CylinderModel_1per.soln';

% /Grid/ attribute nGrids=1,units='m',type='xy',method='finite difference'
% /Grid/Subgrid01/ attributes Nx,Ny,Nz,Nza,origin,rotation=0
% /Grid/Subgrid01/dx
% /Grid/Subgrid01/dy
% /Grid/Subgrid01/dz
% /Efield/ attributes: txType='MT' or e.g., txType='MT,CSEM' (comma-separated)
% /Efield/MT/ attributes: nPeriods, nModes
% /Efield/MT/Periods
% /Efield/MT/Modes
% /Efield/MT/Tx01/Mode01/ attribute Period=value, Mode='X'
% /Efield/MT/Tx01/Mode01/Subgrid01/x attribute dim=(Nx,Ny+1,Nz+1)
% /Efield/MT/Tx01/Mode01/Subgrid01/y attribute dim=(Nx+1,Ny,Nz+1)
% /Efield/MT/Tx01/Mode01/Subgrid01/z attribute dim=(Nx+1,Ny+1,Nz)
% /Efield/MT/Tx01/Mode02/ attribute Period=value, Mode='Y'
% etc
E = readE(fname_binary,nPer); % plotEMsoln(E);

% For now we're assuming only one grid, but allow for multigrid
nGrids = 1; igrid = 1;

% /Grid
datagroup = sprintf('/Grid/Subgrid%02d/',igrid);
h5create(fname_hdf5,[datagroup 'dx'],[E.grid.Nx 1])
h5write(fname_hdf5,[datagroup 'dx'],E.grid.Dx)
h5create(fname_hdf5,[datagroup 'dy'],[E.grid.Ny 1])
h5write(fname_hdf5,[datagroup 'dy'],E.grid.Dy)
h5create(fname_hdf5,[datagroup 'dz'],[E.grid.Nz 1])
h5write(fname_hdf5,[datagroup 'dz'],E.grid.Dz)
h5writeatt(fname_hdf5,datagroup,'Nx', string(E.grid.Nx));
h5writeatt(fname_hdf5,datagroup,'Ny', string(E.grid.Ny));
h5writeatt(fname_hdf5,datagroup,'Nz', string(E.grid.Nz));
h5writeatt(fname_hdf5,datagroup,'Nza', string(E.grid.Nza));
h5writeatt(fname_hdf5,datagroup,'origin', string(E.grid.origin));
h5writeatt(fname_hdf5,datagroup,'rotation', string(E.grid.rotation));
h5writeatt(fname_hdf5,'/Grid/','nGrids', string(nGrids));
h5writeatt(fname_hdf5,'/Grid/','units', E.grid.units);
h5writeatt(fname_hdf5,'/Grid/','type', E.grid.type);
h5writeatt(fname_hdf5,'/Grid/','method', 'finite difference');

% /Efield
for j=1:length(E.Periods)
    for i=1:length(E.Modes)
        Efield = E.E{i,j};
        datagroup = sprintf('/Efield/MT/Tx%02d/Mode%02d/Subgrid%02d/',j,i,igrid);
        h5create(fname_hdf5,[datagroup 'xreal'],Efield.NdX,'Datatype','double');
        h5write(fname_hdf5,[datagroup 'xreal'],real(Efield.x));
        h5create(fname_hdf5,[datagroup 'ximag'],Efield.NdX,'Datatype','double');
        h5write(fname_hdf5,[datagroup 'ximag'],imag(Efield.x));
        h5create(fname_hdf5,[datagroup 'yreal'],Efield.NdY,'Datatype','double');
        h5write(fname_hdf5,[datagroup 'yreal'],real(Efield.y));
        h5create(fname_hdf5,[datagroup 'yimag'],Efield.NdY,'Datatype','double');
        h5write(fname_hdf5,[datagroup 'yimag'],imag(Efield.y));
        h5create(fname_hdf5,[datagroup 'zreal'],Efield.NdZ,'Datatype','double');
        h5write(fname_hdf5,[datagroup 'zreal'],real(Efield.z));
        h5create(fname_hdf5,[datagroup 'zimag'],Efield.NdZ,'Datatype','double');
        h5write(fname_hdf5,[datagroup 'zimag'],imag(Efield.z));
        h5writeatt(fname_hdf5,datagroup,'Period', string(E.Periods(j)));
        h5writeatt(fname_hdf5,datagroup,'Mode', string(E.Modes{i}));
    end
end
h5writeatt(fname_hdf5,'/Efield/','txType', 'MT');
h5writeatt(fname_hdf5,'/Efield/MT/','nPeriods', string(length(E.Periods)));
h5writeatt(fname_hdf5,'/Efield/MT/','nModes', string(length(E.Modes)));
h5create(fname_hdf5,'/Efield/MT/Periods',size(E.Periods))
h5write(fname_hdf5,'/Efield/MT/Periods',E.Periods)
h5create(fname_hdf5,'/Efield/MT/Modes',size(E.Modes),'Datatype','string')
h5write(fname_hdf5,'/Efield/MT/Modes',string(E.Modes))

% /
h5writeatt(fname_hdf5,'/','creation_date', date);
h5writeatt(fname_hdf5,'/','original_file', fname_binary);

h5disp(fname_hdf5)