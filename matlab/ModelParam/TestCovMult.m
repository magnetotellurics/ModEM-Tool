%  script used for initial testing of CovMult
%
gridFile = 'Test1.grd';
grid = readGrid2D(gridFile);
condFile = 'Test1.cpr';
[cond] = readCond2D(condFile);
[NyCond,NzCond] = size(cond.v);
grid.Nza = grid.Nz - NzCond;
z = [0 cumsum(grid.Dz(grid.Nza+1:end))];
y = [0 cumsum(grid.Dy)];
sig = [cond.v cond.v(:,end)];
sig = [sig; sig(end,:)];
pcolor(y,z,log10(sig)')
%  smoothing parameters:
rho = .2;
zMax = 10000;
nOuter = 3;
nYinner = 3;
nZinner = 2;
smthParams = struct('rho',rho,'zMax',zMax,'nOuter',nOuter,...
        'nYinner',nYinner,'nZinner',nZinner);
