%  Script sensTest is used to test sensitivity calculation in fortran 90

path(path,'/home/server/pi/homes/egbert/EM3D/Modular2D/matlab/SensMatrix')
path(path,'/home/server/pi/homes/egbert/EM3D/Modular2D/matlab/DataSpace')
path(path,'/home/server/pi/homes/egbert/EM3D/Modular2D/matlab/EMsoln')
path(path,'/home/server/pi/homes/egbert/EM3D/Modular2D/matlab/ModelParam')

gridFile = 'Test1.grd';

PLOTfields = 0;
epsilon = 0.01;
LOGCOND = 1;
logtest = 0;
MODE = 'TE';

if(LOGCOND)
  condFile = 'TestLog.cpr';
  dCondFile = 'dTestLog.cpr';
  TEZFile = ['TestPredLog' MODE '.imp'];
  dTEZFile = ['TestDeltaPredLog' MODE '.imp'];
  SensMultFile = ['TestSensMultLog' MODE '.imp'];
  SensTransMultFile = ['TestSensMultTransLog' MODE '.cpr'];
  SensMatFile = ['TestLog' MODE '.sns'];
else
  condFile = 'Test1.cpr';
  dCondFile = 'dTest1.cpr';
  TEZFile = 'TestPred.imp';
  dTEZFile = 'TestDeltaPred.imp';
  SensMultFile = 'TestSensMult.imp';
  SensMatFile = 'Test1.sns';
end

%  load in grid
[gridDef] = readGrid2D(gridFile);
Ny = gridDef.Ny;
Nz = gridDef.Nz;
Nza = gridDef.Nza;

%  conductivity
[S] = readCond2D(condFile);

%  perturbation to conductivity
[dS] = readCond2D(dCondFile);

% measured impedances for sigma0
[Z0] = readZ(TEZFile);

% measured impedances for sigma0+dsigma
[dZ] = readZ(dTEZFile);

% impedance linearizatino
% [dZimp] = readZ(LinPredFile);

% Compute perturbation of solution ...
%    (divide by epsilon to compare to JxDs)
nPer = length(Z0);
for iPer = 1:nPer
   dZ{iPer}.Z =  (dZ{iPer}.Z-Z0{iPer}.Z)/epsilon;
end

% forward sensitivity calculation
[JxdS] = readZ(SensMultFile);

% backward sensitivity calculation
[JTxd] = readCond2D(SensTransMultFile);

% load  sensitvity matrix
[J,header] = readCondMTX(SensMatFile);

% compute JxdS from sensitivity matrix ...
JxdSComp = JxdS;
ii = 0;
for iPer = 1:nPer
  nSites = length(JxdS{iPer}.Z);
  for iSite = 1:nSites
     ii = ii + 1;
     Zr =  sum(sum(J{ii}.v .* dS.v ));
     ii = ii + 1;
     Zi =  sum(sum(J{ii}.v .* dS.v ));
     JxdSComp{iPer}.Z(iSite) = Zr+i*Zi;
  end
end

%   compute JxdS from sensitivity matrix, using
%   --> sensitivity matrix for linear conductivity
%   --> perturbation dS of log conductivity
%   --> background conductivity ... 
%       should match JxdS from log case (to within rounding)

if logtest & (1-LOGCOND)
%  perturbation to log conductivity
[dS1] = readCond2D('dTestLog.cpr');
temp = dS1.v.*S.v;
JxdSCompLog = JxdS;
ii = 0;
for iPer = 1:nPer
  nSites = length(JxdS{iPer}.Z);
  for iSite = 1:nSites
     ii = ii + 1;
     Zr =  sum(sum(J{ii}.v .* temp ));
     ii = ii + 1;
     Zi =  sum(sum(J{ii}.v .* temp ));
     JxdSCompLog{iPer}.Z(iSite) = Zr+i*Zi;
  end
end
end 

%  plot comparison
plotBands = [1 4 7 10];
for k = 1:length(plotBands)
   T = Z0{plotBands(k)}.T;
   figure('Position',[100,100,400,700])
   axes('Position',[.1,.55,.8,.35])
   plot(real(JxdS{plotBands(k)}.Z),'ks')
   hold on
   plot(real(JxdSComp{plotBands(k)}.Z),'bo')
   plot(real(dZ{plotBands(k)}.Z),'r+')
   title(['Period = ', num2str(T)])
   axes('Position',[.1,.1,.8,.35])
   plot(imag(JxdS{plotBands(k)}.Z),'ks')
   hold on
   plot(imag(JxdSComp{plotBands(k)}.Z),'bo')
   plot(imag(dZ{plotBands(k)}.Z),'r+')
end

% compute JTxd from sensitivity matrix ...
JTxdComp = JTxd;
JTxdComp.v = 0;
ii = 0;
for iPer = 1:nPer
  nSites = length(Z0{iPer}.Z);
  for iSite = 1:nSites
     Z = Z0{iPer}.Z(iSite);
     temp = J{ii+1}.v - i*J{ii+2}.v;
     JTxdComp.v = JTxdComp.v + real(Z*temp);
     ii = ii + 2;
  end
end
