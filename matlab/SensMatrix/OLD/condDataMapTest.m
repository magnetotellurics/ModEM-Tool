%   script CondDataMapTest : used to test fortran90 program CondDataMap

path(path,'/home/server/pi/homes/egbert/EM3D/Modular2D/matlab/SensMatrix')
path(path,'/home/server/pi/homes/egbert/EM3D/Modular2D/matlab/DataSpace')
path(path,'/home/server/pi/homes/egbert/EM3D/Modular2D/matlab/EMsoln')
path(path,'/home/server/pi/homes/egbert/EM3D/Modular2D/matlab/ModelParam')

global J

gridFile = 'Test1.grd';
MODE = 'TE';
LOGCOND = 1;

if(LOGCOND)
  m0_File = 'TestLog.cpr';
  m_File = 'dTestLog.cpr';
  d_File = ['TestPredLog' MODE '.imp'];
  J_File = ['TestLog' MODE '.sns'];
else
  m0_File = 'Test1.cpr';
  m_File = 'dTest1.cpr';
  d_File = ['TestPred.imp'];
  J_File = ['Test1.sns'];
end

%  load in grid
[gridDef] = readGrid2D(gridFile);
Ny = gridDef.Ny;
Nz = gridDef.Nz;
Nza = gridDef.Nza;

% background model parameter
[m0] = readCond2D(m0_File);

%  test model parameter
[m1] = readCond2D(m_File);

% test data vector
[d1] = readZ(d_File);

% sensitvity matrix
[J,header] = readCondMTX(J_File);

%   test of fwdPred
[d2] = fwdPred(m0,d1);

%   test of makeSens matrix
[J1,header] = makeSens(m0,d1);

%  test of internal/external multiplication of m1 by J
[d_int] = J_times_m(m1,d1);
[d_ext] = Jmult(m1,m0,d1);

%  test of internal/external multiplication of d1 by J
[m_int] = JT_times_d(d1,m1)
[m_ext] = JmultT(m0,d1)

%  plot comparison
plotBands = [1 4 7 10];
for k = 1:length(plotBands)
   T = d_int{plotBands(k)}.T;
   figure('Position',[100,100,400,700])
   axes('Position',[.1,.55,.8,.35])
   plot(real(d_int{plotBands(k)}.Z),'ks')
   hold on
   plot(real(d_ext{plotBands(k)}.Z),'bo')
   title(['Period = ', num2str(T)])
   axes('Position',[.1,.1,.8,.35])
   plot(imag(d_int{plotBands(k)}.Z),'ks')
   hold on
   plot(imag(d_ext{plotBands(k)}.Z),'bo')
end

pclr(m_int.v); title('Matlab computation from J')
pclr(m_ext.v); title('Direct CondDataMaps computation')

%  test of internal/external multiplication of d1 by J: MTX case
[M_int] = JT_times_d_MTX(d1,m1);
[M_ext] = JmultT_MTX(m0,d1);
for k = 1:length(plotBands)
   pclr(M_int{1,k}.v); title('Matlab computation from J')
   pclr(M_ext{k}.v); title('Direct CondDataMaps computation')
end

%  test of internal/external multiplication of m1 by J: MTX case
[d_int] = J_times_m_MTX(M_ext,d1);
[d_ext] = Jmult_MTX(M_ext,m0,d1);

%  plot comparison
plotBands = [1 4 7 10];
for k = 1:length(plotBands)
   T = d_int{plotBands(k)}.T;
   figure('Position',[100,100,400,700])
   axes('Position',[.1,.55,.8,.35])
   plot(real(d_int{plotBands(k)}.Z),'ks')
   hold on
   plot(real(d_ext{plotBands(k)}.Z),'bo')
   title(['Period = ', num2str(T)])
   axes('Position',[.1,.1,.8,.35])
   plot(imag(d_int{plotBands(k)}.Z),'ks')
   hold on
   plot(imag(d_ext{plotBands(k)}.Z),'bo')
end
