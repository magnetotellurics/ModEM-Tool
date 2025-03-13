

%%
DataFile = 'dTest30.imp';
condFile = 'mTrue.cpr';
priorCond = .02;

m_True = MT3DmodelParam();
d = MT3DZ();

%   load model file (m_True ... a template is needed)
[m_True] = readVec(m_True,condFile);
%   plot input true conductivity
plotCond(m_True);

%  load data vector
d = readVec(d,DataFile);

%  prior conductivity: uniform half space
[m_0] = InitHalfSpace(m_True,priorCond);

%  set smoothing parameters:obj.ITER.RMStolerance = 1.0;
Cm = MT3DmodelCov(extractGrid(m_True));
%Cm = Mod3DMTCov();

%mSm = CovMult(Cm,m_True);
%plotCond(mSm);

%   copmpute prediction from "true model"
pred = fwd(m_True,d);

%   plot data and truth ...
plot(ExtractVec(pred))
hold on;
plot(ExtractVec(d),'r--')

%   see if rms is reasonable (note: dot by default uses data errors in d)
rms = sqrt(dot(d-pred,d-pred)/lengthDat(d))

%   create indirect sensitivity
%J = SensitivityIndirect(d,m_0);

%   make a model parameter with random perturbations
%dm = m_True;
%dm.v = .1*randn(size(dm.v));
%   multiply by J
%Jxdm = J_times_m(J,dm);

%  multiply J^T d
%JTd = JT_times_d(J,d);

%    make full sensitivity matrix ...
%Jfull = MT3Dsens(d,m_0);
J = MT3Dsens(d);

%   now multiply using full sensitivity
%Jfull_xdm = J_times_m(Jfull,dm);
%Jfull_Td = JT_times_d(Jfull,d);

%   test rowJ
%J_k = rowJ(J,20);
%Jfull_k = rowJ(Jfull,20);

nu = 1;
obj = GN(d,J,nu,m_0,Cm);
obj.ITER.PhaseI_Max=3;
obj.ITER.PhaseII_Max = 0;
obj.ITER. RMStolerance=1.1;
%obj.ITER.modifiedOccam=false;
%   do inversion using method HybOcc (this works for DASOCC, BIDIAG,
%   BDORTH, BDMTX; not DCG;  result is in obj.m when this completes.
inverseType = 'DASOCC';
obj.HybOcc(inverseType);
%obj.DCG


