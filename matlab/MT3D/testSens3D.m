%%
DataFile = 'dTest30.imp';
condFile = 'mTrue.cpr';
priorCond = .01;

%   create 3D model parameter, data vectors  <<<< change from 2D <<<<<<<<<<
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

%  set smoothing parameters   <<<< change from 2D <<<<<<<<<<<<<<<<<<<<<<<<<
Cm = MT3DmodelCov(extractGrid(m_True));

%   copmpute prediction from "true model"
pred = fwd(m_True,d);

%   plot data and truth ...
plot(ExtractVec(pred))
hold on;
plot(ExtractVec(d),'r--')

%   see if rms is reasonable (note: dot by default uses data errors in d)
rms = sqrt(dot(d-pred,d-pred)/lengthDat(d))

%   create sensitivity matrix object ...
%   full sensitivity object ...
Jind = SensitivityIndirect(d,m_0);
%     ... and indirect object <<<< change from 2D <<<<<<<<<<<<<<<<<<<<<<<<<
Jfull = MT3Dsens(d,m_0);   
save Jfull.mat Jfull

%   make a model parameter with random perturbations
dm = m_True;
dm.v = .1*randn(size(dm.v));
%   multiply by J
Jfull_dm = J_times_m(Jfull,dm);
Jind_dm = J_times_m(Jind,dm);
%   plot results for comparison -- should be identical (to within roundoff)
figure
plot(ExtractVec(Jfull_dm),'b-');
hold on
plot(ExtractVec(Jind_dm),'r--');

%  multiply J^T d
Jfull_Td = JT_times_d(Jfull,d);
Jind_Td = JT_times_d(Jind,d);
%   plot results for comparison -- should be identical (to within roundoff)
%   <<<<<<<<<<<<<<<< change from 2D <<<<<<<<<<<<<<<<<<<<<<<<<
plotCond(Jfull_Td)
plotCond(Jind_Td)

%%   now do the same computations with normalized J
JtildeInd =  copy(Jind);
JtildeFull = copy(Jfull);

Normalize(JtildeFull,Cm);
Normalize(JtildeInd,Cm);

%   multiply by J
JtildeFull_dm = J_times_m(JtildeFull,dm);
JtildeInd_dm = J_times_m(JtildeInd,dm);
%   plot results for comparison -- should be identical (to within roundoff)
figure
plot(ExtractVec(JtildeFull_dm),'b-');
hold on
plot(ExtractVec(JtildeInd_dm),'r--');

%  multiply J^T d
d = Normalize(d);
Jfull_Td = JT_times_d(JtildeFull,d);
Jind_Td = JT_times_d(JtildeInd,d);
%   plot results for comparison -- should be identical (to within roundoff)
plotCond(Jfull_Td)
plotCond(Jind_Td)


%%   Now consider representer matrix ...  multiplication of a data vector by
%   any of these should have the same effect!
nu = 1;
%   representer matrix built on Jfull
R_full = RepMat(Jfull,Cm,nu);
RtildeFull = RepMat(JtildeFull,Cm,nu);
R_ind = RepMat(Jind,Cm,nu);
RtildeInd = RepMat(JtildeInd,Cm,nu);

Rd_1 = R_full*d;
Rd_2 = RtildeFull*d;
Rd_3 = R_ind*d;
Rd_4 = RtildeInd*d;
plot(ExtractVec(Rd_1),'b-')
hold on
plot(ExtractVec(Rd_2),'r--')
figure
plot(ExtractVec(Rd_1),'b-')
hold on
plot(ExtractVec(Rd_3),'r--')
figure
plot(ExtractVec(Rd_1),'b-')
hold on
plot(ExtractVec(Rd_4),'r--')

%   now form R for full cases   ... Don't do this for 3D!!!!
%R_full.formR
%RtildeFull.formR
%Rd_5 = R_full*d;
%Rd_6 = RtildeFull*d;
%figure
%plot(ExtractVec(Rd_1),'b-')
%hold on
%plot(ExtractVec(Rd_5),'r--')
%figure
%plot(ExtractVec(Rd_1),'b-')
%hold on
%plot(ExtractVec(Rd_6),'r--')
