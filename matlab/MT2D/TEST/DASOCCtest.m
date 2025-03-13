%    Simple driver for 2D inversion, synthetic test data;  See
%    invDriver2D.m in INVclass for a more complete/complex example
%     NEED TO MAKE SURE THERE IS A LINK TO compiled Mod2DMT in the directory you run
%     this from!
inverseType = 'DASOCC';   %   other choices: DCG, BIDIAG, BDORTH, BDMTX

%   file names
%    data file, 40 sites, TE + TM, 5% error added
DataFile = 'Test0_05_JT.imp';
%   true conductivity ... use as template for grid, model parameter object
TrueModelFile = 'Test0.cpr';
priorCond = .01;
%   OPTIONS is used to pass some parameters to conductivity section ploting
%    function p
OPTIONS = struct('nYskip',12,'nZplot',32,...
		'cax',[-3.5,-.5],'title','True Model');
    
%  create model and data vector objects ...
m_True = MT2DmodelParam();
d = MT2DZ();

%   load model file (m_True ... a template is needed)
[m_True] = readVec(m_True,TrueModelFile);
%   plot input true conductivity
[h_True,hCB] =  plotCond(m_True,OPTIONS);

%  load data vector
d = readVec(d,DataFile);

%  prior conductivity: uniform half space
priorCond = .01 ;  % conductivity in S/m
[m_0] = InitHalfSpace(m_True,priorCond);

%  set model covariance
Cm = MT2DmodelCov(extractGrid(m_True));

%   create sensitivity object ...
%     In this case the object stores the full sensitivity matrix ...
%  ... but called this way the actual sensitivity computation is not yet done
Jfull = MT2Dsens(d);

%    create Gauss-Newton inversion object ... inputs include data,
%    sensitivity object (not yet computed), tradeoff paraameter (not used
%    in many cases) prior, covariance (optional
%    starting model perturbation)
nu = 1;
obj = GN(d,Jfull,nu,m_0,Cm);
%   do inversion using method HybOcc (this works for DASOCC, BIDIAG,
%   BDORTH, BDMTX; not DCG;  result is in obj.m when this completes.
obj.HybOcc(inverseType);

%  plot result
OPTIONS.title = ['Modified Occam, \nu = ' num2str(obj.nu)];
plotCond(obj.m,OPTIONS)

%   save result to mat file
DASOCCsoln = obj;
save DASOCCtestResults.mat DASOCCsoln