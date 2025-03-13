%%
clear

PLOT_FITS = true;

%COV = 'laplacian';
COV = 'standard';
%  Script for testing object oriented inversion

%%%%%%   Edit these lines
%inverseType = 'DCG';
%inverseType = 'DASOCC';
inverseType = 'OCCAM';
%inverseType = 'BIDIAG';
%inverseType = 'BIDIAGmtx';
%inverseType = 'BDORTH';
%inverseType = 'BDMTX';
%modelID = 'Test7';
%dataID = '_05';
modelID = 'Test0';
dataID = '_05';
ModelDir = '.';
run = 'v';
MODE = 'JT';
ErrFac = 1.0;
priorCond = .003;
nu = 1;
%   file names
DataFile = [modelID dataID '_' MODE '.imp'];
condFile = [modelID '.cpr'];
%   conductivity section plot parameters
nZplot = 32;
nYskip = 12;
cax = [-3.5,-.5];
%%%%%%   Edit above lines

OPTIONS = struct('nYskip',nYskip,'nZplot',nZplot,...
		'cax',cax,'title','True Model');
%   model and data vector types ...
m_True = MT2DmodelParam();
d = MT2DZ();
  
modelID = [modelID run];
Root = [ MODE '_' modelID dataID '_' inverseType];
plotFile =  [Root '.eps'];
if strcmp(MODE,'JT')
    fitFileTE = ['fitTE_' modelID dataID '_' inverseType '.eps'];
    fitFileTM = ['fitTM_' modelID dataID '_' inverseType '.eps'];
else
    fitFile = ['fit' MODE '_' modelID dataID '_' inverseType '.eps'];
end

%   load model file (m_True ... a template is needed)
[m_True] = readVec(m_True,condFile);
%   plot input true conductivity
[h_True,hCB] =  plotCond(m_True,OPTIONS);

%  load data vector
d = readVec(d,DataFile);

%  prior conductivity: uniform half space
[m_0] = InitHalfSpace(m_True,priorCond);
%   set starting model ... same as prior often
m=m_0;
%m_0 = oneD_average(m_True);

if strcmp(COV,'laplacian')
    %   for Laplacian smoothing use zero prior
    m_0 = 0*m_0;
    %   beta controls relative vertical/horizontal smoothing
    beta  = 1;
    Cm = MT2D_laplace(extractGrid(m_True),beta);
else
    %  set smoothing parameters:obj.ITER.RMStolerance = 1.0;
    Cm = MT2DmodelCov(extractGrid(m_True));
end

J = MT2Dsens(d);
obj = GN(d,J,nu,m_0,Cm,m);
obj.ITER.iterMax = 10;
%% 
switch inverseType
    case 'DCG'
        %    Test of CG solver in class GN
        CGsoln = GN(d,J,nu,m_0,Cm);
        CGsoln.ITER.iterMax = 10;
        IterCG=InvIterControl;
        IterCG.iterMax = 30;
        DCG(CGsoln,IterCG)
        eval(['save ' Root '.mat CGsoln']);
        OPTIONS.title = ['DCG, \nu = ' num2str(obj.nu)];
        plotCond(CGsoln.m,OPTIONS)
        if PLOT_FITS  
            obj =  CGsoln;
        end
    case 'OCCAM'
        %
        %obj.ITER.RMStolerance = 1.0;
        obj.ITER.PhaseI_Max = 12;
        obj.ITER.PhaseII_Max = 1;
        obj.Occam;
        OPTIONS.title = ['Occam, \nu = ' num2str(obj.nu)];
        plotCond(obj.m,OPTIONS)
        
    case 'DASOCC'
        %    Test of HybOcc solver in class GN: case "DASOCC"
        %obj = GN(d,J,nu,m_0,Cm);
        obj.ITER.PhaseI_Max = 10;
        obj.ITER.PhaseII_Max = 1;
        obj.ITER.ModifiedOccam = false;
        Algorithm = 'DASOCC';
        
        HybOcc(obj,Algorithm)
        DASOCCsoln = obj;
        eval(['save ' Root '.mat DASOCCsoln']);
        OPTIONS.title = ['Modified Occam, \nu = ' num2str(obj.nu)];
        plotCond(obj.m,OPTIONS)
    case 'BIDIAG'
        %   Test of HybOcc solver in class GN: case "BIDIAG"
        %obj = GN(d,J,nu,m_0,Cm);
        Iter=InvIterControl;
        Iter.iterMax = 30;
        obj.ITER.PhaseI_Max = 10;
        obj.ITER.ModifiedOccam = false;
        Algorithm = 'BIDIAG';
        HybOcc(obj,Algorithm,Iter)
        BIDIAGsoln = obj;
        eval(['save ' Root '.mat BIDIAGsoln']);
        OPTIONS.title = ['Occam, \nu = ' num2str(obj.nu)];
        plotCond(obj.m,OPTIONS)
    case 'BIDIAGmtx'
        %   Test of HybOcc solver in class GN: case "BIDIAG"
        %obj = GN(d,J,nu,m_0,Cm);
        Iter=InvIterControl;
        Iter.iterMax = 30;
        obj.ITER.PhaseI_Max = 10;
        obj.ITER.ModifiedOccam = false;
        Algorithm = 'BIDIAGmtx';
        HybOcc(obj,Algorithm,Iter)
        BIDIAGsoln = obj;
        eval(['save ' Root '.mat BIDIAGmtx_soln']);
        OPTIONS.title = ['Occam, \nu = ' num2str(obj.nu)];
        plotCond(obj.m,OPTIONS)
    case 'BDORTH'
        %   Test of HybOcc solver in class GN: case "BDORTH"
        %obj = GN(d,J,nu,m_0,Cm);
        Iter=InvIterControl;
        Iter.iterMax = 30;
        obj.ITER.PhaseI_Max = 10;
        Algorithm = 'BDORTH';
        HybOcc(obj,Algorithm,Iter)
        BDORTHsoln = obj;
        eval(['save ' Root '.mat BDORTHsoln']);
        OPTIONS.title = ['Occam, \nu = ' num2str(obj.nu)];
        plotCond(obj.m,OPTIONS)
    case 'BDMTX'
        %    Test of HybOcc solver in class GN: case "BDMTX"
        %obj = GN(d,J,nu,m_0,Cm);
        %obj.ITER.RMStolerance = 1.02;
        Algorithm = 'BDMTX';
        %ITER = cell(6,1);
        %for k = 1:6
        %    ITER{k}=InvIterControl;
        %end
        %ITER{1}.ResidualTolerance = .1;
        %ITER{2}.ResidualTolerance = .1   ;
        %ITER{3}.ResidualTolerance = .04;
        %ITER{4}.ResidualTolerance = .02;    h = plotPseud(obj.pred.d,OPTIONS);

        %ITER{5}.ResidualTolerance = .01;
        %ITER{6}.iterMax = 30;
        %ITER{6}.ResidualTolerance = .005;
        %Iter.iterMax = 30;
        Iter=InvIterControl;
        Iter.iterMax = 30;
        obj.ITER.PhaseI_Max = 10;
        obj.ITER.ModifiedOccam = false;
        HybOcc(obj,Algorithm,Iter)
        BDMTXsoln = obj;
        eval(['save ' Root '.mat BDMTXsoln']);
        OPTIONS.title = ['Hybrid MTX \nu = ' num2str(obj.nu)];
        plotCond(obj.m,OPTIONS)
end
%%

if PLOT_FITS   
    OPTIONS.rho_cax = [0 3.5];
    OPTIONS.phi_cax = [0 90];
    h = plotPseud(obj.pred.d,OPTIONS);
    if length(h) == 2
        figure(h(1))
        eval(['print -depsc ' fitFileTE]);
        figure(h(2))
        eval(['print -depsc ' fitFileTM]);
    else
        figure(h)
        eval(['print -depsc ' fitFile]);
    end
end