%%
%  Script for testing object oriented inversion, 3D MT
clear

%%%%%%   Edit these lines to change input field, inversion type, plotting
%%%%%%   options, etc.
PLOT_FITS = true;
%COV = 'laplacian';
COV = 'standard';
%inverseType = 'DCG';
inverseType = 'DASOCC';
%inverseType = 'OCCAM';
%inverseType = 'BIDIAG';
%inverseType = 'BIDIAGmtx';
%inverseType = 'BDORTH';
%inverseType = 'BDMTX';

DataFile = 'dTest3_30_05.imp';
condFile = 'mNew_3.cpr';
priorCond = .01;
modelID = 'New_3_';
dataID = 'dTest3_30_05';
ModelDir = '.';
run = 'tst';

ErrFac = 1.0;
nu = 1;

%%%%%%   Edit above lines
modelID = [modelID run];
Root = [modelID dataID '_' inverseType];

%   create model and data vector types ...
m_True = MT3DmodelParam();
d = MT3DZ();
  
%   load model file (m_True ... a template is needed)
[m_True] = readVec(m_True,condFile);

%   plot input true conductivity
plotCond(m_True);

% load data vector
d = readVec(d,DataFile);
setenv('NMPIPROC',num2str(2*d.NTX+1))
sites = d.d{1}.siteLoc;
x = sites(:,1);
y = sites(:,2);
hold on
plot(y,x,'ko','linewidth',2)

%  prior conductivity: uniform half space
[m_0] = InitHalfSpace(m_True,priorCond);

%   set starting model ... same as prior
m=m_0;
for iTx = 1:d.NTX
    d.d{iTx}.Zerr = d.d{iTx}.Zerr*ErrFac;
end

%   set covariance
%Cm = Mod3DMTCov()
Cm = MT3DmodelCov(m_True.grid);
CmWt = MT3DmodelCovWt(m_True.grid);

%   use either the full Jacobian ....
J = MT3Dsens(d);
%    or indirect sensitivity calculations
%J = SensitivityIndirect(d,m_0);
obj = GN(d,J,nu,m_0,Cm,m);
obj.ITER.iterMax = 3;
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
        plotCond(CGsoln.m)
        if PLOT_FITS  
            obj =  CGsoln;
        end
    case 'OCCAM'
        %  only run this with full sensitivity computation
        %obj.J = MT3Dsens(d);
        %obj.ITER.RMStolerance = 1.0;
        obj.ITER.PhaseI_Max = 2 ;
        obj.ITER.PhaseII_Max = 1;
        obj.Occam;
        OCCAMsoln = obj;
        %    clear senstivity before saving
        OCCAMsoln.J = [];
        OCCAMsoln.R.J = [];
        OCCAMsoln.R.R = [];
        OCCAMsoln.R.U = [];
        eval(['save ' Root '.mat OCCAMsoln']);
        OPTIONS.title = ['Occam, \nu = ' num2str(obj.nu)];
        plotCond(obj.m)
        
    case 'DASOCC'
        %    Test of HybOcc solver in class GN: case "DASOCC"
        %obj = GN(d,J,nu,m_0,Cm);
        obj.ITER.PhaseI_Max = 2;
        obj.ITER.PhaseII_Max = 1;
        obj.ITER.ModifiedOccam = false;
        Algorithm = 'DASOCC';
        
        HybOcc(obj,Algorithm)
        DASOCCsoln = obj;
        eval(['save ' Root '.mat DASOCCsoln']);
        OPTIONS.title = ['Modified Occam, \nu = ' num2str(obj.nu)];
        plotCond(obj.m)
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
        plotCond(obj.m)
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
        plotCond(obj.m)
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
        plotCond(obj.m)
    case 'BDMTX'
        %    Test of HybOcc solver in class GN: case "BDMTX"
        Algorithm = 'BDMTX';
        Iter=InvIterControl;
        Iter.iterMax = 30;
        obj.ITER.PhaseI_Max = 10;
        obj.ITER.ModifiedOccam = false;
        HybOcc(obj,Algorithm,Iter)
        BDMTXsoln = obj;
        eval(['save ' Root '.mat BDMTXsoln']);
        OPTIONS.title = ['Hybrid MTX \nu = ' num2str(obj.nu)];
        plotCond(obj.m)
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