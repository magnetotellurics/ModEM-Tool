%%
%  Script for testing object oriented inversion, 3D MT
for iRun = 1:length(RUNS)
    
    %%%%%%   Edit these lines to change input field, inversion type, plotting
    %%%%%%   options, etc.
    COV = RUNS{iRun}.Cov;
    inverseType = RUNS{iRun}.inverseType;
    DataFile = RUNS{iRun}.DataFile;
    condFile = RUNS{iRun}.condFile;
    priorCond = RUNS{iRun}.priorCond;
    modelID = RUNS{iRun}.modelID;
    dataID = RUNS{iRun}.dataID;
    ModelDir = RUNS{iRun}.ModelDir;
    run = RUNS{iRun}.run;
    ErrFac = RUNS{iRun}.ErrFac;
    nu = RUNS{iRun}.nu;
    SENS = RUNS{iRun}.SENS;
    IterMax = RUNS{iRun}.IterMax;
    PhaseI_Max = RUNS{iRun}.PhaseI_Max;
    PhaseII_Max = RUNS{iRun}.PhaseII_Max;
    ModifiedOccam = RUNS{iRun}.ModifiedOccam;
    DCGmax = RUNS{iRun}.DCGmax;
    DCGtol = RUNS{iRun}.DCGtol;
    
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
    switch COV
        case 'MT3DmodelCov'
            Cm = MT3DmodelCov(m_True.grid);
        case 'MT3DmodelCovWt'
            Cm = MT3DmodelCovWt(m_True.grid);
    end
    
    %   use either the full Jacobian ....
    switch SENS
        case 'FULL'
            J = MT3Dsens(d);
        case 'Indirect'
            J = SensitivityIndirect(d,m_0);
    end
    obj = GN(d,J,nu,m_0,Cm,m);
    obj.ITER.iterMax = IterMax;
    obj.ITER.PhaseI_Max = PhaseI_Max;
    obj.ITER.PhaseII_Max = PhaseII_Max;
    obj.ITER.ModifiedOccam = ModifiedOccam;
    %
    tic
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
            Algorithm = 'DASOCC';
            HybOcc(obj,Algorithm)
                time = toc;
                nTIMES = length(obj.TIMES);
                obj.TIMES{nTIMES+1} = ...
                    struct('step','total','T',time,'iter',obj.ITER.niter);
            DASOCCsoln = obj;
            eval(['save ' Root '.mat DASOCCsoln']);
            OPTIONS.title = ['Modified Occam, \nu = ' num2str(obj.nu)];
            plotCond(obj.m)
        case 'BIDIAG'
            %   Test of HybOcc solver in class GN: case "BIDIAG"
            %obj = GN(d,J,nu,m_0,Cm);
            Iter=InvIterControl;
            Iter.iterMax = DCGmax;
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
            Iter.ResidualTolerance = DCGtol;
            Iter.iterMax = DCGmX;
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
            Iter.ResidualTolerance = DCGtol;
            Iter.iterMax = DCGmax;
            Algorithm = 'BDORTH';
            HybOcc(obj,Algorithm,Iter)
                time = toc;
                nTIMES = length(obj.TIMES);
                obj.TIMES{nTIMES+1} = ...
                    struct('step','total','T',time,'iter',obj.ITER.niter);
            BDORTHsoln = obj;
            eval(['save ' Root '.mat BDORTHsoln']);
            OPTIONS.title = ['Occam, \nu = ' num2str(obj.nu)];
            plotCond(obj.m)
        case 'BDMTX'
            %    Test of HybOcc solver in class GN: case "BDMTX"
            Algorithm = 'BDMTX';
            Iter=InvIterControl;
            Iter.ResidualTolerance = DCGtol;
            Iter.iterMax = DCGmax;
            HybOcc(obj,Algorithm,Iter)
                time = toc;
                nTIMES = length(obj.TIMES);
                obj.TIMES{nTIMES+1} = ...
                    struct('step','total','T',time,'iter',obj.ITER.niter);
            BDMTXsoln = obj;
            eval(['save ' Root '.mat BDMTXsoln']);
            OPTIONS.title = ['Hybrid MTX \nu = ' num2str(obj.nu)];
            plotCond(obj.m)
    end
                   
                   
end