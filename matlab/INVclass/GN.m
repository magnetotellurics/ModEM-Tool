
classdef GN < handle
    % class for for GN-type inversion algorithms
    %
    %   Gary D. Egbert, 2010
    %   College of Oceanic and Atmospheric Sciences
    %   Oregon State University
    
    properties
        J       %   sensitivity matrix handle
        Cm      %   model covariance handle (nominally sqrt)
        R       %   representer matrix object handle
        m_0     %   prior model parameter (value object)
        m       %   current model parameter (value object)
        m_n     %   model parameter deviation from prior: m_n = m-m_0
        mTilde  %   m_n = Cm*mTilde
        d       %   data (data space value object)
        b       %   representer coefficients (data space value object)
        mPerpTilde  %   component of current normalized model not in span
                %       representers
        pred    %   f(m)  (data space value object)
        res     %   residual d - f(m)  (data space value object)
        dHat    %   d-f(m)+Jm  (data space value object)
        %nIter   %   iteration number
        iterInner    %  array of InvIterControl object for inner loop
        ITER    %   InvIterControl object for outer loop
        nu      %   damping parameter
        mu      %   parameter associated with mPerpTilde
        TIMES    %   for timing variaous steps
        
    end
    
    methods
        %******************************************************************
        function obj = GN(d,J,nu,m_0,Cm,m)
        % class constructor
            obj.d = d;
            obj.J = J;
            obj.nu = nu;
            obj.m_0 = m_0;
            obj. Cm = Cm;
            obj.ITER = InvIterControl();
            obj.iterInner = InvIterControl();
            if nargin == 5
                obj.m = m_0;
                obj.m_n = m_0;
                obj.m_n = ZeroVec(obj.m_n);
            else
                obj.m = m;
                obj.m_n = obj.m-m_0;
            end
        end
        %******************************************************************
        function RMSmisfit  = ForwardModel(obj)
        %   using current m_n update pred, res, dHat
            obj.pred = fwd(obj.m,obj.d);
            obj.res = obj.d-obj.pred;
            RMSmisfit = sqrt(real(dot(obj.res,obj.res))/length(obj.res));
        end
        %******************************************************************
        function SetRepMat(obj,J)
        %   create/update representer matrix object (new sensitivity object)
            if obj.ITER.niter == 1 ;
                %   if first iteration, create RepMat object
                obj.R = RepMat(J,obj.Cm,obj.nu);
            else
                %   otherwise just update RepMat
                updateRepMat(obj.R,J);
            end
        end
        %******************************************************************
        function UpdateModel(obj)
        %   once representer coefficient is determined, update model
        %   parameter, compute new predicted, residual, etc.
        %   Need to use J^T (or multiplication by J^T) ... several variants
        %   of this are supported, and we need to keep straight if J is
        %   already "normalized"
        
            if obj.J.normalized
                if ~obj.b.normalized
                    error('In GN.UpdateModel: J is normalized but coefficients b are not')
                else
              %      nCov = 1;
                end
            else
                %   if the coefficients are for the "normalized" problem
                %     they need to be *multiplied* (yes, that's right) by error
                %     standard deviations  ... in this case where J is not
                %     normalized
                if obj.b.normalized
                    obj.b = MultCdInv(obj.b);
                    obj.b.normalized = false;
                end
               % nCov = 2;
            end
        
            %  compute linear combination of representers
            obj.mTilde = JT_times_d(obj.J,obj.b);
            if ~obj.J.normalized
                obj.mTilde = CovMult(obj.Cm,obj.mTilde,1);
            end
            %   smooth result with covariance: only once for J normalized;
            %   otherwise twice (i.e., with sqrt(Cm) in first instance, Cm
            %   otherwise)
            obj.m_n = CovMult(obj.Cm,obj.mTilde,1);
            %  add result to prior
            obj.m = obj.m_0+obj.m_n;
        end
        %******************************************************************
        function DCG(obj,IterCG)
        %    DCG solver
        %    Usage: DCG(obj,IterControl)
        
            if nargin ==  1
                %   use default iter control objectd
                IterCG=InvIterControl;
            end
            obj.ITER.niter =  1;
            obj.ITER.RMS(1) = ForwardModel(obj);
            display(['Initial Misfit = ' num2str(obj.ITER.RMS(1))])
            DONE = FitToTolerance(obj.ITER);
            while  ~DONE
                %  itercontrol object for CG (inner-loop) iterations
                obj.iterInner(obj.ITER.niter) = IterCG;
                %   update sensitivity matrix object (new model parameter)
                updateSens(obj.J,obj.m);
                %  update dHat
                obj.dHat =  obj.res+J_times_m(obj.J,obj.m_n);
                %   create/update representer matrix object (new sensitivity object)
                SetRepMat(obj,obj.J);
                %   normalize dHat
                dTilde = Normalize(obj.dHat,1);
                %   Use CG to solve for representer coefficients
                %      (everything normalized by data error sd)
                [obj.b,obj.iterInner(obj.ITER.niter)] = ...
                    CG(obj.R,dTilde,obj.iterInner(obj.ITER.niter));   
                %  update model parameter
                obj.UpdateModel;
                %   prepare for next iteration ...
                obj.ITER.niter = obj.ITER.niter+1;       
                obj.ITER.RMS(obj.ITER.niter) = ForwardModel(obj);
                n = obj.ITER.niter;
                display(['Iteration # ' num2str(n)])
                display(['Misfit = ' num2str(obj.ITER.RMS())])
                %   ... then check for convergence 
                DONE = FitToTolerance(obj.ITER);
            end
        end
        %******************************************************************
        function Occam(obj) 
        %   Occam solver    
        %   Usage: after creation of GM object obj, obj.Occam
        
            obj.ITER.niter =  1;
            obj.ITER.iterMax = obj.ITER.PhaseI_Max;
            obj.ITER.RMS(1) = ForwardModel(obj);
            display(['Initial Misfit = ' num2str(obj.ITER.RMS(1))])
            %   Occam Phase I
            PHASE1_DONE = FitToTolerance(obj.ITER);
            while ~PHASE1_DONE
                %  itercontrol object for CG (inner-loop) iterations
                obj.iterInner(obj.ITER.niter) = InvIterControl;
                %   update sensitivity matrix object (new model parameter)
                updateSens(obj.J,obj.m);
                %  update dHat
                obj.dHat =  obj.res+J_times_m(obj.J,obj.m_n);
                %   create/update representer matrix object (new sensitivity object)
                SetRepMat(obj,obj.J);
                %   form and factor cross-product matrix for normal
                %   equations
                obj.R.formR;
                obj.R.factorR;
                %   select trade-off parameter
                selectTradeOff(obj);
                obj.iterInner(obj.ITER.niter).OccamPhase = 1;
                %  update model parameter
                obj.UpdateModel;
                %  increment iteration count
                obj.ITER.niter = obj.ITER.niter+1; 
                %  run forward model, compute misfit ....
                %          ... this could be made more efficient
                obj.ITER.RMS(obj.ITER.niter) = ForwardModel(obj);
                % check fit to tolerance
                PHASE1_DONE = FitToTolerance(obj.ITER);
                %   some display ...
                n = obj.ITER.niter-1;
                display(['Phase I : Iteration # ' num2str(n)])
                display(['Misfit = ' num2str(obj.ITER.RMS(n+1))])
            end
            display('Done with Phase I')
            %    Occam Phase II
            for k = 1:obj.ITER.PhaseII_Max
            %   update sensitivity matrix object (new model parameter)
                updateSens(obj.J,obj.m);
                %  update dHat
                obj.dHat =  obj.res+J_times_m(obj.J,obj.m_n);
                %   create/update representer matrix object
                SetRepMat(obj,obj.J);
                %   form and factor cross-product matrix for normal
                %   equations
                obj.R.formR;
                obj.R.factorR;
                %   select trade-off parameter
                selectTradeOff(obj);
                obj.iterInner(obj.ITER.niter).OccamPhase = 2;
                %  update model parameter
                obj.UpdateModel;
                %  increment iteration count
                obj.ITER.niter = obj.ITER.niter+1; 
                %  run forward model, compute misfit ....
                %          ... this could be made more efficient
                obj.ITER.RMS(obj.ITER.niter) = ForwardModel(obj);
                n = obj.ITER.niter-1;
                display(['Phase II : Iteration # ' num2str(n)])
                display(['Misfit = ' num2str(obj.ITER.RMS(n+1))])
            end
            obj.ITER.niter = obj.ITER.niter-1;         
        end
        %******************************************************************
        function HybridOcDCG(obj,K) 
        %   Hybrid DCG-Occam scheme
        %     to start testing input and use a fixed truncation level K
            obj.ITER.niter =  1;
            obj.ITER.iterMax = obj.ITER.PhaseI_Max;
            obj.ITER.RMS(1) = ForwardModel(obj);
            display(['Initial Misfit = ' num2str(obj.ITER.RMS(1))])
            %   Occam Phase I
            PHASE1_DONE = FitToTolerance(obj.ITER);
            while ~PHASE1_DONE
                %  itercontrol object for CG (inner-loop) iterations
                obj.iterInner{obj.ITER.niter} = InvIterControl;
                %   update sensitivity matrix object (new model parameter)
                %     after update J will be unnormalized ...
                updateSens(obj.J,obj.m);
                %  update dHat
                obj.dHat =  obj.res+J_times_m(obj.J,obj.m_n);
           %----->    Here are differences with Occam:
           %    1) run J_bidiag for a fixed number of steps  ... initial
           %         debugging!
                [Jproj] = J_bidiag(obj.J,obj.Cm,obj.dHat,K);
                Jprog.Orthogonalize;
           %    2) createrepresenter matrix object from projected
           %       sensitivity matrix
                SetRepMat(obj,Jproj);
           %----< end of difference with Occam
                %   form and factor cross-product matrix for normal
                %   equations
                obj.R.formR;
                obj.R.factorR;
                %   select trade-off parameter
                selectTradeOff(obj);
                obj.iterInner(obj.ITER.niter).OccamPhase = 1;
                %  update model parameter
                obj.UpdateModel;
                %  increment iteration count
                obj.ITER.niter = obj.ITER.niter+1; 
                %  run forward model, compute misfit ....
                %          ... this could be made more efficient
                obj.ITER.RMS(obj.ITER.niter) = ForwardModel(obj);
                % check fit to tolerance
                PHASE1_DONE = FitToTolerance(obj.ITER);
                %   some display ...
                n = obj.ITER.niter-1;
                display(['Phase I : Iteration # ' num2str(n)])
                display(['Misfit = ' num2str(obj.ITER.RMS(n+1))])
            end
            display('Done with Phase I')
            %    Occam Phase II
            for k = 1:obj.ITER.PhaseII_Max
            %   update sensitivity matrix object (new model parameter)
                updateSens(obj.J,obj.m);
                %  update dHat
                obj.dHat =  obj.res+J_times_m(obj.J,obj.m_n);
           %----->    Here are differences with Occam:
           %    1) run J_bidiag for a fixed number of steps  ... initial
           %         debugging!
                [Jproj] = J_bidiag(obj.J,obj.Cm,obj.dHat,K);
                Jprog.Orthogonalize;
           %    2) createrepresenter matrix object from projected
           %       sensitivity matrix
                SetRepMat(obj,Jproj);
           %----<       end of difference with Occam
                %   form and factor cross-product matrix for normal
                %   equations
                obj.R.formR;
                obj.R.factorR;
                %   select trade-off parameter
                selectTradeOff(obj);
                obj.iterInner(obj.ITER.niter).OccamPhase = 2;
                %  update model parameter
                obj.UpdateModel;
                %  increment iteration count
                obj.ITER.niter = obj.ITER.niter+1; 
                %  run forward model, compute misfit ....
                %          ... this could be made more efficient
                obj.ITER.RMS(obj.ITER.niter) = ForwardModel(obj);
                n = obj.ITER.niter-1;
                display(['Phase II : Iteration # ' num2str(n)])
                display(['Misfit = ' num2str(obj.ITER.RMS(n))])
            end
            obj.ITER.niter = obj.ITER.niter-1;         
        end
        %******************************************************************
        function [NU,RMS] = selectTradeOff(obj)
        %   selects tradeoff parameter in Occam scheme ...   at present
        %        does not do an efficient line search, but rather computes
        %        misfit for a range of tradeoff parameters, and returns the
        %        tradeoff curve, along with the optimai
            NU = SETnu(obj.R);
            RMS = zeros(size(NU));
            
            for k = 1:length(NU)      
                obj.b = solveR(obj.R,obj.dHat,NU(k));
                obj.UpdateModel;
                RMS(k) = ForwardModel(obj);
            end
            indMin = RMS==min(RMS);
            i2 = find(RMS < obj.ITER.RMStolerance,1,'last');
            if isempty(i2)
                %  use nuMin
                obj.nu = NU(indMin); 
            else
                %   use NU(i2)
                obj.nu = NU(i2);
            end
            %   recompute b for optimal trade-off parameter
            obj.b = solveR(obj.R,obj.dHat,obj.nu);
            
            %   save results: NU, X2
            obj.iterInner(obj.ITER.niter).NU = NU;
            obj.iterInner(obj.ITER.niter).X2 = RMS.^2;
        end 
        %******************************************************************
        function [NU,RMS] = selectNuMu(obj)
        %   variant on selectTradeOff, based on new idea of computing
        %    delta m instead of m (still in data space) and projecting
        %    previous m into components in the representer span and
        %    perpendicular to this space
            
            nMU = 10;
            [mPerpTilde,bTilde,mParTilde] = perpMpar(obj,obj.mTilde);
            mPerp = CovMult(obj.Cm,mPerpTilde,1);
            r =  MultCdInv(obj.res);
            NU = SETnu(obj.R);
            MU = 0:1/(nMU-1):1;
            RMS = zeros(size(NU));
            RMSmu = zeros(size(MU));
            RMSmin = inf;
            
            for k = 1:length(NU) 
                %   c_lambda is in "coefficient space"  divided by data
                %       error variance
                d = r-NU(k)*bTilde;
                c_lambda = solveR(obj.R,d,NU(k));
                dm_lambda = JT_times_d(obj.J,c_lambda);
                m1 = obj.m+CovMult(obj.Cm,dm_lambda,2);
                r1 = obj.d-fwd(m1,obj.d);
                RMS(k) = sqrt(dot(r1,r1)/length(r1));
                if RMS(k) < RMSmin
                    NuMin = NU(k);
                    mSave = m1;
                    RMSmin = RMS(k);
                    dmTildeSave = CovMult(obj.Cm,dm_lambda,1);
                    %    after this bTilde is in  "coefficient space"
                    bSave = MultCdInv(bTilde);
                    bSave.normalized = 0;
                    bSave = c_lambda+bSave;
   %             else
   %                 break
                end
            end
            RMSmu(1) = RMSmin;
            kMu = nMU;
            for k = 2:nMU
                m1 = mSave-MU(k)*mPerp;
                r1 = obj.d-fwd(m1,obj.d);
                RMSmu(k) = sqrt(dot(r1,r1)/length(r1));
                if RMSmu(k) > RMSmu(k-1) && ...
                        RMSmu(k) > obj.ITER.RMStolerance
                    %   will update obj, using NuMin, MU(k-1)
                    kMu = k-1;
                    break
                end
            end
             %   do the update of obj, using NuMin, MU(kMu)
             obj.nu = NuMin;
             obj.mu = 1-MU(kMu);
             obj.mPerpTilde = (1-MU(kMu))*mPerpTilde;
             obj.m = mSave-MU(kMu)*mPerp;
             obj.m_n = obj.m-obj.m_0;
             obj.mTilde = obj.mTilde + dmTildeSave-MU(kMu)*mPerpTilde;
             %   perhaps mostly need to do the following, then call
             %   UpdateModel ... need to check this (or maybe obj.b is
             %   now wrong!)
             obj.b = bSave;
             obj.iterInner(obj.ITER.niter).NU = NU;
             obj.iterInner(obj.ITER.niter).X2 = RMS.^2;
             obj.iterInner(obj.ITER.niter).X2mu = RMSmu(1:k).^2;
             
             obj.iterInner(obj.ITER.niter).MU = MU;
        end 

        %******************************************************************
        function [RMS] = rmsNuMu(obj)
        %  
            nMU = 10;
            [mPerp,bTilde,mPar] = perpMpar(obj,obj.mTilde);
            mPerp = CovMult(obj.Cm,mPerp,1);
            mPar = CovMult(obj.Cm,mPar,1);
            r =  MultCdInv(obj.res);
            
            NU = SETnu(obj.R);
            MU = 0:1/(nMU-1):1;
            RMS = zeros(length(NU),length(MU));
            %    try keeing nu fixed at last optimal value ... then as
            %  lambda -> inf the model change goes to zero
            d = r-obj.nu*bTilde;
            for k = 1:length(NU) 
                %d = r-NU(k)*bTilde;
                c_lambda = solveR(obj.R,d,NU(k));
                dm_lambda = JT_times_d(obj.J,c_lambda);
                m1 = obj.m_0+CovMult(obj.Cm,dm_lambda,2);
                m1 = m1 + mPar;
                for j = 1:length(MU)
                    m = m1+MU(j)*mPerp;
                    r1 = obj.d-fwd(m,obj.d);
                    RMS(k,j) = dot(r1,r1);
                end
            end
            RMS = sqrt(RMS/length(r1));
        end
        
         %******************************************************************
        function [mTildePerp,bTilde,mTildePar] = perpMpar(obj,mTilde) 
        %   using J and R from GN object obj, divide mTilde (a normalized
        %       model parameter object) into orthogonal components in the
        %       span of, and orthogonal to, the normalized sensitivities.
        
           %  assuming obj.R is formed and factored, compute truncated SVD of J
           %   this sets truncation level
           lambda0  = obj.nu/100;
           [S,V,U] = svdJ(obj.R,lambda0);
           
           mTemp = reshape(mTilde.v,obj.R.J.M,1);
           vTm = V*mTemp;
           bVec = diag(1./S)*vTm;
           bVec = U*bVec;
           mTildePar = V'*vTm;
           mTildePerp = mTemp - mTildePar;
           %   return outputs as model parameter and data vector objects
           normalized =  true;
           bTilde = SetVec(obj.d,bVec,normalized);
           %    need to sort out correct way to do this for more general
           %     case ... so far only works for DASSOC
           %bTilde = SetNormVec(obj.d,bVec,normalized);
           mTildePerp = SetVec(mTilde,mTildePerp);
           mTildePar = SetVec(mTilde,mTildePar);
        end           
        %******************************************************************
        function HybOcc(obj,Algorithm,ITER) 
        %   Hybrid Lanczos-Occam scheme
        %    Algorithm can be:   
        %         DASOCC : Data space Occam, based on full sensitivity
        %         BIDIAG : Apply Bidiag 2 to scaled sensitivity, use data
        %                   and model space vectors to solve projected
        %                   problem
        %         BDORTH : like BIDIAG, but with solution and
        %                  orthogonalization at each Lanczos step
        %         BDMTX  : like BDORTH, but using component sensitivities
        %                  for each transmitter
        %
            obj.TIMES = cell(50,1);
            nTIMES = 0;
            obj.ITER.niter =  1;
            obj.ITER.iterMax = obj.ITER.PhaseI_Max;
            obj.ITER.RMS(1) = ForwardModel(obj);
            display(['Initial Misfit = ' num2str(obj.ITER.RMS(1))])
            %   Occam Phase I
            PHASE1_DONE = FitToTolerance(obj.ITER);
            while ~PHASE1_DONE
                if nargin <3
                    %  itercontrol object for inner-loop iterations
                    obj.iterInner(obj.ITER.niter) = InvIterControl;
                else
                    
                    %   use input ITER ... could be a cell array of
                    %   iteration control objects, allowing different
                    %   parameters for each iteration.  If this case the
                    %   last cell of ITER is used for Phase II; the second
                    %   to last is used for all Phase I iterations beyond
                    %   length(ITER)-1
                    if length(ITER) == 1
                        obj.iterInner(obj.ITER.niter) = ITER;
                    else
                        n = min(length(ITER)-1,obj.ITER.niter);
                        obj.iterInner(obj.ITER.niter) = ITER{n};
                    end
                end
                tic
                %   update sensitivity matrix object (new model parameter)
                updateSens(obj.J,obj.m);
                
                   time = toc;
                   nTIMES = nTIMES+1;
                   obj.TIMES{nTIMES} = struct('step','J','T',time,'iter',obj.ITER.niter);
                   tic
                
                %  update dHat
                obj.dHat = obj.res+J_times_m(obj.J,obj.m_n);
                
                   time = toc;
                   nTIMES = nTIMES+1;
                   obj.TIMES{nTIMES} = struct('step','dHat','T',time,'iter',obj.ITER.niter);
               
                %----->    Here is where different algorithms differ ...
                switch Algorithm
                    case 'DASOCC'
                        %   create/update representer matrix object 
                            tic
                            
                        SetRepMat(obj,obj.J);
                                        
                           time = toc;
                           nTIMES = nTIMES+1;
                           obj.TIMES{nTIMES} = struct('step','set R','T',time,'iter',obj.ITER.niter);
                           
                    case 'BIDIAG'
                        %    1) run J_bidiag for a fixed number of steps  
                        %        ... initial debugging!
                        K = obj.iterInner(obj.ITER.niter).iterMax;
                        saveMTX = true;
                        kk  = obj.ITER.niter;
                        [JprojMTX{kk}] = J_bidiag(obj.J,obj.Cm,obj.dHat,K,saveMTX);
                        dHatMTX{kk} = obj.dHat;
                        J_MTX{kk} = copy(obj.J);
                        [Jproj] = collapseMTX(JprojMTX{kk});
                        Jproj.Orthogonalize;
                        %    2) createrepresenter matrix object from 
                        %       projected sensitivity matrix
                        SetRepMat(obj,Jproj);
                    case 'BIDIAGmtx'
                        %    1) run J_bidiag for a fixed number of steps  
                        %        ... initial debugging!
                        saveMTX = true;
                        K = obj.iterInner(obj.ITER.niter).iterMax;
                        [Jproj] = J_bidiag(obj.J,obj.Cm,obj.dHat,K,saveMTX);
                        Jproj.Orthogonalize;
                        %    2) createrepresenter matrix object from 
                        %       projected sensitivity matrix
                        SetRepMat(obj,Jproj);
                    case 'BDORTH'
                        %    1) run orthogonalized Lanczos to generate
                        %    projected Jacobian
                        
                           tic
                           
                        saveMTX = false;
                        K  = obj.iterInner(obj.ITER.niter).iterMax;
                        NTX = obj.d.NTX;
                        Jproj = ProjSensMTX(NTX,K,saveMTX);
                        Jproj.d = obj.dHat;
                        [~,itTemp] = Lanczos(Jproj,obj.J,obj.dHat,obj.Cm,...
                            obj.iterInner(obj.ITER.niter),obj.nu);
    
                           obj.iterInner(obj.ITER.niter)=itTemp;
                        %    2) createrepresenter matrix object from 
                        %       projected sensitivity matrix
                        SetRepMat(obj,Jproj);
                                      
                           time = toc;
                           nTIMES = nTIMES+1;
                           obj.TIMES{nTIMES} = struct('step','lANCZOS','T',time,'iter',obj.ITER.niter);
                       
                    case 'BDMTX'
                        %    1) run orthogonalized Lanczos to generate
                        %    projected Jacobian  ... but now for MTX case
                         
                           tic
                           
                        saveMTX = true;     
                        K  = obj.iterInner(obj.ITER.niter).iterMax;
                        NTX = obj.d.NTX;
                        Jproj = ProjSensMTX(NTX,K,saveMTX);
                        Jproj.d = obj.dHat;
                        [~,itTemp] = Lanczos(Jproj,obj.J,obj.dHat,obj.Cm,...
                            obj.iterInner(obj.ITER.niter),obj.nu);
                        obj.iterInner(obj.ITER.niter)=itTemp;
                        %    2) createrepresenter matrix object from 
                        %       projected sensitivity matrix
                        SetRepMat(obj,Jproj);
                        
                           time = toc;
                           nTIMES = nTIMES+1;
                           obj.TIMES{nTIMES} = struct('step','lANCZOS','T',time,'iter',obj.ITER.niter);
                end
                
                %   form and factor cross-product matrix for normal
                %   equations
                   tic
                   
                obj.R.formR;
                               
                   time = toc;
                   nTIMES = nTIMES+1;
                   obj.TIMES{nTIMES} = struct('step','form R','T',time,'iter',obj.ITER.niter);
                           
                   tic
                   
                obj.R.factorR;
                                               
                   time = toc;
                   nTIMES = nTIMES+1;
                   obj.TIMES{nTIMES} = struct('step','factor R','T',time,'iter',obj.ITER.niter);
                           
                   tic
                   
                %   select trade-off parameter
                if obj.ITER.ModifiedOccam && (obj.ITER.niter > 1)
                    obj.selectNuMu;
                else
                    obj.selectTradeOff;
                    %  update model parameter
                    obj.UpdateModel;
                end
                
                   time = toc;
                   nTIMES = nTIMES+1;
                   obj.TIMES{nTIMES} = struct('step','select mu','T',time,'iter',obj.ITER.niter);
                   
                obj.iterInner(obj.ITER.niter).OccamPhase = 1;
                %  increment iteration count
                obj.ITER.niter = obj.ITER.niter+1; 
                %  run forward model, compute misfit ....
                %          ... this could be made more efficient
                obj.ITER.RMS(obj.ITER.niter) = ForwardModel(obj);
                % check fit to tolerance
                PHASE1_DONE = FitToTolerance(obj.ITER);
                %   some display ...
                n = obj.ITER.niter-1;
                display(['Phase I : Iteration # ' num2str(n)])
                display(['Misfit = ' num2str(obj.ITER.RMS(n+1))])
            end
            display('Done with Phase I')
            %    Occam Phase II
            for k = 1:obj.ITER.PhaseII_Max
                if nargin <3
                    %  itercontrol object for inner-loop iterations
                    obj.iterInner(obj.ITER.niter) = InvIterControl;
                else
                    %   use input ITER ... could be a cell array of
                    %   iteration control objects, allowing different
                    %   parameters for each iteration.  If this case the
                    %   last cell of ITER is used for Phase II; the second
                    %   to last is used for all Phase I iterations beyond
                    %   length(ITER)-1
                    if length(ITER) == 1
                        obj.iterInner(obj.ITER.niter) = ITER;
                    else
                        obj.iterInner(obj.ITER.niter) = ITER{end};
                    end
                end
            %   update sensitivity matrix object (new model parameter)
                updateSens(obj.J,obj.m);
                %  update dHat
                obj.dHat =  obj.res+J_times_m(obj.J,obj.m_n);
                %----->    Here is where different algorithms differ ...
                switch Algorithm
                    case 'DASOCC'
                        %   create/update representer matrix object 
                        SetRepMat(obj,obj.J);
                    case 'BIDIAG'
                        % saving JprojMTX objects ... temporary 
                         save JprojMTXtest.mat JprojMTX J_MTX dHatMTX
                        %    1) run J_bidiag for a fixed number of steps  
                        %        ... initial debugging!
                        K = obj.iterInner(obj.ITER.niter).iterMax;
                        [Jproj] = J_bidiag(obj.J,obj.Cm,obj.dHat,K);
                        Jproj.Orthogonalize;
                        %    2) createrepresenter matrix object from 
                        %       projected sensitivity matrix
                        SetRepMat(obj,Jproj);
                    case 'BIDIAGmtx'
                        %    1) run J_bidiag for a fixed number of steps  
                        %        ... initial debugging!
                        saveMTX = true;
                        K = obj.iterInner(obj.ITER.niter).iterMax;
                        [Jproj] = J_bidiag(obj.J,obj.Cm,obj.dHat,K,saveMTX);
                        Jproj.Orthogonalize;
                        %    2) createrepresenter matrix object from 
                        %       projected sensitivity matrix
                        SetRepMat(obj,Jproj);
                    case 'BDORTH'
                        %    1) run orthogonalized Lanczos to generate
                        %    projected Jacobian
                        saveMTX = false;
                        K  = obj.iterInner(obj.ITER.niter).iterMax;
                        NTX = obj.d.NTX;
                        Jproj = ProjSensMTX(NTX,K,saveMTX);
                        Jproj.d = obj.dHat;
                        [~,itTemp] = Lanczos(Jproj,obj.J,obj.dHat,obj.Cm,...
                            obj.iterInner(obj.ITER.niter),obj.nu);
                        obj.iterInner(obj.ITER.niter)=itTemp;
                        %    2) createrepresenter matrix object from 
                        %       projected sensitivity matrix
                        SetRepMat(obj,Jproj);
                    case 'BDMTX'
                        %    1) run orthogonalized Lanczos to generate
                        %    projected Jacobian  ... but now for MTX case
                        saveMTX = true;
                        K  = obj.iterInner(obj.ITER.niter).iterMax;
                        NTX = obj.d.NTX;
                        Jproj = ProjSensMTX(NTX,K,saveMTX);
                        Jproj.d = obj.dHat;
                        [~,itTemp] = Lanczos(Jproj,obj.J,obj.dHat,obj.Cm,...         
                            obj.iterInner(obj.ITER.niter),obj.nu);
                        obj.iterInner(obj.ITER.niter)=itTemp;
                        %    2) createrepresenter matrix object from 
                        %       projected sensitivity matrix
                        SetRepMat(obj,Jproj);
                end 

                %   form and factor cross-product matrix for normal
                %   equations
                obj.R.formR;
                obj.R.factorR;
                %   select trade-off parameter
                if obj.ITER.ModifiedOccam && (obj.ITER.niter > 1)
                    obj.selectNuMu;
                else
                    obj.selectTradeOff;
                    %  update model parameter
                    obj.UpdateModel;
                end
                obj.iterInner(obj.ITER.niter).OccamPhase = 2; 
                %  increment iteration count
                obj.ITER.niter = obj.ITER.niter+1; 
                %  run forward model, compute misfit ....
                %          ... this could be made more efficient
                obj.ITER.RMS(obj.ITER.niter) = ForwardModel(obj);
                n = obj.ITER.niter-1;
                display(['Phase II : Iteration # ' num2str(n)])
                display(['Misfit = ' num2str(obj.ITER.RMS(n+1))])
            end
            obj.ITER.niter = obj.ITER.niter-1;   
        end
    end 
end
    
