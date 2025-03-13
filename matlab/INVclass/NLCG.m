classdef NLCG < handle
    % class for for NLCG-type inversion algorithms
    %
    %   Gary D. Egbert, 2011
    %   College of Oceanic and Atmospheric Sciences
    %   Oregon State University
    
    properties
        J       %   sensitivity matrix handle
        Cm      %   model covariance handle
        m_0     %   prior model parameter (value object)
        m       %   current model parameter (value object)
        m_n     %   model parameter deviation from prior: m_n = m-m_0
        mTilde  %    normalized deviation from prior
        d       %   data (data space value object)
        h       %   search direction in normalized model space
        mTrial  %    trial model parameter (normalized)
        g       %   gradient in normalized model space
        pred    %   f(m)  (data space value object)
        nCG     %   number of steps since last restart
        PEN     %   array of penalty function values
        RMS     %   array of RMS misfits
        ITER    %   InvIterControl object for outer loop
        nu      %   damping parameter
        beta    %   used to compute new search direction
    end
    properties (Dependent=true)
        res     %   residual d - f(m)  (data space value object)
    end
    
    methods
        %******************************************************************
        function obj = NLCG(d,J,nu,m_0,Cm,m)
        % class constructor
            obj.d = d;
            obj.J = J;
            obj.nu = nu;
            obj.m_0 = m_0;
            obj. Cm = Cm;
            obj.ITER = InvIterControl();
            if nargin == 5
                obj.m = m_0;
                obj.m_n = m_0;
                obj.m_n = ZeroVec(obj.m_n);
                obj.mTilde = obj.m_n;
            else
                obj.m = m;
                obj.m_n = m;
                obj.m_n = obj.m-m_0;
                obj.mTilde = ZeroVec(obj.m_n);
            end
        end
        %******************************************************************
        function res = get.res(obj)
            res = obj.d-obj.pred;
        end
        %******************************************************************
        function RMSmisfit  = ForwardModel(obj)
        %   using current m_n update pred, res, dHat
            obj.pred = fwd(obj.m,obj.d);
%             obj.res = obj.d-obj.pred;
            RMSmisfit = sqrt(real(dot(obj.res,obj.res))/obj.res.length);
        end
        %******************************************************************
        function [P,m1] = PenaltyFunctional(obj)
        %   evaluate penalty functional for trial model parameter
        %   obj.mTrial
            m1 = obj.m_0+CovMult(obj.Cm,obj.mTrial);
            obj.pred = fwd(m1,obj.d);
%             r = obj.d-obj.pred;
            r = obj.res;
            P = dot(r,r) + obj.nu*dot(obj.mTrial,obj.mTrial);
        end
        %******************************************************************
        function [] = gradient(obj)
        %   evaluate penalty functional for trial model parameter
        %   obj.mTrial
        %  initialize/update sensitvity
            updateSens(obj.J,obj.m);
            %  compute initial gradient
            r = Normalize(obj.res,2);
            obj.g = JT_times_d(obj.J,r);
            obj.g = CovMult(obj.Cm,obj.g,1);
            %  gradient in preconditioned space
            nu2 = obj.nu*2;
            obj.g = linComb(-2,obj.g,nu2,obj.mTilde);
        end
        %******************************************************************
        function OK = TestOrthog(obj)
            %   my old approach
            g_orth_h = dot(obj.g,obj.h)/dot(obj.g,obj.g);
            OK = obj.beta*g_orth_h < .5;
            %   Anya's approach
%             OK = obj.
        end
        %******************************************************************
        function [pen] = evalPEN(obj)
        %   evaluate penalty functional.
            pen = dot(obj.res,obj.res) + obj.nu * dot(obj.mTilde,obj.mTilde);
        end
        %******************************************************************
        function NLCGsolve(obj)
        %  NLCG solver routine.  Derived from old INV/NLCG
            
            %  initialize models, residuals
            obj.RMS(1) = ForwardModel(obj);
            obj.ITER.PEN(1) = evalPEN(obj);
            obj.ITER.niter = 1;
            obj.gradient;
            obj.h = -1*obj.g;   
            obj.nCG = 0;
            while  true
                obj.lineSearch;
                obj.RMS(obj.ITER.niter) = sqrt(dot(obj.res,obj.res)/obj.res.length);
                if obj.ITER.niter > 2
					if (abs(obj.RMS(obj.ITER.niter)-obj.RMS(obj.ITER.niter-1)) < 2e-3)
						obj.nu = obj.nu/2;
						disp(['Damping parameter reduced to: ' num2str(obj.nu)])
					end
				end
				obj.ITER.PEN(obj.ITER.niter) = evalPEN(obj);
                if obj.NLCGconverged
                    break
                end
                %   Save last gradient
                gOld = obj.g;
                %  compute next gradient
                obj.gradient; 
                obj.beta = (dot(obj.g,obj.g)-dot(gOld,obj.g))/dot(gOld,gOld);
                %   if orthogonality breaks down, restart
                if obj.TestOrthog  && obj.nCG < obj.ITER.nCGmax
                    obj.h = linComb(-1.,obj.g,obj.beta,obj.h);
                    obj.nCG = obj.nCG + 1;
                else
                    %  restart
                    obj.h = obj.g;
                    obj.h = -1*obj.h;
                    obj.nCG = 0;
                end
                disp('===================================================');
                disp(['Iter: ' num2str(obj.ITER.niter) ' RMS: ' num2str(obj.RMS(obj.ITER.niter))]);
                disp('===================================================');
                obj.ITER.niter = obj.ITER.niter+1;
            end
        end
        
        %******************************************************************
        function TorF = NLCGconverged(obj)
            %   test convergence 
% `           RMS = sqrt(obj.SS(obj.ITER.niter)/obj.d.length);
            TorF = obj.RMS(obj.ITER.niter) < obj.ITER.RMStolerance || obj.ITER.niter > obj.ITER.iterMax;
            %   I also had this additional condition
            %   TorF = TorF || obj.ITER.penChange < obj.ITER.minchange
        end
        %******************************************************************
        function lineSearch(obj)
        %
        %   Line search for conjugate gradients, following simple scheme
        %    described in Newman and Alumbaugh, using model transformation
        %    approach for preconditioning
            
            delta = 1e-4;
            
            niter = obj.ITER.niter;
            pen = obj.ITER.PEN(niter);
            
            beta = sqrt(dot(obj.h,obj.h));
            obj.h = (1./beta)*obj.h;
            gamma = dot(obj.h,obj.g);
            alpha_Trial = 1./ModelParamMax(obj.h);
            alphaMin = alpha_Trial*.1;
            
            obj.mTrial = linComb(1.,obj.mTilde,alpha_Trial,obj.h);
            PenTrial = PenaltyFunctional(obj);
            
            while PenTrial > pen+delta*gamma && alpha_Trial > alphaMin
                alphaMin = alpha_Trial*.1;
                alpha_Trial = - gamma*alpha_Trial*alpha_Trial/...
                    (2*(PenTrial-pen-gamma*alpha_Trial));
                alpha_Trial = max(alpha_Trial,alphaMin);
                obj.mTrial = linComb(1.,obj.mTilde,alpha_Trial,obj.h);
                PenTrial = PenaltyFunctional(obj);
            end
            
            b = ((PenTrial - pen) - gamma*alpha_Trial)/(alpha_Trial^2);
            alpha = -gamma/(2*b);
            if b < 0
                obj.mTilde = obj.mTrial;
                obj.m = obj.m_0+CovMult(obj.Cm,obj.mTrial);
                obj.ITER.PEN(niter+1) = PenTrial;
                return
            end
            obj.mTilde =  linComb(1.,obj.mTilde,alpha,obj.h);
            [Pen1,m_1] = PenaltyFunctional(obj);
           
            
            if Pen1 > PenTrial
                obj.mTilde = obj.mTrial;
                obj.m = obj.m_0+CovMult(obj.Cm,obj.mTrial);
                obj.ITER.PEN(niter+1) = PenTrial;
                return
            else
                obj.m = m_1;
                obj.ITER.PEN(niter+1) = Pen1;
            end
            
        end
    end
end
