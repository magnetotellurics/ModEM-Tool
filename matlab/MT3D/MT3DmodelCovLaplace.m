classdef MT3DmodelCovLaplace < ModelCovariance
    % model covariance class for 3DMT using R^{-1} = \nabla^2.
    % Then solve Poisson's equation as R m = \tilde{m} for m.
    
    properties (SetAccess = public)
        opLaplace   %  anonymous function handle ... if we can make this work
        R   %   storage for Incomplete Cholesky decomp
        tol = 1e-6;
        maxit = 100;
        eps = 1000;
        O
    end
    properties %(Access = private)
        
    end
    
    methods
        function obj = MT3DmodelCovLaplace(grid)
            %   class constructor ... attach grid handle, set default
            %        covariance parameter values, and setup horizontal and
            %        verical smoothing matrices
            obj.grid = grid;
            CovSet(obj);
        end
        %******************************************************************
        function  CovSet(obj)
            %   construct the opLaplacian operator
            % construct the gradient operator first.
            %   note that model parameter grid is not a TGrid3D object
            gr = TGrid3D(obj.grid.dx,obj.grid.dy,obj.grid.dz,obj.grid.NzAir);
            op = TOperatorTopology_SG(gr);
            op.setGrad('cells2faces');
            % construct the standard Laplacian operator.
            faceObj = TVector3D_SG(gr,'face');
            [ii,~] = faceObj.int_bdry_indices;
            obj.opLaplace = op.G(ii,:)'*op.G(ii,:);
            obj.O = ones(gr.NCells,1)/sqrt(gr.NCells);
            %   compute and save incomplete Cholesky decomp with no fill-in
            %     this is what would is coded in ModEM for divergence
            %     correction now ... nopt sure how to do this here
            %obj.R = cholinc(obj.GG+eps*speye(gr.NCells),'0');
            % set the Nuemann's boundary condition
        end
        %******************************************************************
        function mOut = CovMult(obj,mIn,nTimes)
            %  Usage  : mout = CovMult(sIn,smthParams)
            %  Input : mIn = input conductivity structure
            %          obj = model covariance object
            %          nTimes = optional argument: number of times to apply
            %		           the smoother (default is 1)
            %  Output : sOut = smoothed output conductivity structure
            
            if nargin < 3
                nTimes = 1;
            end
            
            %initialize output
            mOut = mIn;
            %   smooth nTimes ...  normally 1 or 2
            for n = 1:nTimes
                rhs = mOut.ExtractVec;
                % solve Poisson's equation.
                [x,flag,relres,iter] = pcg(@modLap,rhs,obj.tol,obj.maxit);
                mOut = mOut.SetVec(x);
            end
            function y = modLap(x)
                y = obj.opLaplace*x + obj.eps*obj.O*(obj.O'*x);
            end
        end
    end    % methods
    methods (Static)
        
    end
end   % classdef
