classdef MT3DmodelCovPoisson < ModelCovariance
        % model covariance class for 3DMT using R^{-1} = \nabla^2.
		% Then solve Poisson's equation as R m = \tilde{m} for m.
        
    properties (SetAccess = public)
        opLaplace
    end
    properties %(Access = private) 
        
    end
    
    methods
        function obj = MT3DmodelCovPoisson(grid)
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
			gr = TGrid3D(obj.grid.dx,obj.grid.dy,obj.grid.dz,obj.grid.NzAir);
			op = TOperatorTopology_SG(gr);
			op.setGrad('cells2innerfaces');
			% construct the standard Laplacian operator.
			obj.opLaplace = - op.G'*op.G;	
			% set the Nuemann's boundary condition
			% as \partial{\sigma} / \partial n = 0.
			nCells = gr.NCells;
			nx     = gr.Nx;
			ny     = gr.Ny;
			nz     = gr.Nz;
			% diag elements.
			for i = 1:nCells
				obj.opLaplace(i,i) = -6;
			end % i.
			% loop over all boundary cells.
			% i = 1, nx;
			for J = 1:ny
				for K = 1:nz
					I = 1;
					nRow = gr.vectorIndex(I,J,K,'cell');
					nCol = gr.vectorIndex(I+1,J,K,'cell');
					obj.opLaplace(nRow,nCol) = 2;
					I = nx;
					nRow = gr.vectorIndex(I,J,K,'cell');
					nCol = gr.vectorIndex(I-1,J,K,'cell');
					obj.opLaplace(nRow,nCol) = 2;
				end % k.
			end % j.
			% j = 1, ny;
			for I = 1:nx
				for K = 1:nz
					J = 1;
					nRow = gr.vectorIndex(I,J,K,'cell');
					nCol = gr.vectorIndex(I,J+1,K,'cell');
					obj.opLaplace(nRow,nCol) = 2;
					J = ny;
					nRow = gr.vectorIndex(I,J,K,'cell');
					nCol = gr.vectorIndex(I,J-1,K,'cell');
					obj.opLaplace(nRow,nCol) = 2;
				end % k.
			end % i.
			% k = 1, nz;
			for I = 1:nx
				for J = 1:ny
					K = 1;
					nRow = gr.vectorIndex(I,J,K,'cell');
					nCol = gr.vectorIndex(I,J,K+1,'cell');
					obj.opLaplace(nRow,nCol) = 2;
					K = nz;
					nRow = gr.vectorIndex(I,J,K,'cell');
					nCol = gr.vectorIndex(I,J,K-1,'cell');
					obj.opLaplace(nRow,nCol) = 2;
				end % J.
			end % I.
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
			% form the rhs.
			% A UGLY solution for getting the GRID size without changing 
			% the inconsistent definition.
			try
				nx = mIn.grid.nx;
				ny = mIn.grid.ny;
				nz = mIn.grid.nzEarth;
			end
			try
				nx = mIn.grid.Nx;
				ny = mIn.grid.Ny;
				nz = mIn.grid.NzEarth;
			end
			nCell = nx * ny * nz;
			rhs = reshape(mIn.v,nCell,1);
			% save Poisson's equation.
			tol = 1e-5;
			x = qmr(obj.opLaplace,rhs,tol);
			mOut.v = reshape(x,nx,ny,nz);
		end
    end    % methods
end   % classdef
