classdef MT2D_laplace < ModelCovariance
        % model covariance class for MT2D; simple, quick implmentation!
    properties %(SetAccess = private)
        beta = 1;    %  beta is relative scaling between vertical and 
                     %    horizontal derivatives in Laplacian
        gridMetric=false;    %set gridMetric = true to use actual variable 
                     %   grid spacing for derivatives in Laplacian; else
                     %   ignore actual grid spacing
        w = 1e-6;      %   weight for mean
    end
    properties (SetAccess = private)
       Coeff;
       R;
    end
    
    methods
        function obj = MT2D_laplace(grid,beta,gridMetric)
            %   class constructor ... just attach grid handle, set default
            %        smoothing parameter values
            if nargin > 1
                obj.beta = beta;
            end
            if nargin == 3;
                obj.gridMetric = gridMetric;
            end
            obj.grid = grid;

            %   construct coefficient matrix for Poisson equation
            Nz = grid.Nz-grid.Nza;
            Ny = grid.Ny;
            I=reshape(ones(Nz,1)*[1:Ny],Ny*Nz,1);
            J=reshape((ones(Ny,1)*[1:Nz])',Ny* Nz,1);
            if obj.gridMetric
                error('Not coded for this case yet')
            else
                dy2 = 1;
                dz2 = dy2*obj.beta;
                Nyz = length(I);
                d = [-Nz -1 0 1 Nz];
                C = [(1./dy2)*ones(Nyz,1) (1./dz2)*ones(Nyz,1)  ...
                    -2*(1./dy2+1./dz2)*ones(Nyz,1)...
                    (1./dz2)*ones(Nyz,1) (1./dy2)*ones(Nyz,1)];
                C(1:Nz,5) = 0;
                C(1,4) = 0;
                C(end,2) = 0;
                C(end-Nz+1:end,1) = 0;
                
                %   left edge
                ind = find(I-1 == 0);
                C(ind,3) = C(ind,3)+1./dy2;
                
                %   right edge
                ind = find(I+1 > Ny);
                C(ind,3) = C(ind,3)+1./dy2;
                
                %   top edge
                ind = find(J-1 ==0);
                C(ind,3) = C(ind,3)+1./dz2;
                ind = ind-1;
                
                ind = ind(find(ind>0));
                C(ind,2) = 0;
                
                %   bottom edge
                ind = find(J+1 > Nz);
                C(ind,3) = C(ind,3)+1./dz2;
                ind = ind+1;
                ind = ind(find(ind<Nyz));
                C(ind,4) = 0;
                
            end
            obj.Coeff = spdiags(-C,d,Nyz,Nyz);
            %   cholesky decomp of singular matrix A; R is sparse, and
            %          R'*R  = -A
            [obj.R,p] = chol(obj.Coeff);
        end
        %******************************************************************
        function mOut = CovMult(obj,mIn,nTimes)
        %  Usage  : mout = CovMult(obj,mIn,nTimes)
        %  Input : mIn = input conductivity structure
        %          obj = model covariance object
        %          nTimes = optional argument: number of times to apply 
        %		           the smoother (default is 1)
        %  Output : mOut = smoothed output conductivity structure
        %     IN THIS CLASS: the model parameter is smoothed by solving
        %     Poisson's equation with Neumann BC; the mean is extracted from
        %       the input; the mean of the output is this value, times a
        %       scale factor obj.w
            
            if nargin < 3
                nTimes = 1;
            end
            
            %initialize output
            mOut = mIn;
            
            %  solve Poisson equation
            [ny,nz] = size(mIn.v);
            x = reshape(mIn.v',ny*nz,1);
            xBar = mean(x);
            y = obj.R\(obj.R'\x);
            if nTimes > 1
                y = obj.R\(obj.R'\y);
                y = y - mean(y) + ((1./obj.w)^2)*xBar;
            else
                y = y - mean(y) + (1./obj.w)*xBar;
            end
            mOut.v = (reshape(y,nz,ny))';
        end
  
    end    % methods
end   % classdef