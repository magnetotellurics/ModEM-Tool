classdef MT2DmodelCov < ModelCovariance
        % model covariance class for MT2D; simple, quick implmentation!
    properties %(SetAccess = private)
        rho = .2;
        zeta = .25;
        zMax = 10000;
        nOuter = 4;   % 4
        nYinner = 4;  % 4
        nZinner = 2;   %2
    end
    properties %(SetAccess = private)
        Z0;
        Zm;
        Zp
        W0;
        Wm;
        Wp;
    end
    
    methods
        function obj = MT2DmodelCov(grid)
            %   class constructor ... just attach grid handle, set default
            %        smoothing parameter values
            obj.grid = grid;
            %   construct horizontal smoothing weights
            Dy = [obj.grid.Dy(1) obj.grid.Dy obj.grid.Dy(end)];
            Cy = (Dy(1:end-1)+Dy(2:end))/2;
            Cmin = min(Cy);
            alphaMax = Cmin^2/6;
            alphaMin = obj.rho*alphaMax;
            z = [0 cumsum(obj.grid.Dz(obj.grid.Nza+1:end))];
            z = (z(1:end-1)+z(2:end))/2;
            zWt = min(1,z/obj.zMax);
            alpha = zWt*alphaMax+(1-zWt)*alphaMin;
            %  set up weights for smoothing horizontally ... initially based on diffusion
            obj.Wm = (2./(Cy(1:end-1).*(Cy(1:end-1)+Cy(2:end))))'*alpha;
            obj.Wm(1,:) = 0;
            obj.Wp = (2./(Cy(2:end).*(Cy(1:end-1)+Cy(2:end))))'*alpha;
            obj.Wp(end,:) = 0;
            obj.W0 = 1 - obj.Wm-obj.Wp;
            %   but let's symmetrize!
            obj.Wm = (obj.Wm(2:end,:)+obj.Wp(1:end-1,:))/2;
            obj.Wp = obj.Wm;
            %  set up weights for smoothing vertically
            obj.Zm = obj.zeta*ones(grid.Ny,grid.Nz-grid.Nza);
            obj.Zp = obj.Zm;
            obj.Zm(:,1) = 0;
            obj.Zp(:,end) = 0;
            obj.Z0 = 1-obj.Zm-obj.Zp;
            obj.Zm = obj.Zm(:,2:end);
            obj.Zp = obj.Zp(:,1:end-1);
        end
        %******************************************************************
        function mOut = CovMult(obj,mIn,nTimes)
        %  Usage  : mout = CovMult(mIn,smthParams)
        %  Input : mIn = input conductivity structure
        %          obj = model covariance object
        %          nTimes = optional argument: number of times to apply 
        %		           the smoother (default is 1)
        %  Output : mOut = smoothed output conductivity structure
            
            if nargin < 3
                nTimes = 1;
            end
            
            %initialize output
            mOut = mIn;
            nSmooth = obj.nOuter*nTimes;
            %  outer loop: smooth horizontally, then vertical nSmooth times
            for k = 1:nSmooth
                for j = 1:obj.nYinner
                    mOut.v = obj.W0.*mIn.v;
                    mOut.v(2:end,:) = mOut.v(2:end,:)+obj.Wm.*mIn.v(1:end-1,:);
                    mOut.v(1:end-1,:) = mOut.v(1:end-1,:)+obj.Wp.*mIn.v(2:end,:);
                    mIn = mOut;
                end
                % skip vertical smoothing on last step to insure symmetry
                if k < nSmooth
                    mIn = mOut;
                    for j = 1:obj.nZinner
                        mOut.v = obj.Z0.*mIn.v;
                        mOut.v(:,2:end) =  ...
                            mOut.v(:,2:end)+obj.Zm.*mIn.v(:,1:end-1);
                        mOut.v(:,1:end-1) = ...
                            mOut.v(:,1:end-1)+obj.Zp.*mIn.v(:,2:end);
                    end
                    mIn = mOut;
                end
            end
        end
  
    end    % methods
end   % classdef
