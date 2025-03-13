classdef MT3DmodelCov < ModelCovariance
        % model covariance class for 3DMT ... just a quick and easy example

    properties (SetAccess = public)
        rho = .40;
        zMax = 10000;                 %   what do these parameters mean???     
        nOuter = 3;
        nYinner = 2;
        nXinner = 2;
        nZinner = 1;
        zeta = .05;
    end
    properties %(Access = private) 
        Z0;
        Zm;
        Zp
        
        W0_Y;
        WmY;
        WpY;
        W0_X;
        WmX;
        WpX;
    end
    
    methods
        function obj = MT3DmodelCov(grid)
        %   class constructor ... attach grid handle, set default
        %        covariance parameter values, and setup horizontal and
        %        verical smoothing matrices
            obj.grid = grid;
            CovSet(obj);
        end
        %******************************************************************
        function  CovSet(obj)       
        %   construct horizontal and vertical smoothing matrices, using
        %   current values of model covariance parameters
            dy = [obj.grid.dy(1) obj.grid.dy' obj.grid.dy(end)];
            Cy = (dy(1:end-1)+dy(2:end))/2;
            dx = [obj.grid.dx(1) obj.grid.dx' obj.grid.dx(end)];
            Cx = (dx(1:end-1)+dx(2:end))/2;
            CyMin = min(Cy);
            CxMin = min(Cx);
            alphaMax_y = CyMin^2/6;
            alphaMax_x = CxMin^2/6;         
            alphaMin_y = obj.rho*alphaMax_y;
            alphaMin_x = obj.rho*alphaMax_x;
            z = [0 cumsum(obj.grid.dz)'];
            z = (z(1:end-1)+z(2:end))/2;
            zWt = min(1,z/obj.zMax);
            alpha_y = zWt*alphaMax_y+(1-zWt)*alphaMin_y;
            alpha_x = zWt*alphaMax_x+(1-zWt)*alphaMin_x;
            %  set up depth-dependent weights for smoothing horizontally ... 
            %             initially based on diffusion
            obj.WmY = (2./(Cy(1:end-1).*(Cy(1:end-1)+Cy(2:end))))'*alpha_y;
            obj.WmY(1,:) = 0;
            obj.WpY = (2./(Cy(2:end).*(Cy(1:end-1)+Cy(2:end))))'*alpha_y;
            obj.WpY(end,:) = 0;
            %   but let's symmetrize!
            obj.WmY = (obj.WmY(2:end,:)+obj.WpY(1:end-1,:))/2;
            obj.WpY = obj.WmY;
            obj.W0_Y = 1 - [zeros(1,obj.grid.NzEarth) ; obj.WmY]- ...
                [obj.WpY; zeros(1,obj.grid.NzEarth)];
            
            %     now smoothing in x-direction
            obj.WmX = (2./(Cx(1:end-1).*(Cx(1:end-1)+Cx(2:end))))'*alpha_x;
            obj.WmX(1,:) = 0;
            obj.WpX = (2./(Cx(2:end).*(Cx(1:end-1)+Cx(2:end))))'*alpha_x;
            obj.WpX(end,:) = 0;
            
            %   but let's symmetrize!
            obj.WmX = (obj.WmX(2:end,:)+obj.WpX(1:end-1,:))/2;
            obj.WpX = obj.WmX;
            obj.W0_X = 1 - [zeros(1,obj.grid.NzEarth) ; obj.WmX]- ...
                [obj.WpX; zeros(1,obj.grid.NzEarth)];
            %  set up weights for smoothing vertically
            obj.Zm = obj.zeta*ones(obj.grid.Nx,obj.grid.NzEarth);
            obj.Zp = obj.Zm;
            obj.Zm(:,1) = 0;
            obj.Zp(:,end) = 0;
            obj.Z0 = 1-obj.Zm-obj.Zp;
            obj.Zm = obj.Zm(:,2:end);
            obj.Zp = obj.Zp(:,1:end-1);
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
            nSmooth = obj.nOuter*nTimes;
            %  outer loop: smooth horizontally, then vertical nSmooth times
            for k = 1:nSmooth
                %   1D smooth in x-direction
                for j = 1:obj.nXinner
                    for  iy = 1:obj.grid.Ny
                        mOut.v(:,iy,:) = obj.W0_X.*squeeze(mIn.v(:,iy,:));
                        mOut.v(2:end,iy,:) = squeeze(mOut.v(2:end,iy,:))...
                            +obj.WmX.*squeeze(mIn.v(1:end-1,iy,:));
                        mOut.v(1:end-1,iy,:) = squeeze(mOut.v(1:end-1,iy,:))...
                            +obj.WpX.*squeeze(mIn.v(2:end,iy,:));
                    end
                    mIn = mOut;
                end
                for j = 1:obj.nYinner
                    for  ix = 1:obj.grid.Nx
                        mOut.v(ix,:,:) = obj.W0_Y.*squeeze(mIn.v(ix,:,:));
                        mOut.v(ix,2:end,:) = squeeze(mOut.v(ix,2:end,:))...
                            +obj.WmY.*squeeze(mIn.v(ix,1:end-1,:));
                        mOut.v(ix,1:end-1,:) = squeeze(mOut.v(ix,1:end-1,:))...
                            +obj.WpY.*squeeze(mIn.v(ix,2:end,:));
                    end
                    mIn = mOut;
                end
                for j = 1:obj.nZinner
                    for iy = 1:obj.grid.Ny
                        mOut.v(:,iy,:) = obj.Z0.*squeeze(mIn.v(:,iy,:));
                        mOut.v(:,iy,2:end) = squeeze(mOut.v(:,iy,2:end))...
                            +obj.Zm.*squeeze(mIn.v(:,iy,1:end-1));
                        mOut.v(:,iy,1:end-1) = squeeze(mOut.v(:,iy,1:end-1))...
                            +obj.Zp.*squeeze(mIn.v(:,iy,2:end));
                    end
                    mIn = mOut;
                end
                for j = 1:obj.nZinner
                    for iy = 1:obj.grid.Ny
                        mOut.v(:,iy,:) = obj.Z0.*squeeze(mIn.v(:,iy,:));
                        mOut.v(:,iy,2:end) = squeeze(mOut.v(:,iy,2:end))...
                            +obj.Zm.*squeeze(mIn.v(:,iy,1:end-1));
                        mOut.v(:,iy,1:end-1) = squeeze(mOut.v(:,iy,1:end-1))...
                            +obj.Zp.*squeeze(mIn.v(:,iy,2:end));
                    end
                    mIn = mOut;
                end
                for j = 1:obj.nYinner
                    for  ix = 1:obj.grid.Nx
                        mOut.v(ix,:,:) = obj.W0_Y.*squeeze(mIn.v(ix,:,:));
                        mOut.v(ix,2:end,:) = squeeze(mOut.v(ix,2:end,:))...
                            +obj.WmY.*squeeze(mIn.v(ix,1:end-1,:));
                        mOut.v(ix,1:end-1,:) = squeeze(mOut.v(ix,1:end-1,:))...
                            +obj.WpY.*squeeze(mIn.v(ix,2:end,:));
                    end
                    mIn = mOut;
                end
                for j = 1:obj.nXinner
                    for  iy = 1:obj.grid.Ny
                        mOut.v(:,iy,:) = obj.W0_X.*squeeze(mIn.v(:,iy,:));
                        mOut.v(2:end,iy,:) = squeeze(mOut.v(2:end,iy,:))...
                            +obj.WmX.*squeeze(mIn.v(1:end-1,iy,:));
                        mOut.v(1:end-1,iy,:) = squeeze(mOut.v(1:end-1,iy,:))...
                            +obj.WpX.*squeeze(mIn.v(2:end,iy,:));
                    end
                    mIn = mOut;
                end
            end  
        end
    end    % methods
end   % classdef
