classdef xygrid < latlontools
    %   class to define a regular x/y grid
    properties
        dx
        dy
        dz
        origin = [0; 0; 0];
        rotation = 0;
        units = 'km';
        nzAir = 12;
        zAir = []; % +ve height above ground, length nzAir+1
        xpadding;
        ypadding;
        mstruct; % save map projection here
    end
    properties (SetAccess = public)
        % 2D mesh of Y/X values; regular if computed from dx & dy;
        % irregular if obtained by interpolation from a lat/lon grid
        %Y
        %X
    end
    properties (Dependent = true)
        nzEarth
        nx
        ny
    end
    properties (Dependent = true)
        % for backward compatibility
        NzEarth
        Nx
        Ny
    end
    methods
        function [obj] = xygrid(dx,dy,dz,origin,rotation,units,nzAir,zAir)
            %   [obj] = xygrid(dx,dy,dz,origin,rotation,units,nzAir,zAir)
            %
            %   class constructor; if the origin is exactly zero assume
            %   that it's not supplied and use the grid center as origin
            % 
            %   to convert an xy-type ModEMM TGrid3D to xygrid, use the
            %   function modemm2xygrid
            if nargin == 0
                obj.dx = [];
                obj.dy = [];
                obj.dz = [];
                return
            end
            obj.dx = dx;
            obj.dy = dy;
            obj.dz = dz;
            if nargin < 4
                obj.origin(1) = - sum(dx)/2;
                obj.origin(2) = - sum(dy)/2;
                obj.origin(3) = 0;
            else
                obj.origin = origin;
            end
            if nargin > 4
                obj.rotation = rotation;
            end
            if nargin > 5
                obj.units = units;
            end
            if nargin > 6
                obj.nzAir = nzAir;
            end
            if nargin > 7
                obj.zAir = zAir;
            end
            %[obj.Y,obj.X] = meshgrid(obj.y,obj.x);
            obj.xpadding = 0;
            obj.ypadding = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % map X/Y grid to lat/lon locations at centers of the cells
        function grid = llgrid(obj,mstruct)
            % grid = llgrid(obj,mstruct)
            
            lat0 = mstruct.origin(1);
            lon0 = mstruct.origin(2);
                    
            % convert all distances to meters for formal map projections
            inUnits = obj.units;
            obj = obj.km2m;
            
            top = min(obj.z);
            bottom = max(obj.z);
            gridx = obj.x;
            gridy = obj.y;   
            %[X,Y]=meshgrid(obj.x,obj.y);
            
            % compute grid lat/lon for the limits and dlat/dlon
            switch strtrim(mstruct.mapprojection)
                case 'xy2latlon'
                    [LAT,LON] = latlon(obj,lat0,lon0); % can deal with km
                    
                case 'utm'
                    % convert projection origin to Cartesian
                    [y0,x0] = mfwdtran(mstruct,lat0,lon0);
                    oX = repmat(x0,length(gridy),1);
                    oY = repmat(y0,length(gridx),1);
                    
                    % convert grid nodes to Spherical
                    [~,LON] = minvtran(mstruct,gridy,oX);
                    [LAT,~] = minvtran(mstruct,oY,gridx);

                    %[LAT,LON] = minvtran(mstruct,Y,X);
                    
                otherwise % e.g.,'eqdcylin','eqacylin','lambertstd','eqaazim'                   
                    % convert projection origin to Cartesian
                    [y0,x0] = projfwd(mstruct,lat0,lon0);
                    oX = repmat(x0,length(gridy),1);
                    oY = repmat(y0,length(gridx),1);
                    
                    % convert grid nodes to Spherical
                    [~,LON] = projinv(mstruct,gridy,oX);
                    [LAT,~] = projinv(mstruct,oY,gridx);

                    %[LAT,LON] = projinv(mstruct,Y,X);
                    
            end
            
            lims.latmin = min(min(LAT));
            lims.latmax = max(max(LAT));
            lims.lonmin = min(min(LON));
            lims.lonmax = max(max(LON));
            lims.depthmin = top;
            lims.depthmax = bottom;

            grid = llgrid(lims,top,bottom);

            if  strcmp(mstruct.mapprojection,'xy2latlon')
                %THIS METHOD IS VERY INACCURATE IN SOME CIRCUMSTANCES - DO NOT
                %USE - MUCH BETTER TO USE THE INTERPOLATION BELOW
                %grid.lat = LAT(:,1);
                %ind = find(LAT(:,1)<=mean(grid.lat), 1, 'last'); % use lon for the average lat
                % to map back to the
                % original lat/lon grid
                %grid.lon = LON(ind,:).';
                
                grid.lat = LAT(:,1);
                %THIS IS ALSO TERRIBLE! CAN USE SPLINE BUT DEFINITELY NOT
                %LINEAR INTERPOLATION. ANYWAY, MUCH BETTER TO JUST COMPUTE IT.
                %grid.lon = interp1(LAT(:,1),LON(:,:),lat0,'spline','extrap'); grid.lon = grid.lon';
                grid.lon = obj.lonatorigin(lat0,lon0);
            else
                % this is used for the Matlab Mapping Toolbox options
%                 for i = 1:obj.nx+1
%                     grid.lat(i,1) = mean(LAT(:,i));
%                 end
%                 for j = 1:obj.ny+1                    
%                     grid.lon(j,1) = mean(LON(j,:));
%                 end
                grid.lat = LAT;
                grid.lon = LON;
            end
            
            %sprintf('%.16f\n',(grid.lon - lonatlat0))
            grid.depth = [0; cumsum(obj.dz)] + obj.origin(3);
            grid.units = obj.units;
            
            grid.dlat = diff(grid.lat);
            grid.dlon = diff(grid.lon);
            grid.dz = obj.dz;
            
            grid.nlat = length(grid.dlat);
            grid.nlon = length(grid.dlon);
            grid.nzEarth = length(grid.dz);
                        
            % fix air layers
            grid.zAir = obj.zAir;
            grid.nzAir = length(grid.zAir) - 1;
            
            % save map projection
            grid.mstruct = mstruct;
            
            % convert back to km
            if contains(inUnits,'km')
                grid = grid.m2km;
            end
            
            % padding (use xpadding for both)
            grid = grid.update_padding(obj.xpadding,obj.ypadding);
                
          
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = newspacing(obj,deltax,deltay)
            
            % obj = newspacing(obj,deltax,deltay)
            %
            % take the original grid and modify it so that the new x/y
            % spacing is now deltax & deltay
            %
            % this clears the padding which needs to be created again
            
            if nargin<3
                deltay = deltax;
            end
            
            if isempty(obj.xpadding)
                nxpad = 0;
            else
                nxpad = obj.xpadding;
            end
            
            if isempty(obj.ypadding)
                nypad = 0;
            else
                nypad = obj.ypadding;
            end
            
            xpad1 = obj.dx(1:nxpad);
            xpad2 = obj.dx(end-nxpad+1:end);
            ypad1 = obj.dy(1:nypad);
            ypad2 = obj.dy(end-nypad+1:end);
            
            oldx = cumsum(obj.dx(nxpad+1:end-nxpad)) + obj.origin(1);
            oldy = cumsum(obj.dy(nypad+1:end-nypad)) + obj.origin(2);
            
            newx = oldx(1):deltax:oldx(end);
            newy = oldy(1):deltay:oldy(end);
            
            %obj.dx = [xpad1 diff(newx) xpad2];
            obj.dx = diff(newx);
            obj.nx = length(obj.dx);
            
            %obj.dy = [ypad1 diff(newy) ypad2];
            obj.dy = diff(newy);
            obj.ny = length(obj.dy);

        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = pad(obj,padding,nPad,increaseFactor)
            
            % obj = pad(obj,padding,nPad,increaseFactor)
            %
            % pad the Earth grid by nPad in the directions specified
            % by padding; optional increaseFactor defines the logarithmic
            % increase in cell sizes (except in the air - that's fixed)
            % padding = 'EWSNUD' or any subset thereof.
            %
            % use the trick obj.pad('NSEW',0) to set nPad=0 when it is
            % necessary to overwrite it, such as for obj.select.
            %
            % padding is a string that contains letters W,E,N,S,U,D
            
            if nargin < 3
                nPad = 5;
            elseif nPad <= 0
                warning('no padding required but setting obj.nPad = 0');
                obj.xpadding = 0;
                obj.ypadding = 0;
                return
            else
                obj.xpadding = nPad;
                obj.ypadding = nPad;
            end
            
            if nargin < 4
                increaseFactor = 1.4;
            end
            
            pad = cumsum(exp(log(increaseFactor)*(1:nPad)'));

            if contains(padding,'W')
                padWest = obj.dy(1)*pad; 
                obj.dy = [flipud(padWest); obj.dy];
            end
            
            if contains(padding,'E')
                padEast = obj.dy(end)*pad; 
                obj.dy = [obj.dy; padEast];
            end
            
            if contains(padding,'S')
                padSouth = obj.dx(1)*pad; 
                obj.dx = [flipud(padSouth); obj.dx];
            end
            
            if contains(padding,'N')
                padNorth = obj.dx(end)*pad ; 
                obj.dx = [obj.dx; padNorth];
            end
            
            if contains(padding,'U') % default [32 16 8 4 2]
                padUp = 2.^(obj.nzAir:-1:1);  % km
                obj.zAir = -cumsum(padUp)'; % cell nodes not cell centers?
            end
            
            if contains(padding,'D')
                padDown = obj.dz(end)*increaseFactor.^(1:nPad);
                obj.dz = [obj.dz; padDown];
            end
                        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,ii,jj,kk] = trim(obj,nTrim,direction)
            
            % [obj,ii,jj,kk] = trim(obj,nTrim,direction)
            %
            % use to trim the grid by nTrim cells in a particular direction 
            % direction is a string that contains letters W,E,N,S,U,D
            %
            % intended use is to trim artifacts at the grid boundaries
            % caused by conversions between spherical and cartesian coords,
            % or to remove padding. By default, removes all padding.
            %
            % optional index output ii,jj,kk defines the indices that are
            % NOT trimmed - used to trim the model values.
            
            if nargin < 2
                nTrim = obj.nPad;
            end
            
            if nargin < 3
                direction = obj.padding;
            end
 
            if isfield(obj,'padding')
            if strcmp(direction,obj.padding)
                if obj.nPad > 0
                    obj.nPad = max(0,obj.nPad-nTrim);
                end
            end
            end
            
            ii = 1:obj.nx;
            jj = 1:obj.ny;
            kk = 1:obj.nzEarth;
            x0 = obj.origin(1);
            y0 = obj.origin(2);
            if length(obj.origin) > 2
                z0 = obj.origin(3);
            else
                z0 = 0;
            end

            if strfind(direction,'W')
                jj = jj(1+nTrim:end);
                y0 = y0 + sum(obj.dy(1:nTrim));
            end
            
            if strfind(direction,'E')
                jj = jj(1:end-nTrim);
            end
            
            if strfind(direction,'S')
                ii = ii(1+nTrim:end);
                x0 = x0 + sum(obj.dx(1:nTrim));
            end
            
            if strfind(direction,'N')
                ii = ii(1:end-nTrim);
            end
            
            if strfind(direction,'U')
                kk = kk(1+nTrim:end);
                z0 = z0 + sum(obj.dz(1:nTrim));
            end
            
            if strfind(direction,'D')
                kk = kk(1:end-nTrim);
            end
            
            obj.dx = obj.dx(ii);
            obj.dy = obj.dy(jj);
            obj.dz = obj.dz(kk);
            obj.origin = [x0 y0 z0];
            
            % ... these are dependent properties
            %obj.nx = length(ii);
            %obj.ny = length(jj);
            %obj.nzEarth = length(kk);
                                    
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [cov] = mask(obj,varargin)
            
            % [cov] = mask(obj,varargin)
            %
            % initialized the covariance mask for an XY grid
            % example usage:
            % [COV] = mask(obj,'layer',[15 100],2)
            % will set up a covariance file with index 2
            % for depths below 15 km and above 100 km
            % [COV] = mask(obj,'block',[X1 X2 Y1 Y2 Z1 Z2],4)
            % will define a block with these coordinates (from center)
            % using integer 4
            % [COV] = mask(obj,'ring',[5 15],[0 100],[70 100])
            % will define a ring between depths 5 and 15 km, with center
            % [0 100] and the radii 70 and 100 km, respectively,
            % using integer 5
            % [COV] = mask(obj,'cylinder',[5 15],[0 100],50)
            % will define a cylinder between depths 5 and 15 km, with
            % center [0 100] and radius 50 km
            % using integer 6
            % [COV] = mask(obj,'checkerboard',[15 100],10,[2 3])
            % will define a checkerboard with ~ 10 km cubed blocks
            % (grid permitting, and excluding padding) with integers 2 & 3
            % 
            % it's ok to use vertical arrays on input
            % optionally, pass the mask integer value at the end;
            % then, optionally, pass input covariance as final argument
            %
            % the covariance mask is 3D, try plot(obj,cov,0) 
            % OR
            % for k = 1:nz+1
            %     figure; pcolor(y,x,cov(:,:,k)); colorbar
            % end

            cov = ones(obj.Nx,obj.Ny,obj.NzEarth);
            ctrz = obj.zctr;
            ctrx = obj.xctr;
            ctry = obj.yctr;
            padding = sum(obj.dx(1:obj.xpadding));
            
            if nargin < 2
                return
            elseif nargin == 2
                error(['Please specify a value for argument ',kind])
            end
            
            if nargin > 2
                n = length(varargin);
                option = lower(varargin{1});
                switch option
                    case 'layer'
                        depths = varargin{2};
                        nlayers = size(depths,1);
                        mask = (1:nlayers)+1;
                        if n>2
                            mask = varargin{3};
                        end
                        if length(mask) ~= nlayers
                            error('The number of masks needs to match the layer array')
                        end
                        if n>3
                            cov = varargin{4};
                        end
                    case 'block'
                        bounds = varargin{2};
                        X1 = bounds(:,1) + padding;
                        X2 = bounds(:,2) + padding;
                        Y1 = bounds(:,3) + padding;
                        Y2 = bounds(:,4) + padding;
                        Z1 = bounds(:,5) + padding;
                        Z2 = bounds(:,6) + padding;
                        mask = 2;
                        if n>2
                            mask = varargin{3};
                        end
                        if length(mask) ~= size(bounds,1)
                            error('Need 1 mask for each block')
                        end
                        if n>3
                            cov = varargin{4};
                        end
                    case 'ring'
                        depths = varargin{2};
                        center = varargin{3} + padding;
                        radius = varargin{4};
                        mask = 5;
                        if n>4
                            mask = varargin{5};
                        end
                        if n>5
                            cov = varargin{6};
                        end
                    case 'cylinder'
                        depths = varargin{2};
                        center = varargin{3} + padding;
                        radius = varargin{4};
                        mask = 5;
                        if n>4
                            mask = varargin{5};
                        end
                        if n>5
                            cov = varargin{6};
                        end
                    case 'checkerboard'
                        depths = varargin{2};
                        nlayers = size(depths,1);
                        blocks = varargin{3};
                        mask = [1 2];
                        if n>3
                            mask = varargin{4};
                        end
                        if length(blocks) ~= nlayers
                            error('Need 1 block size for each layer of checkerboard')
                        elseif size(mask,1) ~= nlayers
                            error('Need 2 masks for each layer of checkerboard')
                        end
                        if n>4
                            cov = varargin{5};
                        end
                    otherwise
                        error('Optional argument not defined')
                end                
            end
            
            switch (option)
                case 'layer'
                    for k=1:nlayers
                        ii = ctrz > depths(k,1) & ctrz < depths(k,2);
                        cov(:,:,ii) = mask(k);
                    end
                case 'block'
                    for k=1:length(X1)
                        ii = ctrx > X1(k) & ctrx < X2(k);
                        ij = ctry > Y1(k) & ctry < Y2(k);
                        ik = ctrz > Z1(k) & ctrz < Z2(k);
                        cov(ii,ij,ik) = mask;
                    end
                case 'ring'
                    for k=1:obj.nzEarth
                        if ctrz(k) < depths(1) || ctrz(k) > depths(2)
                            continue
                        end
                        disp(['Creating a ring in grid layer: ' num2str(k)]);
                        for i=1:obj.nx
                            for j=1:obj.ny
                                gridcell(1) = ctrx(i);
                                gridcell(2) = ctry(j);
                                gridcell(3) = ctrz(k);
                                D = sqrt((gridcell(1)-center(1))^2+(gridcell(2)-center(2))^2);
                                %D = pdist2(center,gridcell);
                                if D >= radius(1) && D <= radius(2)
                                    cov(i,j,k) = mask;
                                end
                            end
                        end
                    end
                case 'cylinder'
                    for k=1:obj.nzEarth
                        if ctrz(k) < depths(1) || ctrz(k) > depths(2)
                            continue
                        end
                        disp(['Creating a cylinder in grid layer: ' num2str(k)]);
                        for i=1:obj.nx
                            for j=1:obj.ny
                                gridcell(1) = ctrx(i);
                                gridcell(2) = ctry(j);
                                gridcell(3) = ctrz(k);
                                D = sqrt((gridcell(1)-center(1))^2+(gridcell(2)-center(2))^2);
                                %D = pdist2(center,gridcell);
                                if D <= radius(1)
                                    cov(i,j,k) = mask;
                                end
                            end
                        end
                    end
                case 'checkerboard'
                    for n=1:nlayers
                        Z1 = depths(n,1);
                        Z2 = depths(n,2);
                        nX = floor((sum(obj.dx)-2*padding)/blocks(n));
                        nY = floor((sum(obj.dy)-2*padding)/blocks(n));
                        %nZ = floor((depths(n,2)-depths(n,1))/blocks(n));
                        %ind = 1;
                        ik = ctrz > Z1 & ctrz < Z2;
                        for i = 1:nX+1
                            for j = 1:nY+1
                                X1 = padding + (i-1)*blocks(n);
                                X2 = padding + i*blocks(n);
                                Y1 = padding + (j-1)*blocks(n);
                                Y2 = padding + j*blocks(n);
                                ii = ctrx > X1 & ctrx < X2;
                                ij = ctry > Y1 & ctry < Y2;
                                if mod(abs(i-j),2)
                                    cov(ii,ij,ik) = mask(n,1);
                                else
                                    cov(ii,ij,ik) = mask(n,2);
                                end
                                %bounds(ind,:) = [X1 X2 Y1 Y2 Z1 Z2];
                                %ind = ind+1;
                            end
                        end
                     end
                otherwise
                    % nothing
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = uiplot(obj,value,padding,cblabel)
            
            if nargin > 2
                obj.xpadding = padding;
                obj.ypadding = padding;
            end
            options.padding = obj.xpadding;
            if nargin > 3
                options.cblabel = cblabel;
            else
                options.cblabel = '';
            end
            options.slice = 'Z';
            options.Np = 1;
            [Nx,Ny,Nz] = size(value);
            options.iXlim(1) = 1;
            options.iXlim(2) = Nx+1;
            options.iYlim(1) = 1;
            options.iYlim(2) = Ny+1;
            options.iZlim(1) = 1;
            options.iZlim(2) = Nz+1;
            CondPlotSet(value,obj,options);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = m2km(obj)
           
            if ~contains(obj.units,'km')
                obj.dx = obj.dx / 1e3;
                obj.dy = obj.dy / 1e3;
                obj.dz = obj.dz / 1e3;
                obj.zAir = obj.zAir / 1e3;
                obj.origin = obj.origin / 1e3;
                obj.units = 'km';
            end
                    
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = km2m(obj)
           
            if contains(obj.units,'km')
                obj.dx = 1e3 * obj.dx;
                obj.dy = 1e3 * obj.dy;
                obj.dz = 1e3 * obj.dz;
                obj.zAir = 1e3 * obj.zAir;
                obj.origin = 1e3 * obj.origin;
                obj.units = 'm';
            end
                    
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [LON] = lonatorigin(obj,lat0,lon0)
            xy(1,:) = zeros(1,obj.ny+1);
            xy(2,:) = obj.y;
            ll = latlontools.xy2ll(xy,lat0,lon0);
            LON = ll(2,:)';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [LAT,LON] = latlon(obj,lat0,lon0)
            % NOTE: the LON is WRONG when computing the limits for an
            % xymodel but I cannot debug this now
            LAT = (obj.x/obj.kmPerDeg(lat0))*ones(1,obj.ny+1) + lat0;
            LON = (ones(obj.nx+1,1)*obj.y'/obj.kmPerDeg(lat0))./cos(LAT*pi/180) + lon0;
%           or, alternatively...
%             xxyy(1,:) = repmat(obj.x',1,obj.ny+1);
%             tmp = repmat(obj.y,1,obj.nx+1);
%             xxyy(2,:) = reshape(tmp',1,(obj.nx+1)*(obj.ny+1));
%             llll = latlontools.xy2ll(xxyy,lat0,lon0);
%             LAT = reshape(llll(1,:),obj.nx+1,obj.ny+1);
%             LON = reshape(llll(2,:),obj.nx+1,obj.ny+1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [LAT,LON] = latlonctr(obj,lat0,lon0)
            LAT = (obj.xctr/obj.kmPerDeg(lat0))*ones(1,obj.ny) + lat0;
            LON = (ones(obj.nx,1)*obj.yctr'/obj.kmPerDeg(lat0))./cos(LAT*pi/180) + lon0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [LAT,LON] = latlonbounds(obj,lat0,lon0)
            LAT = ([min(obj.x) max(obj.x)]/obj.kmPerDeg(lat0))*ones(1,obj.ny) + lat0;
            LON = (ones(obj.nx,1)*[min(obj.y) max(obj.y)]'/obj.kmPerDeg(lat0))./cos(LAT*pi/180) + lon0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = x(obj)
            x = [0 ; cumsum(obj.dx)]+obj.origin(1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = y(obj)
            y = [0 ; cumsum(obj.dy)]+obj.origin(2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function z = z(obj)
            if ~isempty(obj.zAir)
                z = [- sort(obj.zAir,1,'descend'); cumsum(obj.dz)]+obj.origin(3);
            else
                z = [0 ; cumsum(obj.dz)]+obj.origin(3);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = xctr(obj)
            x = (cumsum(obj.dx) - obj.dx/2) + obj.origin(1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = yctr(obj)
            y = (cumsum(obj.dy) - obj.dy/2) + obj.origin(2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function z = zctr(obj)
            tmp = - sort(obj.zAir,1,'descend');
            if ~isempty(obj.zAir)
                z = (tmp(1:end-1) + diff(tmp)/2);
            else
                z = [];
            end            
            z = [z; (cumsum(obj.dz) - obj.dz/2)] + obj.origin(3);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = get.Nx(obj)
            res = obj.nx;
        end
        function res = get.nx(obj)
            res = length(obj.dx);
        end
        function obj = set.nx(obj,~) %#ok<MANU>
            error('The value of nx is computed from dx, can''t set it explicitly');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = get.Ny(obj)
            res = obj.ny;
        end
        function res = get.ny(obj)
            res = length(obj.dy);
        end
        function obj = set.ny(obj,~) %#ok<MANU>
            error('The value of ny is computed from dy, can''t set it explicitly');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = get.NzEarth(obj)
            res = obj.nzEarth;
        end
        function res = get.nzEarth(obj)
            res = length(obj.dz);
        end
        function obj = set.nzEarth(obj,~) %#ok<MANU>
            error('The value of nzEarth is computed from dz, can''t set it explicitly');
        end
    end
    
    methods(Static)
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dz] = logz(top,bottom,nz,logdz)
            
            % [dz] = logz(top,bottom,nz,logdz)
            % 
            % generated a logarithmically spaced set of grid layers
            % logdz makes depth distribution logarithmic;
            % (top - bottom) = sum_{k=0}^nz (startdz * logdz^k)
            %                = startdz * (logdz^(nz+1) - 1)/(logdz - 1)
            % therefore startdz = (top - bottom)*(logdz - 1)/(logdz^(nz+1) - 1)
            % and n = 0:nz-1; z = cumsum(startdz*(logdz.^n)).          
            
            if (logdz == 1) && (nz>0)
                % logdz = 1;
                dz(1:nz) = (bottom-top)/nz;
                disp(['Spacing: ' num2str(dz(1)) ' km']);
            elseif nz>0
                startdz = (bottom-top)*(logdz - 1)/(logdz^nz - 1);
                n = 0:nz-1; dz = startdz*(logdz.^n);
                disp(['Spacing: starts with ' num2str(dz(1)) ' km']);
            else
                dz = [];
                disp('Number of layers not a positive integer; no layers created');
            end
            
            dz = dz';
                                         
        end
    end
end