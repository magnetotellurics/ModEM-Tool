classdef llgrid < latlontools
    %   class to define a regular lat/lon grid; uses ETOPO2 database
    %   to create an array that masks air & sea water - used for
    %   covariance and setting up a prior model for inversion
    % 
    %   all calculations including merging models, interpolation, etc
    %   are significantly streamlined if the grid lat/lon/depth coordinates
    %   are always defined at CELL NODES; locations of cell centers can
    %   always be straighforwardly computed using dlat/dlon/dz
    %   [e.g. ctrlat = lat(1:end-1) + dlat/2 ].
    %   this is the convention we are going to adopt throughout.
    %   note that models can still be defined at either CELLS or NODES.
    %   can go in the backwards direction (cell centers to nodes) using
    %   latlontools.delta, but this is ambiguous, not for routine use.
    %   
    properties   
        lat % always defined at cell nodes
        lon
        depth  
        dlat % define cell sizes
        dlon
        dz
        units = 'km';
        zAir = []; % +ve height above ground, length nzAir+1
        nzAir = 12;
        nzCrust = 0; % ignore this for MT; used for global
        nzEarth = 0;
        nlat = 0; % number of latitude cells 
        nlon = 0; % number of longitude cells
        limits;
        mstruct;
       
    end
    properties (Constant)
        % covariance mask as defined in ModEM
        AIR = 0;
        EARTH = 1;
        OCEAN = 9;
    end
    properties (SetAccess = protected)
        % padding parameters
        increaseFactor = 1.4;
        nPad = 0;
        padding = '';
        logdz = 1;
    end
    properties (SetAccess = public)
        % 2D mesh of LAT/LON values; regular if computed from lat & lon;
        % irregular if obtained by interpolation from an x/y grid
        %LAT
        %LON
    end
    properties (SetAccess = protected)
        % a non-regular mesh of X/Y values computed using a lat0/lon0
        % origin in geographic coords and used for interpolation to a
        % regular x/y grid (class xygrid object)
        %X
        %Y
        %Z
        % indices for air / water / Earth regions obtained from ETOPO2
        %iAir
        %iOcean
        %iEarth
    end
    methods
        function [obj] = llgrid(lims,top,bottom,nlat,nlon,nz,varargin)
            % [obj] = llgrid(lims,top,bottom,nlat,nlon,nz,varargin)
            %   class constructor; lims must contain lon/lat min/max
            if nargin<1
                % allow an empty constuctor
                obj.limits.latmin = -90;
                obj.limits.latmax = 90;
                obj.limits.lonmin = 0;
                obj.limits.lonmax = 360;
                obj.limits.depthmin = 0;
                obj.limits.depthmax = 3500; % km
                return
            elseif ~isstruct(lims)
                error('First argument must be a structure that contains min and max values in lon & lat')
            end
            obj.limits.latmin = lims.latmin;
            obj.limits.latmax = lims.latmax;
            obj.limits.lonmin = obj.lon360(lims.lonmin);
            obj.limits.lonmax = obj.lon360(lims.lonmax);
            if isfield(lims,'depthmin')
                obj.limits.depthmin = lims.depthmin;
            end
            if isfield(lims,'depthmax')
                obj.limits.depthmax = lims.depthmax;
            end
            if nargin == 1
                % define the grid with limits and exit
                return
            end
            obj.limits.depthmin = top;
            obj.limits.depthmax = bottom;
            if ~isempty(varargin)
                obj.logdz = varargin{1};
            end
            if length(varargin)>1
                % no other arguments at present
            end
            if nargin < 6
                nz = 0;
            end
            if nargin >= 5
                [gridlat,gridlon,gridz] = obj.mesh(lims,top,bottom,nlat,nlon,nz,obj.logdz);
                obj.dlat = diff(gridlat);
                obj.dlon = diff(gridlon);
                obj.dz = diff(gridz);
                obj.lat = gridlat;
                obj.lon = gridlon;
                if length(obj.dz)>1
                    obj.depth = gridz;
                end
                obj.nlat = length(obj.dlat);
                obj.nlon = length(obj.dlon);
                obj.nzEarth = length(obj.dz);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function grid = xygrid(obj,mstruct)
            
            % convert lat/lon to a regular grid in km.
            % set the origin to match that of the data set;
            % we want to use this origin for all grid distance computations
            % to ensure that these match the distances for the data set.
            % e.g., for the Yellowstone data,
            % [lat0, lon0] = [42.016 -112.524 top]; % smaller area
            % [lat0, lon0] = [42.016 -116.477 top];
            %
            % generalized based on Ph.D. work of Han Qi to use the Matlab
            % Plotting Toolbox, which uses the mstruct structure to set up
            % the coordinate projection. Had to completely rewrite it 
            % to use the true grid origin for conversion.
            %
            % This will not work without the toolbox, but some basic methods 
            % will still be available.
            %
            % At a minimum, need to specify:
            %   mstruct.origin = [lat0,lon0,0]
            %   mstruct.mapprojection = 'eqdcylin' 
            %
            % Some projection options: 
            % eqdcylin, eqacylin, lambertstd, eqaazim, utm
            % type 'maps' for all projections
            % To trigger the homemade projection option, use 
            %   mstruct.mapprojection = 'latlon2xy'
            % 
            % Set up lambertstd like this:
            % mstruct = defaultm('lambertstd');
            % mstruct.origin = [lat0,lon0,0];
            % mstruct.mapparallels = [lat0-3,lat0+3]; % degrees vary
            %
            % Set up UTM like this:
            % z1 = utmzone([latm,lonm]); % e.g., z1 = '15N'
            % ellipsoid = wgs84Ellipsoid;
            % mstruct = defaultm('utm');
            % mstruct.zone = z1;
            % mstruct.geoid = ellipsoid;
            %
            % If using UTM, note that the binding Northing and Easting are
            % added to your grid origin on output.
            % If using UTM, make sure that you add the binding Northing and
            % Easting to the grid origin before applying this routine.
                        
            center = [mean(obj.lat) mean(obj.lon)];
            
            % by default, use the homemade projection
            if nargin < 2
                mstruct.mapprojection = 'latlon2xy';
            end
            
            % use center for origin by default
            if ~isfield(mstruct,'origin')
                mstruct.origin(1) = center(1);
                mstruct.origin(2) = center(2);
            end
            
            % by definition, the grid always defines cell edges/nodes
            % so the number of these is always one greater than the number
            % of cells. If it is not so, issue an error to fix this before
            % we get started.
            if length(obj.lon) == length(obj.dlon)+1 && length(obj.lat) == length(obj.dlat)+1
                gridlon = obj.lon;
                gridlat = obj.lat;
            else
                error('Grid size mismatch: make sure that lon & lat define cell edges, not cell centers!');
            end
            
            
            % origin defines the output origin in km, too
            lat0 = mstruct.origin(1);
            lon0 = mstruct.origin(2);
            olat = repmat(lat0,length(gridlon),1);
            olon = repmat(lon0,length(gridlat),1);
            
            % convert the grid to meters
            inUnits = obj.units;
            obj = obj.km2m;
            
            if strcmp(mstruct.mapprojection,'latlon2xy')
                
                y = obj.ll2xy([olat gridlon].',lat0,lon0); dY = diff(y(2,:).');
                x = obj.ll2xy([gridlat olon].',lat0,lon0); dX = diff(x(1,:).');
                
                origin = - [sum(dX)/2 sum(dY)/2 0]'; % origin is grid center by default
                km = obj.ll2xy(center',lat0,lon0); % distance of center from origin
                km(3) = 0;
                origin = origin + km; % compute new origin in km
                
            else
                try
                    mstruct = defaultm(mstruct);
                catch
                    try
                        mstruct = defaultm(mstruct.mapprojection);
                        mstruct.origin = [lat0,lon0,0];
                        mstruct = defaultm(mstruct);
                    catch
                        error(['Unknown map projection: ',mstruct.mapprojection]);
                    end
                end
                
                switch strtrim(mstruct.mapprojection)
                    
                    case 'utm'
                        
                        % convert grid nodes to Cartesian
                        [E,~] = mfwdtran(mstruct,olat,gridlon); %E = rad2km(E);
                        [~,N] = mfwdtran(mstruct,gridlat,olon); %N = rad2km(N);
                        
                        % create grid distances
                        dY = diff(E);
                        dX = diff(N);
                        
                        % convert grid origin to Cartesian - BUT output
                        % origin is ALWAYS the lower left corner, as
                        % currently set up
                        [oY,oX] = mfwdtran(mstruct,lat0,lon0);
                        
                        % compute lower left corner in Cartesian
                        minY = min(E); minX = min(N);
                        
                        % compute distance to origin
                        origin = [minX,minY,0]'; % - [oX,oY,0]';
                        
                    otherwise
                        %case {'eqdcylin','eqacylin','lambertstd','eqaazim'}
                        
                        % convert grid nodes to Cartesian
                        [E,~] = projfwd(mstruct,olat,gridlon); %E = rad2km(E);
                        [~,N] = projfwd(mstruct,gridlat,olon); %N = rad2km(N);
                        
                        % create grid distances
                        dY = diff(E);
                        dX = diff(N);
                        
                        % convert grid origin to Cartesian - BUT output
                        % origin is ALWAYS the lower left corner, as
                        % currently set up
                        [oY,oX] = projfwd(mstruct,lat0,lon0);
                        
                        % compute lower left corner in Cartesian
                        minY = min(E); minX = min(N);
                        
                        % compute distance to origin
                        origin = [minX,minY,0]'; % - [oX,oY,0]';
                        
                end
                
                % now convert to meters
                %dY = 1000.0*dY;
                %dX = 1000.0*dX;
            end
                        
               
            % at this point, all distances are in meters
            dZ = obj.dz;
            grid = xygrid(dX,dY,dZ,origin,0,'m');
            
            % fix air layers
            grid.zAir = obj.zAir;
            grid.nzAir = obj.nzAir;
            
            % AND convert all distances back to km
            if contains(inUnits,'km')
                grid = grid.m2km;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function z = z(obj)
            if ~isempty(obj.zAir)
                z = [- sort(obj.zAir,1,'descend'); cumsum(obj.dz)];
            else
                z = [0 ; cumsum(obj.dz)];
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function z = zctr(obj)
            tmp = - sort(obj.zAir,1,'descend');
            if ~isempty(obj.zAir)
                dzAir = diff(tmp);
            else
                dzAir = [];
            end            
            z = [(tmp(1:end-1) + dzAir/2); (cumsum(obj.dz) - obj.dz/2)];
        end
             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = m2km(obj)
           
            if ~contains(obj.units,'km')
                obj.depth = obj.depth / 1e3;
                obj.dz = obj.dz / 1e3;
                obj.zAir = obj.zAir / 1e3;
                obj.units = 'km';
            end
                    
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = km2m(obj)
           
            if contains(obj.units,'km')
                obj.depth = 1e3 * obj.depth;
                obj.dz = 1e3 * obj.dz;
                obj.zAir = 1e3 * obj.zAir;
                obj.units = 'm';
            end
                    
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = insert(obj,newobj,direction)
            
            % obj = insert(obj,newobj,direction)
            %
            % inserts e.g., a finer grid into a large & coarse grid
            % might fail for global grids - no error checking yet
            %
            % if present, direction = 'depth', 'latitude' or 'longitude'
            
            if nargin < 3
                direction = 'latitude,longitude,depth';
            end
            
            if contains(direction,'lat')
                i = obj.lat<newobj.limits.latmin | obj.lat>newobj.limits.latmax;
                obj.lat = obj.lat(i);
                obj.lat = sort([obj.lat; newobj.lat]);
                obj.dlat = diff(obj.lat);
                obj.nlat = length(obj.dlat);
            end
            
            if contains(direction,'lon')
                j = obj.lon<newobj.limits.lonmin | obj.lon>newobj.limits.lonmax;
                obj.lon = obj.lon(j);
                obj.lon = sort([obj.lon; newobj.lon]);
                obj.dlon = diff(obj.lon);
                obj.nlon = length(obj.dlon);
            end
            
            if contains(direction,'depth')
                k = obj.depth<newobj.limits.depthmin | obj.depth>newobj.limits.depthmax;
                obj.depth = obj.depth(k);
                obj.depth = sort([obj.depth; newobj.depth]);
                obj.dz = diff(obj.depth);
                obj.nzEarth = length(obj.dz);
            end
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = plot(obj,value,depth,sitelat,sitelon)
            
            % plot(obj,value,depth,sitelat,sitelon)
            %
            % a simple pcolor plot of a value of the grid
            % at a specified depth in km
            
            if size(value,3)>1 && nargin > 2
                k = find(obj.depth>=depth, 1, 'first');
            else
                k = 1;
            end
            ctrlon = obj.lon(1:end-1) + obj.dlon/2;
            ctrlat = obj.lat(1:end-1) + obj.dlat/2;
            if size(value,1) ~= length(ctrlat) ...
                || size(value,2) ~= length(ctrlon)
                error('The value to plot should have dimensions (nlat,nlon,:)');
            end
            if size(value,3)==1
                disp('Plotting the 2D map');
            elseif obj.depth(k)==0
                disp('Plotting the value at Earth''s surface');
            elseif obj.depth(k)<0
                disp(['Plotting the value at elevation ' num2str(-obj.depth(k)) ' km']);
            else
                disp(['Plotting the value at depth ' num2str(obj.depth(k)) ' km']);
            end
            figure; 
            pcolor(ctrlon,ctrlat,squeeze(value(:,:,k)));
            colorbar; hold on;
            if nargin > 4
                plot(sitelon,sitelat,'k^','MarkerSize',10,'LineWidth',2);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [cov,elev] = mask(obj)
            
            % [cov,elev] = mask(obj)
            %
            % computes the covariance mask for air, ocean and Earth
            % for an initialized lat/lon grid
            %
            % plot the ETOPO2 elevation using
            % pcolor(mlon,mlat,melev); shading interp; colorbar; caxis([-1000 1000])
            %
            % the covariance mask is 3D, try plot(obj,cov,0) 
            % OR
            % for k = 1:nz+1
            %     figure; pcolor(y,x,cov(:,:,k)); colorbar
            % end

            minlat = obj.limits.latmin;
            maxlat = obj.limits.latmax;
            minlon = min(llgrid.lon180(obj.lon));
            maxlon = max(llgrid.lon180(obj.lon));
            
            [melev mlon mlat]=m_etopo2([minlon maxlon minlat maxlat]);
            
            % interpolate for plotting
            %[y,x] = meshgrid(obj.lon,obj.lat);
            %elev = interp2(mlon,mlat,melev,y,x);
            %pcolor(y,x,elev); colorbar;
            
            % interpolate for elevation comparisons at cell centers
            ctrlon = obj.lon(1:end-1) + obj.dlon/2;
            ctrlat = obj.lat(1:end-1) + obj.dlat/2;
            [y,x] = meshgrid(ctrlon,ctrlat);
            elev = interp2(mlon,mlat,melev,y,x);
            
            % compute grid cell centers (+ve up), all in meters
            if strcmp(obj.units,'km')
                ctrz = - 1000*(obj.depth(1:end-1) + obj.dz/2);
            else
                ctrz = - (obj.depth(1:end-1) + obj.dz/2);                
            end
            
            cov = zeros(length(ctrlat),length(ctrlon),length(ctrz));
            v   = zeros(length(ctrlat),length(ctrlon));
            
            for k = 1:length(ctrz)
                if ctrz(k)>0
                    iAir = elev<ctrz(k); v(iAir) = llgrid.AIR;
                    iEarth = elev>=ctrz(k); v(iEarth) = llgrid.EARTH;
                else
                    iOcean = elev<ctrz(k); v(iOcean) = llgrid.OCEAN;
                    iEarth = elev>=ctrz(k); v(iEarth) = llgrid.EARTH;
                end
                cov(:,:,k) = v;
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [idx] = index(obj,mindepth,maxdepth)
            
            % [idx] = index(obj,mindepth,maxdepth)
            %
            % the output idx is a nlat x nlon x nzEarth sized matrix with 
            % 1's corresponding to cells whose centers are vertically below 
            % mindepth and above maxdepth, and 0's otherwise
            %
            % mindepth and maxdepth could be nlat x nlon matrices, or
            % scalars
            %
            % example usage: index all cells that are below the
            % bathymetry level and correspond to the depths of the sediments
            % [COV,ELEV] = grid.mask;
            % idx = grid.index(ELEV,ELEV-SED); 
            % COV(idx) = 3;
                        
            if isscalar(mindepth)
                top(1:obj.nlon,1:obj.nlat) = mindepth;
            elseif size(mindepth,1) == obj.nlat && size(mindepth,2) == obj.nlon
                top = mindepth;
            else
                error('The size of mindepth must be match the grid (nlat x nlon)');
            end
            
            if isscalar(maxdepth)
                bottom(1:obj.nlon,1:obj.nlat) = maxdepth;
            elseif size(maxdepth,1) == obj.nlat && size(maxdepth,2) == obj.nlon
                bottom = maxdepth;
            else
                error('The size of maxdepth must be match the grid (nlat x nlon)');
            end
            
            % compute grid cell centers (+ve up), all in meters
            if strcmp(obj.units,'km')
                ctrz = - 1000*(obj.depth(1:end-1) + obj.dz/2);
            else
                ctrz = - (obj.depth(1:end-1) + obj.dz/2);                
            end
           
            idx = zeros(obj.nlat,obj.nlon,obj.nzEarth);
            
            for j = 1:obj.nlon
                for i = 1:obj.nlat
                    idx(i,j,:) = (ctrz < top(i,j)) & (ctrz > bottom(i,j));
                end
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function newobj = select(obj,lims)
            %
            % newobj = select(obj,lims)
            %
            % select a 3D area based on lims which are either a string
            % (from the options in modelplot) or a 4 or 6 value vector
            % (lonmin lonmax latmin latmax [depthmin depthmax])
            
            if nargin < 2
                error('To select an llgrid based on new limits, need to specify lims');
            elseif ischar(lims)
                region = lims;
                lims = latlontools.getLimits(region);
                lims.depthmin = obj.limits.depthmin;
                lims.depthmax = obj.limits.depthmax;
            elseif length(lims) >= 4
                temp = lims; clear lims;
                lims.lonmin = temp(1);
                lims.lonmax = temp(2);
                lims.latmin = temp(3);
                lims.latmax = temp(4);
                if length(temp) == 6
                    lims.depthmin = temp(5);
                    lims.depthmax = temp(6);
                else
                    lims.depthmin = obj.limits.depthmin;
                    lims.depthmax = obj.limits.depthmax;
                end
            end
            k1 = find(obj.depth >= lims.depthmin, 1, 'first');
            k2 = find(obj.depth <= lims.depthmax, 1, 'last');
            i1 = find(obj.lat >= lims.latmin, 1, 'first');
            i2 = find(obj.lat <= lims.latmax, 1, 'last');
            if issorted(llgrid.lon360(obj.lon))
                j1 = find(llgrid.lon360(obj.lon) >= llgrid.lon360(lims.lonmin), 1, 'first');
                j2 = find(llgrid.lon360(obj.lon) <= llgrid.lon360(lims.lonmax), 1, 'last');
            else
                j1 = find(llgrid.lon180(obj.lon) >= llgrid.lon180(lims.lonmin), 1, 'first');
                j2 = find(llgrid.lon180(obj.lon) <= llgrid.lon180(lims.lonmax), 1, 'last');
            end
                        
            newobj = obj;
            newobj.depth = obj.depth(k1:k2);
            newobj.lat = obj.lat(i1:i2);
            newobj.lon = obj.lon(j1:j2);
            newobj.dz = diff(newobj.depth);
            newobj.dlat = diff(newobj.lat);
            newobj.dlon = diff(newobj.lon);
            newobj.nzEarth = length(k1:k2-1);
            newobj.nlat = length(i1:i2-1);
            newobj.nlon = length(j1:j2-1);
            newobj.limits = lims;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = extend(obj,direction,newvalue,cellsize)
            
            % obj = extend(obj,direction,newvalue,cellsize)
            %
            % use to extend the grid in a particular direction 
            % (depth, lon, lat); use obj.pad to add padding
            % by default, extends using the last cell size available.
            % don't use this to expand into the air! use obj.pad instead.
            
            lims = obj.limits;
                        
            if strcmp(direction,'depth')
                nodez = obj.depth;
                if newvalue > lims.depthmax
                    lims.depthmax = newvalue;
                    if nargin < 4; cellsize = obj.dz(end); end
                    extend = nodez(end)+cellsize:cellsize:lims.depthmax;
                    obj.depth = [obj.depth; extend'];
                    obj.limits.depthmax = newvalue;
                elseif newvalue < lims.depthmin
                    lims.depthmin = newvalue;
                    if nargin < 4; cellsize = obj.dz(1); end
                    extend = fliplr(nodez(1)-cellsize:-cellsize:lims.depthmin);
                    obj.depth = [extend'; obj.depth];
                    obj.limits.depthmin = newvalue;
                end
                obj.dz = diff(obj.depth);
                obj.nzEarth = length(obj.dz);
                
            elseif strfind(direction, 'lon')
                nodelon = obj.lon;
                if latlontools.lon360(newvalue) > lims.lonmax
                    lims.lonmax = newvalue;
                    if nargin < 4; cellsize = obj.dlon(end); end
                    extend = nodelon(end)+cellsize:cellsize:lims.lonmax;
                    obj.lon = [obj.lon; extend'];
                    obj.limits.lonmax = latlontools.lon360(newvalue);
               elseif latlontools.lon360(newvalue) < lims.lonmin
                    lims.lonmin = newvalue;
                    if nargin < 4; cellsize = obj.dlon(1); end
                    extend = fliplr(nodelon(1)-cellsize:-cellsize:lims.lonmin);
                    obj.lon = [extend'; obj.lon];
                    obj.limits.lonmin = latlontools.lon360(newvalue);
                    
                end
                obj.dlon = diff(obj.lon);
                obj.nlon = length(obj.dlon);
                
            elseif strfind(direction, 'lat')                    
                nodelat = obj.lat;
                if newvalue > lims.latmax
                    lims.latmax = newvalue;
                    if nargin < 4; cellsize = obj.dlat(end); end
                    extend = nodelat(end)+cellsize:cellsize:lims.latmax;
                    obj.lat = [obj.lat; extend'];
                    obj.limits.latmax = newvalue;
              elseif newvalue < lims.latmin
                    lims.latmin = newvalue;
                    if nargin < 4; cellsize = obj.dlat(1); end
                    extend = fliplr(nodelat(1)-cellsize:-cellsize:lims.latmin);
                    obj.lat = [extend'; obj.lat];
                    obj.limits.latmin = newvalue;
                end
                obj.dlat = diff(obj.lat);
                obj.nlat = length(obj.dlat);
                
            else
                error('Please specify direction: depth, longitude or latitude');
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = addz(obj,newbottom,startdz,logdz)
            
            % obj = addz(obj,newbottom,startdz,logdz)
            %
            % update the number of vertical layers in the grid using a
            % logarithmic approach. In practice this is usually used 
            % reduce the number of vertical layers to something manageable
            
            if nargin < 2
                error('Need to specify new bottom for the grid');
            end
            
            top = obj.limits.depthmax;
            bottom = newbottom;
            
            if nargin < 3
                startdz = obj.dz(end);
            end
            
            if nargin < 4
                logdz = 1.05;
            end
            
            nz = ceil(log((bottom-top)*(logdz - 1)/(startdz) + 1)/log(logdz));
            n = 0:nz-1; newdz = startdz*(logdz.^n);

            obj.dz = [obj.dz; newdz'];
            obj.nzEarth = obj.nzEarth + length(newdz);
            obj.depth = [obj.depth; top+cumsum(newdz')];
            obj.padding = setdiff(obj.padding,'UD'); % remove UD padding
            obj.limits.depthmax = max(obj.depth);

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
            % note that by convention, lat, lon, depth refer to cell
            % centers so need to be careful when padding.
            %
            % padding is a string that contains letters W,E,N,S,U,D
            
            if nargin < 3
                nPad = 5;
            elseif nPad <= 0
                warning('no padding required but setting obj.nPad = 0');
                obj.padding = setdiff(obj.padding,padding)';
                obj.nPad = 0;
                return
            end
            
            if nargin < 4
                increaseFactor = 1.4;
            end
            
            pad = cumsum(exp(log(increaseFactor)*(1:nPad)));

            nodelat = obj.lat;
            nodelon = obj.lon;
            nodez = obj.depth;

            if contains(padding,'W')
                padWest = obj.dlon(1)*pad; % degrees
                nodelon = [nodelon(1)-fliplr(cumsum(padWest))'; nodelon];
                obj.dlon = diff(nodelon);
                obj.lon = nodelon;
            end
            
            if contains(padding,'E')
                padEast = obj.dlon(end)*pad; % degrees
                nodelon = [nodelon; nodelon(end)+cumsum(padEast)'];
                obj.dlon = diff(nodelon);
                obj.lon = nodelon;
            end
            
            if contains(padding,'S')
                padSouth = obj.dlat(1)*pad; % degrees
                nodelat = [nodelat(1)-fliplr(cumsum(padSouth))'; nodelat];
                obj.dlat = diff(nodelat);
                obj.lat = nodelat;
            end
            
            if contains(padding,'N')
                padNorth = obj.dlat(end)*pad ; % degrees
                nodelat = [nodelat; nodelat(end)+cumsum(padNorth)'];
                obj.dlat = diff(nodelat);
                obj.lat = nodelat;
            end
            
            if contains(padding,'U') % default [32 16 8 4 2]
                padUp = 2.^(obj.nzAir:-1:1);  % km
                obj.zAir = nodez(1)-cumsum(padUp)'; % cell nodes not cell centers?
            end
            
            if contains(padding,'D')
                padDown = obj.dz(end)*increaseFactor.^(1:nPad); % km
                nodez = [nodez; nodez(end)+cumsum(padDown)'];
                obj.dz = diff(nodez);
                obj.depth = nodez;
            end
            
            obj.limits.latmin = min(obj.lat);
            obj.limits.latmax = max(obj.lat);
            obj.limits.lonmin = min(obj.lon);
            obj.limits.lonmax = max(obj.lon);
            obj.limits.depthmin = min(obj.depth);
            obj.limits.depthmax = max(obj.depth);
            obj.nlat = length(obj.dlat);
            obj.nlon = length(obj.dlon);
            obj.nzEarth = length(obj.dz);
            obj.padding = union(obj.padding,padding);
            if size(obj.padding,1)>1; obj.padding = obj.padding'; end
            obj.increaseFactor = increaseFactor;
            obj.nPad = nPad;
            
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
            
            if strcmp(direction,obj.padding)
                if obj.nPad > 0
                    obj.nPad = max(0,obj.nPad-nTrim);
                end
            end
            
            ii = 1:obj.nlat;
            jj = 1:obj.nlon;
            kk = 1:obj.nzEarth;

            if strfind(direction,'W')
                jj = jj(1+nTrim:end);
            end
            
            if strfind(direction,'E')
                jj = jj(1:end-nTrim);
            end
            
            if strfind(direction,'S')
                ii = ii(1+nTrim:end);
            end
            
            if strfind(direction,'N')
                ii = ii(1:end-nTrim);
            end
            
            if strfind(direction,'U')
                kk = kk(1+nTrim:end);
            end
            
            if strfind(direction,'D')
                kk = kk(1:end-nTrim);
            end
            
            obj.dlat = obj.dlat(ii);
            obj.dlon = obj.dlon(jj);
            obj.dz = obj.dz(kk);
            obj.lat = obj.lat(ii(1)) + [0; cumsum(obj.dlat)];
            obj.lon = obj.lon(jj(1)) + [0; cumsum(obj.dlon)];
            obj.depth = obj.depth(kk(1)) + [0; cumsum(obj.dz)];
            obj.nlat = length(ii);
            obj.nlon = length(jj);
            obj.nzEarth = length(kk);
            
            % update limits
            obj.limits.latmin = min(obj.lat);
            obj.limits.latmax = max(obj.lat);
            obj.limits.lonmin = min(obj.lon);
            obj.limits.lonmax = max(obj.lon);
            obj.limits.depthmin = min(obj.depth);
            obj.limits.depthmax = max(obj.depth);
                                    
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = newspacing(obj,deltalat,deltalon)
            
            % obj = newspacing(obj,deltalat,deltalon)
            %
            % take the original grid and modify it so that the lat/lon
            % spacing of the main grid is now delta
            %
            % this clears the padding which needs to be created again
            
            if nargin<3
                deltalon = deltalat;
            end
            
            if strfind(obj.padding,'S') 
                newlat = obj.lat(1+obj.nPad):deltalat:obj.lat(end-obj.nPad);
            else
                newlat = obj.lat(1):deltalat:obj.lat(end);
            end
            if strfind(obj.padding,'W')
                newlon = obj.lon(1+obj.nPad):deltalon:obj.lon(end-obj.nPad);
            else
                newlon = obj.lon(1):deltalon:obj.lon(end);
            end
            
            obj.lat = newlat';
            obj.dlat = diff(obj.lat);
            obj.nlat = length(obj.dlat);
            
            obj.lon = newlon';
            obj.dlon = diff(obj.lon);
            obj.nlon = length(obj.dlon);
            
            obj.padding = '';
            obj.nPad = 0;

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = newz(obj,top,bottom,startdz,logdz)
            
            % obj = newz(obj,top,bottom,startdz,logdz)
            %
            % update the number of vertical layers in the grid using a
            % logarithmic approach. In practice this is usually used 
            % reduce the number of vertical layers to something manageable
            
            if nargin < 2
                top = obj.limits.depthmin;
            end
            
            if nargin < 3
                bottom = obj.limits.depthmax;
            end
            
            if nargin < 4
                startdz = obj.dz(1);
            end
            
            if nargin < 5
                logdz = 1.05;
            end
            
            nz = ceil(log((bottom-top)*(logdz - 1)/(startdz) + 1)/log(logdz));
            n = 0:nz-1; newdz = startdz*(logdz.^n);

            obj.dz = newdz';
            obj.nzEarth = length(newdz);
            obj.depth = [0; cumsum(newdz')];
            obj.padding = setdiff(obj.padding,'UD'); % remove vertical padding
            obj.limits.depthmax = max(obj.depth);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make a copy of a handle object, but store as a child object.
        function new = copy(this,new)
            
            if nargin == 1
                % Instantiate new object of the same class.
                new = feval(class(this));
            elseif isempty(strfind(class(new),'grid'))
                warning('Trying to hard copy an llgrid object into an inacceptible class');
            end
            
            meta = metaclass(this);
 
            % Copy all non-hidden properties.
            % p = properties(this);
            for ii = 1:length(meta.Properties)  
                
                % Do not copy Transient or Dependent or Hidden Properties
                if (meta.Properties{ii}.Transient ...
                        || meta.Properties{ii}.Dependent ...
                        || meta.Properties{ii}.Hidden)
                    continue; 
                end
                
                % Do not copy properties that will be read-only in 'new'
                if strcmp(meta.Properties{ii}.SetAccess,'none')
                    continue;
                end
                
                new.(meta.Properties{ii}.Name) = this.(meta.Properties{ii}.Name);

            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update padding information (no actual changes in the grid)
        function obj = update_padding(obj,xpadding,ypadding,zpadding)
            
            if nargin < 4
                zpadding = 0;
            end
            
            % padding (use xpadding)
            if xpadding > 0 && ypadding > 0
                obj.padding = 'EWNS';
                obj.nPad = xpadding;
            elseif xpadding > 0
                obj.padding = 'NS';
                obj.nPad = xpadding;
            elseif ypadding > 0
                obj.padding = 'EW';
                obj.nPad = ypadding;
            end
            
            if zpadding > 0
                obj.padding = [obj.padding 'D'];
            end
            
        end
        
      
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods(Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [lat,lon,depth] = mesh(lims,top,bottom,nlat,nlon,nz,logdz)
            
            % [lat,lon,depth] = mesh(lims,top,bottom,nlat,nlon,nz,logdz)
            % 
            % uses the grid limits to make a regular lat/lon/depth mesh
            % (note - the values correspond to corners and need to be
            % converted to cell centers and distances for further use)
            % longitudes & latitudes are in degrees, vertical spacing in km.
            % padding can be added to any given grid using pad(); not done here.
            % optional logdz, if specified, makes depth distribution logarithmic;
            % (top - bottom) = sum_{k=0}^nz (startdz * logdz^k)
            %                = startdz * (logdz^(nz+1) - 1)/(logdz - 1)
            % therefore startdz = (top - bottom)*(logdz - 1)/(logdz^(nz+1) - 1)
            % and n = 0:nz-1; z = cumsum(startdz*(logdz.^n)).
            
            % [elev lon lat]=m_etopo2([-130 -118 40 50]);
            % minlon = -126;
            % maxlon = -122;
            % minlat = 42;
            % maxlat = 48;
            
            minlon = llgrid.lon180(lims.lonmin);
            maxlon = llgrid.lon180(lims.lonmax);
            minlat = lims.latmin;
            maxlat = lims.latmax;
            %avglon = (lims.lonmin+lims.lonmax)/2;
            avglat = (lims.latmin+lims.latmax)/2;
            
            EWfactor = (maxlon-minlon)/nlon;
            NSfactor = (maxlat-minlat)/nlat;
            UDfactor = (bottom-top)/nz;
                        
            dlat = (maxlat - minlat)/nlat; newlat = minlat:dlat:maxlat;
            dlon = (maxlon - minlon)/nlon; newlon = minlon:dlon:maxlon;
            
            dx = latlontools.kmPerDeg(avglat)*dlat;
            dy = latlontools.kmPerDeg(avglat)*cos(avglat*pi/180)*dlon;
            
            % quick check for consistency
            if min(dx) <= 0 || min(dy) <= 0
                error('Longitude and latitude coordinates are wrong on input to mesh().');
            end
            
            if (nargin <= 6 || logdz == 1) && (nz>0)
                % logdz = 1;
                dz(1:nz) = (bottom-top)/nz;
                disp(['Spacing: lon ' num2str(dy) ' km, lat ' num2str(dx) ...
                    ' km, vertical ' num2str(dz(1)) ' km']);
            elseif nz>0
                startdz = (bottom-top)*(logdz - 1)/(logdz^nz - 1);
                n = 0:nz-1; dz = startdz*(logdz.^n);
                disp(['Spacing: lon ' num2str(dy) ' km, lat ' num2str(dx) ...
                    ' km, vertical starts with ' num2str(dz(1)) ' km']);
            else
                dz = [];
                disp('Number of layers not a positive integer; no layers created');
            end
                              
            z = [top top+cumsum(dz)];
            
            lat = newlat';
            lon = newlon';
            depth = z';            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = merge(array)
            % [obj] = merge(obj,array)
            %   merges an array of llgrid objects into one big grid;
            %   very little error checking at present: only works for
            %   something reasonable
            %   given all the sorting etc, conversion to 0-360 longitude
            %   creates problems here if we cross zero longitude (results
            %   in a jump from 360 to zero that messes up the sorting);
            %   for now, leaving this out obj.lon = obj.lon360(obj.lon)
            %   and implicitly assuming that the arrays are all compatible.
            %   taking care at the grid boundaries to avoid duplicates.
            obj = array(1);
            eps = 1e-6;
            
            nodelat = obj.lat;
            nodelon = obj.lon;
            nodez = obj.depth;
            
            for i=2:length(array)
                newlat = array(i).lat;
                newlon = array(i).lon;
                newz = array(i).depth;
                if abs(newlat(1) - nodelat(end)) < eps
                    nodelat = [nodelat; newlat(2:end)];
                else
                    nodelat = [nodelat; newlat];
                end
                if abs(newlon(1) - nodelon(end)) < eps
                    nodelon = [nodelon; newlon(2:end)];
                else
                    nodelon = [nodelon; newlon];
                end
                if abs(newz(1) - nodez(end)) < eps
                    nodez = [nodez; newz(2:end)];
                else
                    nodez = [nodez; newz];
                end
            end
            
            nodelat = unique(sort(nodelat));
            nodelon = unique(sort(nodelon));
            nodez = unique(sort(nodez));
            
            obj.dlat = diff(nodelat);
            obj.dlon = diff(nodelon);
            obj.dz = diff(nodez);
            obj.lat = nodelat;
            obj.lon = nodelon;
            obj.depth = nodez;
            obj.limits.latmin = min(obj.lat);
            obj.limits.latmax = max(obj.lat);
            obj.limits.lonmin = min(obj.lon);
            obj.limits.lonmax = max(obj.lon);
            obj.limits.depthmin = min(obj.depth);
            obj.limits.depthmax = max(obj.depth);
            obj.nlat = length(obj.dlat);
            obj.nlon = length(obj.dlon);
            obj.nzEarth = length(obj.dz);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function lon = lon360(lon)
            %   given list of longitudes, convert them to [0,360)
            lon(lon<0) = lon(lon<0)+360;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function lon = lon180(lon)
            %   given list of longitudes, convert them to (-180,180]
            lon(lon>180) = lon(lon>180)-360;
        end
    end
end