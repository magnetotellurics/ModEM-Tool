classdef llmodel < modelplot
    %   addd depth and conductivity to the latgrid class, to store
    %   conductivity on a standard cartesian grid with lat/lon info, and
    %   tools for extracting slices in lat/lon/depth coordinates
    
%     properties
%         grid       % llgrid
%         modelType  = 'electrical conductivity'; % conductivity or resistivity
%         modelUnits = 'S/m'; % S/m or Ohm*m
%         paramType  % LOGE, LOG10 or LINEAR
%         v          % 3D array
%         AirCond = 1e-10;
%         SeaWaterCond = 4.8;
%     end
%     
    properties
        fileName
        fileHeader
        modelVariables
        geospatialInfo
    end
%
%     properties (SetAccess = protected)
%         % a non-regular mesh of X/Y values computed using a lat0/lon0
%         % origin in cartesian coords and used for interpolation to a
%         % regular x/y grid (class xygrid object)
%         %X
%         %Y
%         %Z
%         lat0
%         lon0
%         limits
%     end

    
    methods
        function [obj] = llmodel(varargin)
            %   class constructor   
            %   [obj] = llmodel(grid,value,paramType,modelType,modelUnits)
            %     OR
            %   [obj] = llmodel(grid,location,paramType,modelType,modelUnits)
            %         
            %   if no value, initialized with NaN's.
            %   example call to initialize without a value:
            %       obj = modelplot(grid,'CELL');
            %   grid always defined at cell centers

            obj = obj@modelplot(varargin{:});
            
            if nargin == 0
                obj.grid = llgrid;
                return
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = uiplot(obj,padding)
            % If conductivity or resistivity, convert to log10 for plotting
            
            Nx = obj.grid.nlat;
            Ny = obj.grid.nlon;
            Nz = obj.grid.nzEarth;
            obj.v = reshape(obj.v,[Nx Ny Nz]);
            if nargin > 1
                options.padding = padding;
            else
                options.padding = 0;
            end
            if strfind(obj.modelType,'conductivity')
                options.cblabel = 'log_{10} \sigma';
                options.clims = [-4 0];
            elseif strfind(obj.modelType,'resistivity')
                options.cblabel = 'log_{10} \rho';
                options.clims = [0 4];
            end
            % just in case, convert obj to log10
            pltobj = obj.log10;
            options.slice = 'Z';
            options.Np = 1;
            options.iXlim(1) = 1;
            options.iXlim(2) = Nx+1;
            options.iYlim(1) = 1;
            options.iYlim(2) = Ny+1;
            options.iZlim(1) = 1;
            options.iZlim(2) = Nz+1;
            options.latlon = 1; 
            CondPlotSet(pltobj.v,pltobj.grid,options);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function xyobj = xymodel(obj,mstruct,newgrid,varargin)
            % Usage:
            %    xyobj = xymodel(obj,mstruct,newgrid,varargin)
            %
            % convert model parameter to regular lat/lon grid
            %   defaults: largest set of lats and lons possible
            %
            % by default, will pad with two additional cells, then trim
            % these cells upon conversion to X/Y grid - this usually avoids
            % numerical issues at the boundaries
            %
            % generalized based on Ph.D. work of Han Qi to use the Matlab
            % Plotting Toolbox, which uses the mstruct structure to set up
            % the coordinate projection. This will not work without the
            % toolbox, but some basic methods will still be available.
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
            
            if nargin == 1
                if isempty(obj.lat0) || isempty(obj.lon0)
                    error('Please initialize the lat/lon origin first');
                else
                    lat0 = obj.lat0;
                    lon0 = obj.lon0;
                end
            end
            
            % origin in mstruct overrides origin in object
            if nargin > 1
                if isfield(mstruct,'origin')
                    lat0 = mstruct.origin(1);
                    lon0 = mstruct.origin(2);
                end
            end

            % make a default output grid if it is not provided
            if nargin > 2 && isobject(newgrid)
                cartgrid = newgrid;
            else
                cartgrid = xygrid(obj.grid,mstruct);
            end
            inUnits = cartgrid.units;
            
            % by default, set the method to nearest neighbour interpolation
            method = 'nearest';

            if nargin > 3
                n = length(varargin);
                if mod(n,2)
                    error('Optional arguments must occur in pairs')
                end
                for k = 1:2:n
                    option = lower(varargin{k});
                    switch option
                        case 'method'
                            method = varargin{k+1};                                
                        otherwise
                            error('Optional argument not defined')
                    end
                end
            end
            
            % initialize xyobj
            xyobj = xymodel; 
            xyobj.grid = cartgrid;
            xyobj.coords = 'cartesian';
            
            
            if strcmp(method,'nointerp')
                % no interpolation of the value; grid reshaped above
                % works if the grid sizes are compatible, no other checks
                if (cartgrid.nx == obj.grid.nlat) && (cartgrid.ny == obj.grid.nlon) && (cartgrid.nz == obj.grid.nz); then
                    [xyobj] = xymodel(xyobj.grid,obj.v,obj.paramType,obj.modelType,obj.modelUnits);
                    xyobj = setOrigin(xyobj,lat0,lon0);
                else                    
                    error('Unable to convert llmodel to xymodel without interpolation: grid sizes do not match');
                end                           
                                             
            else
                
                % convert all distances to meters, including vertical
                cartgrid = km2m(cartgrid);
                inGrid = km2m(obj.grid);
                
                % ALWAYS convert to linear for interpolation ... otherwise,
                % huge errors at discontinuities such as coast
                paramType = obj.paramType;
                if ~strcmp(paramType,'LINEAR')
                    obj = linear(obj);
                end
                
                % make a regular output mesh at cell centers
                %inGrid = inGrid.pad('NSEW',2);
                if strcmp(obj.location,'CELL')
                    % compute center locations from grid lon/lat/depth
                    ctrlon = inGrid.lon(1:end-1) + inGrid.dlon/2;
                    ctrlat = inGrid.lat(1:end-1) + inGrid.dlat/2;
                    ctrdepth = inGrid.depth(1:end-1) + inGrid.dz/2;
                else
                    % the values are given in the exact locations lon/lat/depth
                    ctrlon = inGrid.lon;
                    ctrlat = inGrid.lat;
                    ctrdepth = inGrid.depth;
                end
                nlon = cartgrid.ny;
                nlat = cartgrid.nx;
                ndepth = cartgrid.nzEarth;
                nTot = nlon*nlat*ndepth;
                [Y,X] = meshgrid(cartgrid.yctr,cartgrid.xctr);
                Z = cartgrid.zctr; %Z = Z(cartgrid.nzAir:end);

                % use map projection to make an irregular lat/lon grid
                if strcmp(mstruct.mapprojection,'latlon2xy')
                    xy = [reshape(X,[1,nlon*nlat]);reshape(Y,[1,nlon*nlat])];
                    ll = llgrid.xy2ll(xy,lat0,lon0);
                else
                    % use km
                    %Y = 1.0e-3*Y;
                    %X = 1.0e-3*X;
                    %Z = 1.0e-3*Z;
                    try
                        mstruct = defaultm(mstruct);
                    catch
                        error(['Unknown map projection: ',mstruct.mapprojection]);
                    end
                    switch strtrim(mstruct.mapprojection)
                        case 'utm'
                            %mstruct.maplatlimit = [min(ctrlat) max(ctrlat)];
                            %mstruct.maplonlimit = [min(ctrlon) max(ctrlon)];
                            [LATs,LONs] = minvtran(mstruct,Y,X);
                            ll(1,:) = reshape(LATs,1,nlon*nlat);
                            ll(2,:) = reshape(LONs,1,nlon*nlat);
                        otherwise  % e.g., 'eqdcylin','eqacylin','lambertstd','eqaazim'
                            [LATs,LONs] = projinv(mstruct,Y,X);
                            ll(1,:) = reshape(LATs,1,nlon*nlat);
                            ll(2,:) = reshape(LONs,1,nlon*nlat);
                    end
                end
                LATI = reshape(ll(1,:)'*ones(1,ndepth),nTot,1);
                LONI = reshape(ll(2,:)'*ones(1,ndepth),nTot,1);
                
                % set up a regular lat/lon mesh
                temp = reshape(Z,1,ndepth);
                ZI = reshape(ones(nlat*nlon,1)*temp,nTot,1);
                [LON,LAT,Z] = meshgrid(ctrlon,ctrlat,ctrdepth);
                
                % interpolate from regular to irregular lat/lon mesh
                xyobj.v = interp3(LON,LAT,Z,obj.v,LONI,LATI,ZI,method);
                xyobj.v = reshape(xyobj.v,[nlat,nlon,ndepth]);
                
                % put together the complete xyobj
                xyobj.modelType  = obj.modelType;
                xyobj.modelUnits = obj.modelUnits;
                xyobj.paramType  = obj.paramType;
                xyobj.AirCond  = obj.AirCond;
                xyobj.SeaWaterCond  = obj.SeaWaterCond;
                xyobj = setOrigin(xyobj,lat0,lon0);
                
                % NOW convert xyobj back to log if needed
                if strcmp(paramType,'LOG10')
                    xyobj = log10(xyobj);
                elseif strcmp(paramType,'LOGE')
                    xyobj = loge(xyobj);
                end
                
                % convert all distances back to km (projection in radians)
                % I don't understand this but that seems to only convert to
                % km, not to meters - there's a bug somewhere but leaving
                % it like this for now because I cannot find it
                xyobj.grid.dx = 1000*rad2km(xyobj.grid.dx);
                xyobj.grid.dy = 1000*rad2km(xyobj.grid.dy);
                xyobj.grid.origin = 1000*rad2km(xyobj.grid.origin);
                if contains(inUnits,'m')
                    xyobj.grid = km2m(xyobj.grid);
                end
                
                % trim the 2 padding cells back
                %xyobj = xyobj.trim('NSEW',2);
            end

        end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = setOrigin(obj,lat,lon)
            
            obj.lat0 = lat;
            obj.lon0 = lon;
            
        end
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = setLimits(obj)
            
            % obj = obj.setLimits;
            %
            % used to update or fix the limits based on obj.grid
            
            obj.limits.lonmin = min(obj.grid.lon);
            obj.limits.lonmax = max(obj.grid.lon);
            obj.limits.latmin = min(obj.grid.lat);
            obj.limits.latmax = max(obj.grid.lat);
            obj.limits.depthmin = min(obj.grid.depth);
            obj.limits.depthmax = max(obj.grid.depth);
            
        end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = setModelVariables(obj,modelvar)
            
            % obj = setModelVariables(obj,modelvar)
            %
            % sets the model variables needed for ModEM to function
            
            if nargin<2
                rotation = 0.0;
            elseif isstruct(modelvar)
                obj.modelVariables = modelvar;
                return
            else
                rotation = modelvar;
            end
            
            % add global variable primary_coords [latlon|xy]
            obj.modelVariables.primary_coords = 'latlon';
            
            % use corner coordinates to exactly compute grid cells from grid centers
            obj.modelVariables.corner_description = 'used to compute exact grid geometry from cell centers';
            obj.modelVariables.corner_location = 'upper_southwest';
            obj.modelVariables.corner_latitude = obj.grid.lat(1);
            obj.modelVariables.corner_longitude = obj.grid.lon(1);
            obj.modelVariables.corner_depth = obj.grid.depth(1);
            
            % model rotation complements data rotation information from data file
            obj.modelVariables.rotation_units = 'degrees';
            obj.modelVariables.rotation_angle = rotation;
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = trim(obj,nTrim,direction)
            
            % obj = trim(obj,nTrim,direction)
            %
            % use to trim the grid by nTrim cells in a particular direction 
            % direction is a string that contains letters W,E,N,S,U,D
            %
            % intended use is to trim artifacts at the grid boundaries
            % caused by conversions between spherical and cartesian coords,
            % or to remove padding. By default, removes padding.
            
            if nargin < 2
                nTrim = obj.grid.nPad;
            end
            
            if nargin < 3
                direction = obj.grid.padding;
            end
                        
            if size(direction,1)>1
                direction = direction';
            end
            
            [newgrid,ii,jj,kk] = trim(obj.grid,nTrim,direction);
            
            obj.v = obj.v(ii,jj,kk);
            obj.grid = newgrid;
            
            % update the limits
            obj.limits = obj.getLimits;
            
            % update model variables
            obj = obj.setModelVariables;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = regrid(obj,newgrid,bg)
            %   map the model to a new grid; based on interpCond_3D
            %
            % Argument bg gives the background absolute conductivity values to
            % be used to replace NaN's (if any) - could be a scalar or a vector
            % changing with depth. Use bg=0 if interpolating conductivity variations.
            % By default, using nearest neighbor average to fill in NaN's.
            %
            % If you want the NaN's, e.g. for mapping onto a bigger grid,
            % hack this by calling newmodel = mymodel.regrid(newgrid,NaN)
            % as is done for the llmodel.insert functionality.
            %
            % Often there is need to first convert to [-180 180) or [0 360] 
            
            oldgrid = obj.grid;
            oldv = obj.v;
            
            if max(newgrid.depth) > max(oldgrid.depth)
                warning('The vertical extent of your new grid greater than that of the original. Think again!');
            elseif min(newgrid.depth) < min(oldgrid.depth)
                warning('The old grid cannot be interpolated to the new near-surface layers. Fill them in meaningfully.');
            end
            
            if strcmp(obj.location,'CELL')
                [oldY,oldX,oldZ] = meshgrid(oldgrid.lon(1:end-1)+oldgrid.dlon/2,...
                    oldgrid.lat(1:end-1)+oldgrid.dlat/2,...
                    oldgrid.depth(1:end-1)+oldgrid.dz/2);
                
                [newY,newX,newZ] = meshgrid(newgrid.lon(1:end-1)+newgrid.dlon/2,...
                    newgrid.lat(1:end-1)+newgrid.dlat/2,...
                    newgrid.depth(1:end-1)+newgrid.dz/2);

            elseif strcmp(obj.location,'NODE')
                [oldY,oldX,oldZ] = meshgrid(oldgrid.lon,...
                    oldgrid.lat,...
                    oldgrid.depth);
                
                [newY,newX,newZ] = meshgrid(newgrid.lon,...
                    newgrid.lat,...
                    newgrid.depth);
                
            end
            
            newv = interp3(oldY,oldX,oldZ,oldv,newY,newX,newZ,'nearest');
           
            obj.grid = newgrid;
            obj.v = newv;
            obj.limits = newgrid.limits;
            % used to set the origin to the lower left corner probably for
            % compatibility with ModEMM - do not do this, bad idea - will
            % deal with it when we convert to ModEMM model format
            %obj.lat0 = obj.grid.lat(1);
            %obj.lon0 = obj.grid.lon(1);
            
            % replace NaN's, if any, with these values... use these for conductivity;
            % use zero for conductivity perturbations
            if nargin < 3
                obj = obj.fillNaNs;
            else
                obj = obj.fillNaNs('background',bg);
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = regrid2depth(obj,depth,bg)
            %   map the model to a new set of depths using 1D interp
            % (note that the depths define new layer boundaries even if the
            % model is defined in the cells - that is taken care of)
            %
            % Argument bg gives the background absolute conductivity values to
            % be used to replace NaN's (if any) - could be a scalar or a vector
            % changing with depth. Use bg=0 if interpolating conductivity variations.
            % By default, using nearest neighbor average to fill in NaN's.
            %
            % If you want the NaN's, e.g. for mapping onto a bigger grid,
            % hack this by calling newmodel = mymodel.regrid2depth(newdepth,NaN)
            % as is done for the llmodel.insert functionality.
            
            oldgrid = obj.grid;
            oldv = obj.v;
            
            if max(depth) > max(oldgrid.depth)
                warning('The vertical extent of your new grid greater than that of the original. Think again!');
            elseif min(depth) < min(oldgrid.depth)
                warning('The old grid cannot be interpolated to the new near-surface layers. Fill them in meaningfully.');
            end
            
            if strcmp(obj.location,'CELL')
                oldY = llgrid.lon360(oldgrid.lon(1:end-1)+oldgrid.dlon/2);
                oldX = oldgrid.lat(1:end-1)+oldgrid.dlat/2;
                oldZ = oldgrid.depth(1:end-1)+oldgrid.dz/2;
                newZ = depth(1:end-1)+diff(depth)/2;

            elseif strcmp(obj.location,'NODE')
                oldY = llgrid.lon360(oldgrid.lon);
                oldX = oldgrid.lat;
                oldZ = oldgrid.depth;
                newZ = depth;
                                
            end
            
            newv = zeros(length(oldX),length(oldY),length(newZ));
            for j = 1:length(oldY)
                for i = 1:length(oldX)
                    tmp = squeeze(oldv(i,j,:));
                    newv(i,j,:) = interp1(oldZ,tmp,newZ);
                end
            end
                                
            obj.grid.depth = depth;
            obj.grid.dz = diff(depth);
            obj.grid.nzEarth = length(obj.grid.dz);
            obj.grid.limits.depthmin = min(depth);
            obj.grid.limits.depthmax = max(depth);
            
            obj.v = newv;
            % used to set the origin to the lower left corner probably for
            % compatibility with ModEMM - do not do this, bad idea - will
            % deal with it when we convert to ModEMM model format
            %obj.lat0 = obj.grid.lat(1);
            %obj.lon0 = obj.grid.lon(1);
            
            % replace NaN's, if any, with these values... use these for conductivity;
            % use zero for conductivity perturbations
            if nargin < 3
                obj = obj.fillNaNs;
            else
                obj = obj.fillNaNs('background',bg);
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = trim2grid(obj)
            %   sets obj limits to grid limits
            
            obj.limits = obj.grid.limits;

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = regrid2latlon(obj,lat,lon,bg)
            %   map the model to a new set of lats and lons using 2D interp
            % (note that the lat & lon define new cell boundaries even if the
            % model is defined in the cells - that is taken care of)
            %
            % Argument bg gives the background absolute conductivity values to
            % be used to replace NaN's (if any) - could be a scalar or a vector
            % changing with depth. Use bg=0 if interpolating conductivity variations.
            % By default, using nearest neighbor average to fill in NaN's.
            %
            % If you want the NaN's, e.g. for mapping onto a bigger grid,
            % hack this by calling newmodel = mymodel.regrid(newgrid,NaN)
            % as is done for the llmodel.insert functionality.
            %
            % Often there is need to first convert to [-180 180) or [0 360] 
            
            oldgrid = obj.grid;
            oldv = obj.v;
            
            if max(lat) > max(oldgrid.lat) || min(lat) < min(oldgrid.lat)
                warning('The latitudes of your new grid extend outside of the model grid. Think again!');
            elseif max(lon) > max(oldgrid.lon) || min(lon) < min(oldgrid.lon)
                warning('The longitudes of your new grid extend outside of the model grid. Fill them in meaningfully.');
            end
            
            if strcmp(obj.location,'CELL')
                [oldY,oldX] = meshgrid(oldgrid.lon(1:end-1)+oldgrid.dlon/2,...
                    oldgrid.lat(1:end-1)+oldgrid.dlat/2);
                
                [newY,newX] = meshgrid(lon(1:end-1)+diff(lon)/2,...
                    lat(1:end-1)+diff(lat)/2);
                
                oldZ = oldgrid.depth(1:end-1)+oldgrid.dz/2;

            elseif strcmp(obj.location,'NODE')
                [oldY,oldX] = meshgrid(oldgrid.lon,...
                    oldgrid.lat);
                
                [newY,newX] = meshgrid(lon,...
                    lat);
                
                oldZ = oldgrid.depth;
                
            end
            
            newv = zeros(size(newY,1),size(newY,2),length(oldZ));
            for k = 1:length(oldZ)
                tmp = squeeze(oldv(:,:,k));
                newv(:,:,k) = interp2(oldY,oldX,tmp,newY,newX,'nearest');
            end
           
            obj.grid.lat = lat;
            obj.grid.dlat = diff(lat);
            obj.grid.nlat = length(obj.grid.dlat);
            obj.grid.limits.latmin = min(lat);
            obj.grid.limits.latmax = max(lat);

            obj.grid.lon = lon;
            obj.grid.dlon = diff(lon);
            obj.grid.nlon = length(obj.grid.dlon);
            obj.grid.limits.lonmin = min(lon);
            obj.grid.limits.lonmax = max(lon);

            obj.v = newv;
            obj.limits = obj.grid.limits;
            if min(lat) == -90 && max(lat) == 90 && min(lon) <= 0 && max(lon) >= 360
                obj.isglobal = 1;
            else
                obj.isglobal = 0;
            end
            
            % used to set the origin to the lower left corner probably for
            % compatibility with ModEMM - do not do this, bad idea - will
            % deal with it when we convert to ModEMM model format
            %obj.lat0 = obj.grid.lat(1);
            %obj.lon0 = obj.grid.lon(1);
            
            % replace NaN's, if any, with these values... use these for conductivity;
            % use zero for conductivity perturbations
            if nargin < 4
                obj = obj.fillNaNs;
            else
                obj = obj.fillNaNs('background',bg);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = node2cell(obj)
            %   for model values defined on NODES of the grid, 
            %   interpolate to CELL CENTERS; grid or origin do not change.
            %   to interpolate to a different grid (keeping the CELL or
            %   NODE convention) use obj.regrid(newgrid,bg)
            
            if ~strcmp(obj.location,'NODE')
                error('The object is not defined on grid nodes, can''t interpolate to cells.');
            end
            
            oldv = obj.v;
            mygrid = obj.grid;
            
            [oldY,oldX,oldZ] = meshgrid(llgrid.lon360(mygrid.lon),...
                mygrid.lat,...
                mygrid.depth);
                        
            [newY,newX,newZ] = meshgrid(llgrid.lon360(mygrid.lon(1:end-1)+mygrid.dlon/2),...
                mygrid.lat(1:end-1)+mygrid.dlat/2,...
                mygrid.depth(1:end-1)+mygrid.dz/2);
                        
            newv = interp3(oldY,oldX,oldZ,oldv,newY,newX,newZ);
                       
            obj.v = newv;
            obj.location = 'CELL';
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cell2node(obj,bg)
            %   for model values defined on CELL CENTERS of the grid, 
            %   interpolate to NODES; grid or origin do not change.
            %   to set the boundaries, use bg 
            %   (defaults to depth averaged value).
            %   can also set bg = NaN.
            %   to interpolate to a different grid (keeping the CELL or
            %   NODE convention) use obj.regrid(newgrid,bg)
            
            if nargin < 2
                bg = squeeze(mean(mean(obj.v)));
            end 
            
            if ~strcmp(obj.location,'CELL')
                error('The object is not defined on grid cells, can''t interpolate to nodes.');
            end
            
            oldv = obj.v;
            mygrid = obj.grid;
            
            [oldY,oldX,oldZ] = meshgrid(llgrid.lon360(mygrid.lon(1:end-1)+mygrid.dlon/2),...
                mygrid.lat(1:end-1)+mygrid.dlat/2,...
                mygrid.depth(1:end-1)+mygrid.dz/2);
                        
            [newY,newX,newZ] = meshgrid(llgrid.lon360(mygrid.lon),...
                mygrid.lat,...
                mygrid.depth);
                        
            newv = interp3(oldY,oldX,oldZ,oldv,newY,newX,newZ);
 
            % replace NaN's, if any, with these values... use these for conductivity;
            % use zero for conductivity perturbations
            nz = length(newgrid.dz);
            if isscalar(bg)
                bg(1:nz) = bg;
            end
            if strcmp(obj.paramType,'LOGE') && (min(bg) > 0)
                bg = log(bg);
            end
            if strcmp(obj.paramType,'LOG10') && (min(bg) > 0)
                bg = log10(bg);
            end
            
            for k = 1:nz
                temp = squeeze(newv(:,:,k));
                temp(isnan(temp)) = bg(k);
                newv(:,:,k) = temp;
            end
            
            obj.v = newv;
            obj.location = 'NODE';
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,dtotal] = croptodata(obj,datalon,datalat,km)
            %
            %   [obj,dtotal] = croptodata(obj,datalon,datalat,distance)
            %
            % crop the llmodel to data locations by replacing all cells
            % greater than a distance (in km) from a data point by NaNs
            % optionally output the number of sites that are within the km
            % distance from each grid cell
            
            % first convert the distance from km to degrees
            %kmPerDeg = latlontools.kmPerDeg(mean(datalat));
            %delta = km/kmPerDeg;
            
            if strcmp(obj.location,'CELL')
                vlon = obj.grid.lon(1:end-1) + obj.grid.dlon/2;
                vlat = obj.grid.lat(1:end-1) + obj.grid.dlat/2;
            elseif strcmp(obj.location,'NODE')
                vlon = obj.grid.lon;
                vlat = obj.grid.lat;
            end
            [vlon,vlat] = meshgrid(vlon,vlat);
            vlon = reshape(vlon,size(obj.v,2)*size(obj.v,1),1);
            vlat = reshape(vlat,size(obj.v,2)*size(obj.v,1),1);

            % make matrix of distances using haversine formula for
            % great-circle distance delta(i,j) = 2*R_earth*arcsin(sqrt(d(i,j))
            R_earth=6372.795;
            h=(sin(km/(2*R_earth)))^2;
            d = zeros(length(datalon),length(vlon));
            for i=1:length(datalon)
                for j=1:length(vlon)
                    phi1=datalon(i)*pi/180; theta1=datalat(i)*pi/180;
                    phi2=vlon(j)*pi/180; theta2=vlat(j)*pi/180;
                    dphi=phi1-phi2; dtheta=theta1-theta2;
                    d(i,j) = (sin(dtheta/2))^2+cos(theta1)*cos(theta2)*(sin(dphi/2))^2;
                end
            end
            
            % set values for which there are no distances within km to NaNs
            ii = d<h;
            dlogical = zeros(length(datalon),length(vlon));
            dlogical(ii) = 1;
            dtotal = sum(dlogical);
            dtotal = reshape(dtotal,size(obj.v,2),size(obj.v,1));
            %figure; imagesc(dtotal); colorbar
            jj = dtotal==0;
            for k = 1:size(obj.v,3)
                value = obj.v(:,:,k);
                value(jj) = NaN;
                obj.v(:,:,k) = value;
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = insert(obj,smallobj)
            %   insert a smaller llmodel into a bigger llmodel:
            %
            % concise and inefficient but OH MY GOD HOW POWERFUL!
            %
            
            if ~strcmp(obj.location,smallobj.location)
                if strcmp(obj.location,'CELL')
                    myobj = smallobj.node2cell;
                elseif strcmp(obj.location,'NODE')
                    myobj = smallobj.cell2node(NaN);
                end
            else
                myobj = smallobj;
            end
            
            temp = myobj.regrid(obj.grid,NaN);
            
            validind = find(~isnan(temp.v));
            N = length(validind);
            if N > 0
                disp(['Inserted ' num2str(N) ' valid ' obj.location ' values from the small model into large']);
            else
                disp('Model insert unsuccessful: no values were changed in the big model. Check for zero longitude problems.');
            end
            obj.v(validind) = temp.v(validind);
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function Data = contour(obj,name,fname)
            %
            % Data = contour(obj,name,fname)
            % creates a 2D contour to save as a shapefile and for plotting.
            % Uses obj.v to find all non-NaN values in layer 1 and contour 
            % around them. Can be easily extended to plot contours at all
            % depths or at a selected depth.
            % If fname is supplied, saves to a shapefile.
            % 
            % Equivalent to croptodata2shp.m script and (roughly) to
            % Q = bwtraceboundary( bwmorph( isnan( YourMatrix ), 'thicken', 1) );
            % in the Image Processing Toolbox.
            
            if nargin<2
                name = 'Crop to Data Contour';
            end
            
            if nargin==3
                data2shp = 1;
            else
                data2shp = 0;
            end
            
            % use layer 1 and the middles of the cells
            v = squeeze(obj.v(:,:,1)); % i = lat; j = lon
            midlon = obj.grid.lon(1:end-1)+obj.grid.dlon/2;
            midlat = obj.grid.lat(1:end-1)+obj.grid.dlat/2;
            
            % find all non-NaN nodes and use the boundary function
            ii=0;
            lon = zeros(obj.grid.nlat*obj.grid.nlon,1);
            lat = zeros(obj.grid.nlat*obj.grid.nlon,1);
            for i=1:obj.grid.nlat
                for j=1:obj.grid.nlon
                    if ~isnan(v(i,j))
                        ii=ii+1;
                        lon(ii) = midlon(j);
                        lat(ii) = midlat(i);
                    end
                end
            end
            lon = lon(1:ii);
            lat = lat(1:ii);
            k = boundary(lon,lat,1);
            
            % save as a shapefile structure; use shapewrite(Data, fname)
            Data.Geometry = 'Polygon' ;
            Data.BoundingBox(:,1) = [obj.grid.lat(1) obj.grid.lat(end)];
            Data.BoundingBox(:,2) = [obj.grid.lon(1) obj.grid.lon(end)];
            Data.X = lon(k) ;
            Data.Y = lat(k) ;
            Data.Name = name ;
            
            % save to shapefile if desired
            if data2shp
                shapewrite(Data,fname);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h1 = xsectionPlot(obj,yline,xline,maxdepth,P,pts,marks)
            % h1 = xsectionPlot(obj,yline,xline,maxdepth,P,pts,marks)
            % 
            % Plot a model cross section at an angle
            % yline defines Y or longitude coordinates of line ends
            % xline defines X or latitude coordinates of line ends
            % dy defines new cell distances in longitude (default grid.dlon)
            % dx defines new cell distances in latitude (default grid.dlat)
            % pts is an optional array of distances from zero which should be marked
            % marks are one-character text strings
            %
            % adapted from my [h1,cb] = CondPlotSection(h0,h1,yline,xline,pts,marks)
            % which was initially used to plot a cross-section using
            % CondPlot interactive plotting. This requires the Image
            % Processing toolbox, so here I've adapted it to merely use the
            % pre-defined corners of a cross-section.
            % Example interactive usage:
            % P  = get(h0,'UserData');
            % ax = get(h0,'CurrentAxes');
            % [yline,xline] = getline(ax);
            % Then, if model is in X/Y, convert to lat/lon. Here, doing
            % everything in lat/lon and assuming an llmodel as input.
       
            ii = (obj.grid.lon >= yline(1)) & (obj.grid.lon <=yline(2));
            jj = (obj.grid.lat >= xline(1)) & (obj.grid.lat <=xline(2));
            nres = max(sum(ii),sum(jj));

            dy = (yline(2) - yline(1))/nres;
            newy = yline(1):dy:yline(2);
            % for a lonSlice
            if isempty(newy)
                newy(1:nres+1) = yline(1);
            end
            lon = yline;
            
            dx = (xline(2) - xline(1))/nres;
            newx = xline(1):dx:xline(2);
            % for a latSlice
            if isempty(newx)
                newx(1:nres+1) = xline(1);
            end
            lat = xline;
            
            if nargin < 4
                maxdepth = obj.grid.depth(end-1);
            end
            
            if nargin < 5
                P.Clims = [-4 0];
                P.cmap = 'jet';
                P.cblabel = 'log_{10} \sigma';
                P.PlotColorBar = 1;
                P.flipud = 0;
                P.Label = 0;
            end
            
            oldy = obj.grid.lon(1:end-1) + obj.grid.dlon/2;
            oldx = obj.grid.lat(1:end-1) + obj.grid.dlat/2;
            oldz = obj.grid.depth(1:end-1) + obj.grid.dz/2;
            oldv = obj.v;
            
            %[oldY,oldX,oldZ] = meshgrid(oldy,oldx,oldz);
            [oldY,oldX] = meshgrid(oldy,oldx);
            
            ii = obj.grid.depth <=maxdepth;
            newz = obj.grid.depth(ii) + obj.grid.dz(ii)/2;
            
            newv(length(newx),length(newz)) = 0;
            for k = 1:length(newz)
                for i = 1:length(newx)
                    temp = interp2(oldY,oldX,squeeze(oldv(:,:,k)),newy(i),newx(i));
                    newv(i,k) = temp;
                end
            end
            %for i = 1:length(newx)
            %    temp = interp3(oldY,oldX,oldZ,oldv,newy(i),newx(i),newz);
            %    newv(i,:) = squeeze(temp)';
            %end            
            [Y,X] = meshgrid(newz,sqrt((newx-xline(1)).^2+(newy-yline(1)).^2));
            
            h1 = figure('Position',[200,400,1400,600],...
                    'Units','normalized','Color','w');
            pcolor(X,Y,newv);caxis(P.Clims);colormap(P.cmap);
            if P.PlotColorBar
                cb=colorbar('FontSize',14);
                ylabel(cb,P.cblabel,'FontWeight','demi','FontSize',14);
            end
            if P.flipud; colormap(flipud(P.cmap)); end
            axis ij; shading flat
            if nargin > 5
                hold on;
                text(pts,-10*ones(size(pts)),marks,'FontWeight','demi','FontSize',14);
            end
            hold off
            set(gca,'FontWeight','demi','FontSize',16)
            if P.Label
                xlabel(['Distance along the profile (' P.units ')'])
                ylabel(['Depth (' P.units ')'])
            end
            
            if nargin < 6
                pt1 = ['( ' num2str(lat(1)) ', ' num2str(lon(1)) ' )'];
                text(-0.05,1.05,pt1,'FontWeight','demi','FontSize',10,...
                    'Color','b','units','normalized');
                pt2 = ['( ' num2str(lat(2)) ', ' num2str(lon(2)) ' )'];
                text(0.95,1.05,pt2,'FontWeight','demi','FontSize',10,...
                    'Color','b','units','normalized');
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [h,X] = pcolor(obj,type,value,lims,P)
            % Sample usage: h = obj.pcolor('depth',16)
            % Optional structure P contains plotting parameters.
            % By default does 1D interpolation which can be switched off
            % by setting P.nointerpbyvalue = 1.
            % Full example:
            % obj = llmodel.read('SEMUM_percent.nc','netcdf');
            % [lims,posn] = latlontools.getLimits('global');
            % obj.pcolor('depth',150,lims);
            % colormap(flipud(jet));
            % obj.plot_boundaries(gcf,'tectonic_boundaries_by_zone')
            % obj.plot_boundaries(gcf,'gshhg_global_country_borders.shp')
            % obj.plot_boundaries(gcf,'gshhg_usa_states_pol.shp')
            % obj.plot_boundaries(gcf,'physio.shp')
            % set(gcf,'PaperPosition',posn);
            %
            % Note that at least the conductivity models are piecewise
            % constant in depth, and the grid specifies the boundaries...
            % while for lat/lon, the global conductivity models are
            % interpolated to a set of nodes in the process of GM to GG
            % conversion... need to be careful with this and should
            % sometime revisit.
            
            if strcmp(obj.location,'CELL')
                lon = obj.grid.lon(1:end-1) + obj.grid.dlon/2;
                lat = obj.grid.lat(1:end-1) + obj.grid.dlat/2;
                depth = obj.grid.depth(1:end-1) + obj.grid.dz/2;
            elseif strcmp(obj.location,'NODE')
                lon = obj.grid.lon;
                lat = obj.grid.lat;
                depth = obj.grid.depth;
            else
                error('Set model location to CELL or NODE before plotting.')
            end
                
            nlon = length(lon);
            nlat = length(lat);
            nzEarth = length(depth);
           
            if strcmp(type,'depth') 
                if issorted(depth)
                    k2 = find(depth>value,1);
                    k = max(k2-1,1);
                else
                    k2 = find(depth>value,1,'last');
                    k = min(k2+1,nzEarth);
                end
                alpha = (depth(k2)-value)/(depth(k2)-depth(k));
                X = obj.v(:,:,k);
                nointerp = 0;
                if nargin >= 5
                    if isfield(P,'nointerpbyvalue')
                        nointerp = P.nointerpbyvalue;
                    end
                end
                if ~nointerp
                    X2 = obj.v(:,:,k2);
                    X = alpha*X + (1-alpha)*X2;
                    mydepth = num2str(value);
                    disp(['Interpolated to depth ' num2str(value) ' km from depths '...
                        num2str(depth(k)) ' and ' num2str(depth(k2)) ' km']);
                else
                    mydepth = num2str(depth(k));
                    disp(['Plotting the value for depth ' num2str(depth(k)) ' km']);
                end
           elseif strcmp(type,'layer')
                k = value;
                X = obj.v(:,:,k);
                if isfield(obj,'depth_pr') % i.e., global model object
                    if length(obj.depth_pr) == size(obj.v,3)
                        disp(['Plotting the value for layer down to depth ' num2str(obj.depth_pr(k)) ' km']);
                    else
                        disp(['Plotting the value for layer index ' num2str(k) '; depth unknown']);
                    end
                else
                    disp(['Plotting the model value for layer index ' num2str(k)]);
                end
           elseif strcmp(type,'value')
                X = value;
                if size(X,1) == 1 && size(X,2) == 1
                    X = value * ones(nlat,nlon);
                end
                if size(X,1) ~= length(lat) || size(X,2) ~= length(lon)
                    error('Format of the value to plot should be X(latitude,longitude)')
                end
            else
                error('Not implemented yet')
            end
 
            if size(X,1) ~= length(lat) || size(X,2) ~= length(lon)
                error('Check that the size of the model (CELL or NODE) is compatible with the grid.')
            end
            
            % deal with Inf and -Inf which sometimes appear during
            % interpolation... not sure why that happens
            ii = find(X < -1e9); X(ii) = NaN;
            ii = find(X > 1e9); X(ii) = NaN;

            proj = 'miller';
            delta = 1e-6;
            clims = [min(min(X))-delta max(max(X))+delta];
            cmap = 'jet';
            nointerp = 0;
            showmesh = 0;
            regional = 1;
            clong = 180;
            titleOff = 0;
            labelsOff = 0;
            latlonOff = 0;
            if nargin >= 5
                if isfield(P,'proj')
                    proj = P.proj;
                end
                if isfield(P,'Clims')
                    clims = P.Clims;
                elseif isfield(P,'clims')
                    clims = P.clims;
                end
                if isfield(P,'cmap')
                    cmap = P.cmap;
                end
                if isfield(P,'nointerp')
                    nointerp = P.nointerp;
                end
                if isfield(P,'showmesh')
                    showmesh = P.showmesh;
                end
                if isfield(P,'clong')
                    clong = P.clong;
                end
                if isfield(P,'titleOff')
                    titleOff = P.titleOff;
                end
                if isfield(P,'labelsOff')
                    labelsOff = P.labelsOff;
                end
                if isfield(P,'latlonOff')
                    latlonOff = P.latlonOff;
                end
            end
            
            if nargin < 4
                lonmin = obj.limits.lonmin;
                lonmax = obj.limits.lonmax;
                latmin = obj.limits.latmin;
                latmax = obj.limits.latmax;
            elseif isstruct(lims)
                lonmin = lims.lonmin;
                lonmax = lims.lonmax;
                latmin = lims.latmin;
                latmax = lims.latmax;
            else
                lonmin = lims(1);
                lonmax = lims(2);
                latmin = lims(3);
                latmax = lims(4);
            end
            
            if strcmp(proj,'hammer')
                regional = 0;
                lonmin = lonmin - 180 + clong;
                lonmax = lonmax - 180 + clong;
            end
                        
            posn = [1,1,12,7];
            h = figure('Position',100*posn,...
                'PaperPosition',posn,...
                'PaperOrientation','Portrait',...
                'Color',[1 1 1]);
            % buggy in m_map - for now, don't plot GM continents
%             if strcmp(obj.coords,'geographic')
%                 m_coord('geographic');
%             else
%                 m_coord('IGRF2000-geomagnetic');
%             end
%             switch proj
%                 case 'miller'
%                     m_proj(proj,...
%                         'long',clong,...
%                         'lat', [latmin latmax]);
%                 case 'hammer'
%                     m_proj(proj,...
%                         'long',[lonmin lonmax],'clong',clong,...
%                         'lat', [latmin latmax]);
%                 otherwise
                    m_proj(proj,...
                        'long',[lonmin lonmax],...
                        'lat', [latmin latmax]);
%             end
            m_pcolor(llgrid.lon360(lon),lat,X); hold on;
            m_pcolor(llgrid.lon360(lon)-360,lat,X); 
            caxis(clims); colormap(cmap);
            if nointerp; shading flat; end
            if ~nointerp; shading interp; end
            if showmesh; shading faceted; end
            if ~labelsOff; colorbar('FontSize',28); end
            hold on;
            if strcmp(obj.coords,'geographic') && ~labelsOff
                cbndry  = [0.1 0.1 0.1];
                %m_gshhs_i('color',cbndry);              % Coastline... higher resolution
                %m_gshhs_i('speckle','color',cbndry);    % with speckle added
                m_coast('color',cbndry);              % Coastline...
                m_coast('speckle','color',cbndry);    % with speckle added
                %plotStates
            end
            if ~isempty(strfind(lower(proj),'miller')) ...
                    || ~isempty(strfind(lower(proj),'equidistant'))
                if ~labelsOff && ~latlonOff
                m_grid('box','fancy','tickdir','in',... %'box','on',...
                    'xtick',2,'linestyle','none',...
                    'XaxisLocation','top',...
                    'fontsize',28,'fontweight','demi');
                else
                m_grid('box','fancy','tickdir','in',... %'box','on',...
                    'xtick',0,'xticklabels',[],'ytick',0,'yticklabels',[],'linestyle','none',...
                    'XaxisLocation','top',...
                    'fontsize',28,'fontweight','demi');
                end
            elseif regional
                m_grid('tickdir','out','yaxislocation','left',...
                    'XaxisLocation','bottom','xlabeldir','end','ticklen',.02,...
                    'fontsize',28,'fontweight','demi');
            else
                m_grid('box','on','tickdir','in',...
                    'xtick',12,'ytick',[-45 -30 0 30 45],...
                    'XaxisLocation','middle','xticklabels','',...
                    'fontsize',28,'fontweight','demi');
            end
            %m_grid('box','fancy','tickdir','in','fontsize',14,'fontweight','demi');
            if strcmp(type,'depth')
                if ~titleOff
                    title([obj.modelName ': ' mydepth ' km'],'fontsize',38,'fontweight','demi');
                end
            end

        end

        function [h] = map(obj,lims,P)
            % Sample usage: h = obj.map(lims,P) or h = obj.map
            % Optional structure P contains plotting parameters.
            % Same as pcolor but used for plotting an empty map with data.
            % Overloads geoplotaddons.map
            %
            % Full example:
            % [lims,posn] = latlontools.getLimits('SouthAmerica');
            % obj = llmodel;
            % obj.map(lims);
            % obj.plot_boundaries(gcf,'tectonic_boundaries_by_zone')
            % obj.plot_boundaries(gcf,'gshhg_global_country_borders.shp')
            % obj.plot_boundaries(gcf,'gshhg_usa_states_pol.shp')
            % obj.plot_boundaries(gcf,'physio.shp')
            % set(gcf,'PaperPosition',posn);                  
                        
            if nargin == 1
                lims = obj.limits;
            end
            
            tmpobj = geoplotaddons;
            tmpobj.coords = obj.coords;
                
            if nargin >= 3
                h = tmpobj.map(lims,P);
            else
                h = tmpobj.map(lims);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function newobj = select(obj,lims)
            %
            % newobj = select(obj,lims)
            %
            % select a 3D area based on lims if present; otherwise, use the
            % limits from obj
            
            if nargin < 2
                lims = obj.limits;
            elseif ischar(lims)
                region = lims;
                lims = latlontools.getLimits(region);
                lims.depthmin = 0;
                lims.depthmax = sum(obj.grid.dz);
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
                    lims.depthmin = 0;
                    lims.depthmax = sum(obj.grid.dz);
                end
            end
            k1 = find(obj.grid.depth >= lims.depthmin, 1, 'first');
            k2 = find(obj.grid.depth <= lims.depthmax, 1, 'last');
            i1 = find(obj.grid.lat >= lims.latmin, 1, 'first');
            i2 = find(obj.grid.lat <= lims.latmax, 1, 'last');
            if issorted(llgrid.lon360(obj.grid.lon))
                j1 = find(llgrid.lon360(obj.grid.lon) >= llgrid.lon360(lims.lonmin), 1, 'first');
                j2 = find(llgrid.lon360(obj.grid.lon) <= llgrid.lon360(lims.lonmax), 1, 'last');
            else
                j1 = find(llgrid.lon180(obj.grid.lon) >= llgrid.lon180(lims.lonmin), 1, 'first');
                j2 = find(llgrid.lon180(obj.grid.lon) <= llgrid.lon180(lims.lonmax), 1, 'last');
            end
            
            value = obj.v(i1:i2-1,j1:j2-1,k1:k2-1);
            
            outgrid = obj.grid;
            outgrid.depth = obj.grid.depth(k1:k2);
            outgrid.lat = obj.grid.lat(i1:i2);
            outgrid.lon = obj.grid.lon(j1:j2);
            outgrid.dz = diff(outgrid.depth);
            outgrid.dlat = diff(outgrid.lat);
            outgrid.dlon = diff(outgrid.lon);
            outgrid.nzEarth = length(k1:k2-1);
            outgrid.nlat = length(i1:i2-1);
            outgrid.nlon = length(j1:j2-1);
            outgrid = outgrid.pad('NSEW',0); % very important! no padding info for selection
            outgrid.limits = lims;
            
            newobj = obj;
            newobj.grid = outgrid;
            newobj.v = value;
            newobj.limits = lims;
            
%             lati = intersect(find(lat >= latlim(1)),find(lat <= latlim(2)));
%             loni = intersect(find(lon >= lonlim(1)),find(lon <= lonlim(2)));
%             area = pvalue(lati,loni);
%             disp(['Average resistivity of ' headers(ip,:) ': ' num2str(1/(10^nanmean(nanmean(area)))) ' Ohm m']);

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function write(obj,fname,format,lims)
            %
            % Usage: write(obj,fname,format,lims)
            %
            % format options: Slices, NetCDF, ModEM (assumes WS lat/lon)
            % if lims are specified, saves a subset of model values
            % using the new limits as model bounds
            %
            % for Slices, both NODE and CELL formats are allowed...
            % other file formats, including NetCDF, currently assume CELLS.
            % This is restrictive, but on the other hand it simplifies
            % sharing. Might want to rethink this later.
            
            if nargin < 3
                format = 'ModEM';
            end
            if nargin < 4
                lims = obj.limits;
            end
            myobj = select(obj,lims);
            myvalue = myobj.v;
            mygrid = myobj.grid;
            switch lower(format)
                case 'slices'
                    if ~exist(fname,'dir'); mkdir(fname); end
                    for k = 1:mygrid.nzEarth
                        if strcmp(mygrid.location,'CELL')
                            str = sprintf('%07.2f',mygrid.depth(k)+mygrid.dz(k)/2);
                            nlat = mygrid.nlat;
                            nlon = mygrid.nlon;
                        elseif strcmp(mygrid.location,'NODE')
                            str = sprintf('%07.2f',mygrid.depth(k));
                            nlat = mygrid.nlat+1;
                            nlon = mygrid.nlon+1;
                        end
                        f = fopen([fname '/' fname '_' str 'km.slice'],'w');
                        fprintf(f,'# %d %d %s\n',nlat,nlon,mygrid.location);
                        if strcmp(mygrid.location,'CELL')
                            for i = 1:mygrid.nlat
                                for j = 1:mygrid.nlon
                                    fprintf(f,'%f\t%f\t%g\n',mygrid.lat(i)+mygrid.dlat(i)/2,mygrid.lon(j)+mygrid.dlon(j)/2,myvalue(i,j,k));
                                end
                            end
                        elseif strcmp(mygrid.location,'NODE')
                            for i = 1:mygrid.nlat+1
                                for j = 1:mygrid.nlon+1
                                    fprintf(f,'%f\t%f\t%g\n',mygrid.lat(i),mygrid.lon(j),myvalue(i,j,k));
                                end
                            end
                        end
                        fclose(f);
                    end
                case 'netcdf'
                    % save information for writing to file; note that we
                    % are keeping the original geospatial info limits even
                    % if we are only writing a subset area of the model
                    % this could be changed if needed
                    if isempty(obj.fileHeader)
                        error('Unable to write the model to NetCDF: please provide a file header. Start with netcdfiris.initHeader');
                    end
                    header = obj.fileHeader;
                    if isempty(obj.modelVariables)
                        obj = obj.setModelVariables; 
                        modelvar = obj.modelVariables; disp(modelvar);
                    else
                        modelvar = obj.modelVariables;
                    end
                    if isempty(obj.geospatialInfo)
                        geospatial = netcdfiris.initGeospatial(obj.limits,[median(obj.grid.dlon) median(obj.grid.dlat)]); disp(geospatial);
                    else
                        geospatial = obj.geospatialInfo;
                    end
                    if ~strfind(obj.location,'CELL')
                        error('NetCDF model format assumes models defined in CELLS. Run obj.node2cell or expand the format.');
                    else
                        latitude = mygrid.lat(1:end-1) + mygrid.dlat/2;
                        longitude = mygrid.lon(1:end-1) + mygrid.dlon/2;
                        depth = mygrid.depth(1:end-1) + mygrid.dz/2;
                    end
                    vars(1).short_name = obj.modelName;
                    vars(1).long_name = obj.modelType;
                    if strcmp(obj.paramType,'LOG10')
                        vars(1).display_name = ['log(10) ' obj.modelType ', in ' obj.modelUnits];
                    elseif strcmp(obj.paramType,'LOGE')
                        vars(1).display_name = ['ln ' obj.modelType ', in ' obj.modelUnits];
                    else
                        vars(1).display_name = [obj.modelType ', in ' obj.modelUnits];
                    end
                    vars(1).units = obj.modelUnits;
                    vars(1).value = permute(myvalue,[2 1 3]);
                    disp('Values permuted to conform to X(longitude,latitude,depth)');
                    % now, open NetCDF file and write - use a new name
                    ncid = netcdfiris.open([fname '.nc']);
                    netcdfiris.putHeader(ncid,header,geospatial,lims);
                    netcdfiris.putModelVariables(ncid,modelvar);
                    netcdfiris.putPoints(ncid,longitude,latitude,depth);
                    netcdfiris.putValue(ncid,vars);
                    netcdf.close(ncid);
                    disp(['Model object written to NetCDF file ' fname '.nc']);
                case 'modem'
                    % first convert from NODE to CELL if necessary
                    if strfind(obj.location,'NODE')
                        obj = obj.node2cell;
                    end
                    % now convert from LOG10 to LOGE at least until LOG10
                    % finally gets implemented in ModEM
                    if strfind(obj.paramType,'LOG10')
                        obj = obj.loge;
                    end                    
                    % writing to WS format modified for lat/lon spacing
                    nzAir = obj.grid.nzAir;
                    dx = obj.grid.dlat;
                    dy = obj.grid.dlon;
                    dz = obj.grid.dz;  
                    % here, we redefine the meaning of "origin" and write
                    % down the lower left corner instead - that's how we
                    % reconstruct the model when we read it in
                    origin = [obj.grid.lat(1) llgrid.lon180(obj.grid.lon(1)) 0];
                    rotation = 0;                    
                    if isfield(struct(obj.grid),'units')
                        if strcmp(obj.grid.units,'km')
                            % convert vertical grid spacing to meters!
                            dz = 1000*dz;
                        end
                    end                    
                    if strfind(obj.modelType,'conductivity')
                        switch obj.paramType
                            case {'LOGE','LOG10'}
                                rho = - obj.v;
                            otherwise
                                rho = 1./(obj.v);
                        end
                    else
                        rho = obj.v;
                    end
                    type = obj.paramType;
                    status = write_WS3d_model(fname,dx,dy,dz,rho,nzAir,type,origin,rotation);
                otherwise
                    error(['Can''t write to lat/lon model format ' format ': method unknown']);
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make a copy of a handle object, but store as a child object.
        function new = copy(this,new)
            
            if nargin == 1
                % Instantiate new object of the same class.
                new = feval(class(this));
            elseif isempty(strfind(class(new),'model'))
                warning('Trying to hard copy an llmodel object into an inacceptible class');
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
    
    end
    
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,vars] = read(cfile,format,isglobal)
            %
            % Usage:  [obj,vars] = read(cfile,format,isglobal)
            %
            % Reads in any llmodel object in formats Slices, NetCDF
                        
            if ~exist(cfile,'file')
                error(['Model file ' cfile ' not found']);
            end
            
            if nargin < 3
                isglobal = 0;
            end
            
            switch lower(format)
                case 'slices'
                    if ~exist(cfile,'dir')
                        error(['Directory ' cfile ' not found']); 
                    end
                    % get dimensions and lat/lon from first file
                    slices = findfiles('slice',cfile);
                    if isempty(slices)
                        slices = findfiles('cond',cfile);
                    end
                    ndepth = length(slices);
                    f = fopen(char(slices(1)),'r');
                    dims = textscan(f,'# %f %f %s',1);
                    nlat = cell2mat(dims(1));
                    nlon = cell2mat(dims(2));
                    if ~isempty(dims(3))
                        location = cell2char(dims(3));
                    else
                        location = 'CELL';
                    end
                    temp = textscan(f,'%f %f %f');
                    fclose(f);
                    latitude = reshape(temp{1},[nlon nlat]).';
                    longitude = reshape(temp{2},[nlon nlat]).';
                    depth = zeros(1,ndepth);
                    value = zeros(nlat,nlon,ndepth);
                    for k = 1:length(slices)
                        str = char(slices(k));
                        ind = findstr(str,'km');
                        depth(k) = str2double(str(ind-7:ind-1));
                        f = fopen(str,'r');
                        textscan(f,'# %f %f',1);
                        temp = textscan(f,'%f %f %f');
                        fclose(f);
                        value(:,:,k) = reshape(temp{3},[nlon nlat]).';
                    end
                    % now fill in the grid details from file
                    lims = llgrid.limits(latitude,longitude,depth);
                    obj = llmodel(llgrid(lims));
                    obj.location = location;
                    if strfind(obj.location,'CELL')
                        [obj.grid.dlat,obj.grid.lat] = latlontools.delta(latitude(:,1));
                        [obj.grid.dlon,obj.grid.lon] = latlontools.delta(longitude(1,:));
                        [obj.grid.dz,obj.grid.depth] = latlontools.delta(depth);
                    else
                        obj.grid.lat = latitude(:,1);
                        obj.grid.lon = longitude(1,:);
                        obj.grid.depth = depth;
                        obj.grid.dlat = diff(obj.grid.lat);
                        obj.grid.dlon = diff(obj.grid.lon);
                        obj.grid.dz = diff(obj.grid.depth);
                    end
                    obj.grid.nlat = length(obj.grid.dlat);
                    obj.grid.nlon = length(obj.grid.dlon);
                    obj.grid.nzEarth = length(obj.grid.dz);
                    obj.v = value;
                case 'netcdf'
                    % open NetCDF file and initialize the object with
                    % grid limits
                    ncid = netcdfiris.open(cfile);
                    try
                        lims = netcdfiris.getLimits(ncid);
                    catch exception
                        if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                            lims = [];
                        end
                    end
                    try
                        [header,geospatial] = netcdfiris.getHeader(ncid);
                    catch exception
                        if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                            header = [];
                            geospatial = [];
                        end
                    end
                    try
                        modelvar = netcdfiris.getModelVariables(ncid);
                    catch exception
                        if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                            modelvar = [];
                        end                        
                    end
                    [longitude,latitude,depth] = netcdfiris.getPoints(ncid);
                    % now keep reading the variables until 
                    [vars,myind] = netcdfiris.getValue(ncid);
                    netcdf.close(ncid);
                    if isglobal && longitude(end)<360
                        longitude(end+1)=360;
                        for i=1:length(vars)
                            vars(i).value(end+1,:,:) = vars(i).value(1,:,:);
                        end
                    end
                    % now fill in the grid details from file - the values
                    % are specified at points which is, for our purposes,
                    % equivalent to cell centers - and that's the way we
                    % want it for any future use (plotting or modeling).
                    obj = llmodel(llgrid);
                    obj.location = 'CELL';
                    obj.grid.nlat = length(latitude);
                    obj.grid.nlon = length(longitude);
                    obj.grid.nzEarth = length(depth);
                    obj.grid.dlat = zeros(obj.grid.nlat,1);
                    obj.grid.dlon = zeros(obj.grid.nlon,1);
                    obj.grid.dz = zeros(obj.grid.nzEarth,1);                 
                    obj.grid.lat = zeros(obj.grid.nlat+1,1);
                    obj.grid.lon = zeros(obj.grid.nlon+1,1);
                    obj.grid.depth = zeros(obj.grid.nzEarth+1,1);                 
                    if isempty(modelvar)                    
                        % we make an educated guess at the edges assuming
                        % equal spacing in the center
                        [obj.grid.dlat,obj.grid.lat] = latlontools.delta(latitude);
                        [obj.grid.dlon,obj.grid.lon] = latlontools.delta(longitude);
                        [obj.grid.dz,obj.grid.depth] = latlontools.deltaz(depth);
                        obj.lat0 = obj.grid.lat(1);
                        obj.lon0 = obj.grid.lon(1);
                        % also, set up the modelVariables
                        % add global variable primary_coords [latlon|xy]
                        modelvar.primary_coords = 'latlon';
                        
                        % use corner coordinates to exactly compute grid cells from grid centers
                        modelvar.corner_description = 'used to compute exact grid geometry from cell centers';
                        modelvar.corner_location = 'upper_southwest';
                        modelvar.corner_latitude = obj.grid.lat(1);
                        modelvar.corner_longitude = obj.grid.lon(1);
                        modelvar.corner_depth = obj.grid.depth(1);
                        
                        % model rotation complements data rotation information from data file
                        modelvar.rotation_units = 'degrees';
                        modelvar.rotation_angle = 0.0;
                    else
                        % use the provided modelVariables to compute
                        % accurate distances on the grid; rotation is
                        % recorded but currently not used to compute the
                        % grid corner coordinates
                        obj.grid.lat(1) = modelvar.corner_latitude;
                        for i=1:obj.grid.nlat
                            obj.grid.dlat(i) = 2*(latitude(i)-obj.grid.lat(i));
                            obj.grid.lat(i+1) = obj.grid.lat(i)+obj.grid.dlat(i);
                        end
                        obj.grid.lon(1) = modelvar.corner_longitude;
                        for i=1:obj.grid.nlon
                            obj.grid.dlon(i) = 2*(longitude(i)-obj.grid.lon(i));
                            obj.grid.lon(i+1) = obj.grid.lon(i)+obj.grid.dlon(i);
                        end
                        obj.grid.depth(1) = modelvar.corner_depth;
                        for i=1:obj.grid.nzEarth
                            obj.grid.dz(i) = 2*(depth(i)-obj.grid.depth(i));
                            obj.grid.depth(i+1) = obj.grid.depth(i)+obj.grid.dz(i);
                        end
                        obj.lat0 = modelvar.corner_latitude;
                        obj.lon0 = modelvar.corner_longitude;
                        % Rotation not currently supported, will do later
                        %obj.rotation = modelvar.rotation_angle;
                    end
                    % if no limits in file, make them now
                    if isempty(lims)
                        lims.lonmin = min(obj.grid.lon);
                        lims.lonmax = max(obj.grid.lon);
                        lims.latmin = min(obj.grid.lat);
                        lims.latmax = max(obj.grid.lat);
                        lims.depthmin = min(obj.grid.depth);
                        lims.depthmax = max(obj.grid.depth);
                    end
                    obj.grid.limits = lims;
                    % and the model values as X(latitude,longitude,depth)
                    disp(['Reading the variable ' vars(myind).short_name ' from NetCDF file into a model object']);
                    obj.modelName = vars(myind).short_name;
                    obj.modelType = vars(myind).long_name;
                    obj.modelUnits = vars(myind).units;
                    obj.displayName = vars(myind).display_name;
                    if contains(vars(myind).display_name,'log(10)') || contains(vars(myind).display_name,'log10')
                        obj.paramType = 'LOG10';
                    elseif contains(vars(myind).display_name,'ln')
                        obj.paramType = 'LOGE';
                    else
                        obj.paramType = 'LINEAR';
                    end
                    obj.v = permute(vars(myind).value,[2 1 3]);
                    disp('Output permuted to conform to X(latitude,longitude,depth)');
                    % finally, save metadata
                    obj.fileHeader = header;
                    obj.modelVariables = modelvar;
                    obj.geospatialInfo = geospatial;
                case 'modem'
                    [dlat,dlon,dz,rho,nzAir,type,origin,rotation] = read_WS3d_model(cfile);
                    obj = llmodel;
                    obj.location = 'CELL';
                    obj.paramType = type;
                    switch type
                        case {'LOGE','LOG10'}
                            obj.v = - rho;
                        otherwise
                            obj.v = 1./(rho);
                    end
                    switch type
                        case 'LOGE'
                            obj.AirCond = log(obj.AirCond);
                            obj.SeaWaterCond = log(obj.SeaWaterCond);
                        case 'LOG10'
                            obj.AirCond = log10(obj.AirCond);
                            obj.SeaWaterCond = log10(obj.SeaWaterCond);
                    end
                    lowerleftlat = origin(1);
                    lowerleftlon = origin(2);
                    obj.grid.dlat = dlat;
                    obj.grid.dlon = dlon;
                    obj.grid.dz = dz/1000;
                    obj.grid.lat = lowerleftlat + [0; cumsum(obj.grid.dlat)];
                    obj.grid.lon = lowerleftlon + [0; cumsum(obj.grid.dlon)];
                    obj.grid.depth = [0; cumsum(obj.grid.dz)];
                    obj.grid.units = 'km';
                    obj.grid.nlat = length(obj.grid.dlat);
                    obj.grid.nlon = length(obj.grid.dlon);
                    obj.grid.nzAir = nzAir;
                    obj.grid.nzEarth = length(obj.grid.dz);
                    obj.grid.limits.latmin = min(obj.grid.lat);
                    obj.grid.limits.latmax = max(obj.grid.lat);
                    obj.grid.limits.lonmin = latlontools.lon360(min(obj.grid.lon));
                    obj.grid.limits.lonmax = latlontools.lon360(max(obj.grid.lon));
                    obj.grid.limits.depthmin = 0;
                    obj.grid.limits.depthmax = max(obj.grid.depth);
                    % finally, convert from LOGE to LOG10 at least until LOG10
                    % finally gets implemented in ModEM
                    if strfind(obj.paramType,'LOGE')
                        obj = obj.log10;
                    end
                    % leave the origin as lower left corner - all we do here is read; 
                    % change it later to grid center if needed for ll2xy conversion
                    obj.lat0 = origin(1);
                    obj.lon0 = origin(2);
                    %obj.lat0 = mean(obj.grid.lat);
                    %obj.lon0 = mean(obj.grid.lon);
                case 'global'
                    warning('Check the CELL/NODE compatibility for global models. Not tested yet.');
                    %[rho3d.v,rho3d.pr,rho3d.lon,rho3d.lat] = ReadModel(cfile,'GM',0,'2013');
                    [rho3d.v,rho3d.pr,rho3d.lon,rho3d.lat] = ReadModel(cfile,'GM',0);
                    disp(['Reading of global model from ' cfile ' completed']);
                    clong=0; res = [min(diff(rho3d.lon)) min(diff(rho3d.lat))];
                    [rho3d.v_gg,rho3d.lon_gg,rho3d.lat_gg] = ...
                        InterpGrid(rho3d.v,rho3d.lon,90-rho3d.lat,clong,res);
                    disp('3D model converted to geographic');
                    [rho3d.pr_gg] = InterpGrid(rho3d.pr,rho3d.lon,90-rho3d.lat,clong,res);
                    disp('Prior model converted to geographic');
                    rho3d.v = (rho3d.v + rho3d.pr);
                    rho3d.v_gg = (rho3d.v_gg + rho3d.pr_gg);
                otherwise
                    error(['Can''t read lat/lon model format ' format ': method unknown']);
            end
            
            obj.fileName = cfile;
            obj.limits = obj.grid.limits;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = prior(grid,priorType,rho,mask)
            %
            % [obj] = prior(grid,priorType,rho,mask)
            % creates a 3D resistivity prior model
            % rho should be linear resistivity
            % priorType may be uniform or radial
            % if radial, specify top and bottom rho values;
            % resistivity will vary logarithmically with depth.
            % optional covariance mask sets air & oceans;
            % use [COV,ELEV] = mask(grid) to compute.
            % output will be *electrical resistivity*.
            
            obj = llmodel(grid,rho(1),'LINEAR','resistivity','Ohm*m');
            
            switch priorType
                case 'uniform'
                    if length(rho)~=1 || rho <= 0
                        error('Use one linear resistivity value for a uniform prior')
                    end
                case 'radial'
                    if length(rho)~=2 || min(rho) <= 0
                        error('Use two linear resistivity values for a radial prior')
                    end
                    totaldepth = abs(max(grid.z) - min(grid.z));
                    depthfactor = (max(abs(grid.z))-abs(grid.z))/totaldepth;
                    radialrho(1:grid.nzEarth) = max(rho(1)*depthfactor(1:grid.nzEarth),rho(2));
                    for k = 1:grid.nzEarth
                        obj.v(:,:,k) = radialrho(k);
                    end
                otherwise
                    error('Prior types may be uniform or radial');
            end
            
            % convert to linear conductivity
            %obj.v = 1./(obj.v);
            %obj.modelType = 'electrical conductivity';
            %obj.modelUnits = 'S/m';
            obj.location = 'CELL';
            
            if nargin > 3
                obj.v(mask == grid.AIR) = 1/obj.AirCond;
                obj.v(mask == grid.OCEAN) = 1/obj.SeaWaterCond;
            end
        end    
        
    end

end

    
                                
   
                
                