classdef xymodel < modelplot
    %   addd depth and conductivity to the latgrid class, to store
    %   conductivity on a standard cartesian grid with lat/lon info, and
    %   tools for extracting slices in lat/lon/depth coordinates
    
%     properties
%         grid       % xygrid
%         modelType  = 'electrical conductivity'; % conductivity or resistivity
%         modelUnits = 'S/m'; % S/m or Ohm*m
%         paramType  % LOGE, LOG10 or LINEAR
%         v          % 3D array
%         AirCond = 1e-10;
%         SeaWaterCond = 4.8;
%     end

    properties
        fileName
        fileHeader
        modelVariables
    end
%     
%     properties (SetAccess = protected)
%         % a non-regular mesh of LAT/LON values computed using a lat0/lon0
%         % origin in geographic coords and used for interpolation to a
%         % regular lat/lon grid (class llgrid object)
%         %LAT
%         %LON
%         %Z
%         lat0
%         lon0
%         limits
%     end

    
    methods
        function [obj] = xymodel(varargin)
            %   class constructor   
            obj = obj@modelplot(varargin{:});
            
            if nargin == 0
                obj.grid = xygrid;
                return
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [status] = write(obj,cfile,format)
            %
            % Usage:  [status] = write(obj,cfile,format)
            %
            %  writes xymodel object;
            %  status is total number of bytes written
            %
            %  format = 'RM' means Randie Mackie's format
            %  format = 'WS' means Weerachai Siripunvaraporn's format
            
            nzAir = obj.grid.nzAir;
            dx = obj.grid.dx;
            dy = obj.grid.dy;
            dz = obj.grid.dz;
            
            % first convert from LOG10 to LOGE at least until LOG10
            % finally gets implemented in ModEM
            %if strfind(obj.paramType,'LOG10')
            %    obj = obj.loge;
            %end
            
            if isfield(struct(obj.grid),'origin')
                origin = obj.grid.origin;
            else
                origin = [0 0 0];
            end
            
            if isfield(struct(obj.grid),'units')
                if strcmp(obj.grid.units,'km')
                    % convert everything to meters!
                    dx = 1000*dx;
                    dy = 1000*dy;
                    dz = 1000*dz;
                    origin = 1000*origin;
                end
            end
            
            if isfield(struct(obj.grid),'rotation')
                rotation = obj.grid.rotation;
            else
                rotation = 0;
            end
            
            if isempty(obj.modelVariables)
                % add global variable primary_coords [latlon|xy]
                modelvar.primary_coords = 'xy';
                
                % origin always defined in meters from upper southwest corner
                modelvar.origin_description = 'defined in meters from upper southwest corner';
                modelvar.origin_location = 'data_zero';
                modelvar.origin_x = origin(1);
                modelvar.origin_y = origin(2);
                modelvar.origin_z = origin(3);
                
                % model rotation complements data rotation information from data file
                modelvar.rotation_units = 'degrees';
                modelvar.rotation_angle = rotation;
            else
                modelvar = obj.modelVariables;
            end
            
            switch format
                case 'WS'
                    if contains(obj.modelType,'conductivity')
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
                    status = write_WS3d_model(cfile,dx,dy,dz,rho,nzAir,type,origin,rotation);
                case 'RM'
                    if contains(obj.modelType,'conductivity')
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
                    status = write_mackie3d_model(cfile,dx,dy,dz,rho,nzAir,type,origin,rotation);
                case 'netcdf'
                    if isempty(obj.modelVariables)
                        error('Unable to write the model to NetCDF: please provide obj.modelVariables needed for forward modeling.');
                    end
                    header = obj.fileHeader;
                    if ~strfind(obj.location,'CELL')
                        error('NetCDF model format assumes models defined in CELLS. Run obj.node2cell or expand the format.');
                    else
                        latitude = obj.grid.xctr;
                        longitude = obj.grid.yctr;
                        depth = obj.grid.zctr;
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
                    vars(1).value = permute(obj.v,[2 1 3]);
                    disp('Values permuted to conform to X(longitude,latitude,depth)');
                    % now, open NetCDF file and write - use a new name
                    ncid = netcdfxyz.open(cfile);
                    netcdfxyz.putHeader(ncid,header);
                    netcdfxyz.putModelVariables(ncid,modelvar)
                    netcdfxyz.putPoints(ncid,longitude,latitude,depth);
                    netcdfxyz.putValue(ncid,vars);
                    netcdf.close(ncid);
                    disp(['Model object written to NetCDF file ' cfile]);
                otherwise
                    error(['Can''t write ' obj.modelType ' model: unknown file format ' format]);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = uiplot(obj,padding)
            
            Nx = obj.grid.Nx;
            Ny = obj.grid.Ny;
            Nz = obj.grid.NzEarth;
            obj.v = reshape(obj.v,[Nx Ny Nz]);
            if nargin > 1
                options.padding = padding;
            else
                options.padding = 0;
            end
            if contains(obj.modelType,'conductivity')
                options.cblabel = '\sigma';
                options.clims = [-4 0];
            else
                options.cblabel = '\rho';
                options.clims = [0 4];
            end
            % just in case, convert to log10
            pltobj = obj.log10;
            options.slice = 'Z';
            options.Np = 1;
            options.iXlim(1) = 1;
            options.iXlim(2) = Nx+1;
            options.iYlim(1) = 1;
            options.iYlim(2) = Ny+1;
            options.iZlim(1) = 1;
            options.iZlim(2) = Nz+1;
            CondPlotSet(pltobj.v,pltobj.grid,options);
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function llobj = llmodel(obj,mstruct,newgrid,varargin)
            % Usage:
            %    llobj = llmodel(obj,mstruct,newgrid,varargin)
            %
            % convert model parameter to regular lat/lon grid
            %   defaults: largest set of lats and lons possible
            %   converts to a lat/lon model defined in CELLS
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
            
            % if mstruct is specified, it overwrites the origin
            if nargin > 1
                lat0 = mstruct.origin(1);
                lon0 = mstruct.origin(2);
            end
            
            % make a default output grid if it is not provided
            if nargin > 2 && isobject(newgrid)
                latgrid = newgrid;
            else
                latgrid = llgrid(obj.grid,mstruct);
            end
            inUnits = latgrid.units;
                                    
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
                        case 'depths' % in meters!!!
                            depth = varargin{k+1};
                            latgrid.depth = depth;
                            latgrid.nzEarth = length(depth)-1;
                        case 'latitudes'
                            lat = varargin{k+1};
                            latgrid.lat = lat;
                            latgrid.nlat = length(lat)-1;
                        case 'longitudes'
                            lon = varargin{k+1};
                            latgrid.lon = llgrid.lon180(lon);
                            latgrid.nlon = length(lon)-1;
                        otherwise
                            error('Optional argument not defined')
                    end
                end
            end
            
            if strcmp(method,'nointerp')
                % imitate xy2latlon function, but for any projection: 
                % no interpolation of the values, simply reshape the grid itself
                % works if the grid sizes are compatible, no other checks
                if (latgrid.nlat == obj.grid.nx) && (latgrid.nlon == obj.grid.ny) && (latgrid.nz == obj.grid.nz); then
                    llobj = llmodel(latgrid,obj.v,obj.paramType,obj.modelType,obj.modelUnits);               
                    llobj = setOrigin(llobj,lat0,lon0);                   
                else                    
                    error('Unable to convert xymodel to llmodel without interpolation: grid sizes do not match');
                end
                
            else
                
                % make sure that everything is in meters for map projection
                latgrid = km2m(latgrid);
                inGrid = km2m(obj.grid);
                
                % ALWAYS convert to linear for interpolation ... otherwise,
                % huge errors at discontinuities such as coast but only if
                % linear interpolation is used
                paramType = obj.paramType;
                if ~strcmp(paramType,'LINEAR')
                    obj = linear(obj);
                end
                
                % initialize llobj
                llobj = llmodel;
                
                % make a regular output mesh at cell centers; use map
                % projection to make an irregular x/y grid
                [LON,LAT] = meshgrid(latgrid.lon(1:end-1)+latgrid.dlon/2,latgrid.lat(1:end-1)+latgrid.dlat/2);
                nlon = length(latgrid.dlon);
                nlat = length(latgrid.dlat);
                ndepth = length(latgrid.dz);
                nTot = nlon*nlat*ndepth;
                temp = reshape(latgrid.depth(1:end-1)+latgrid.dz/2,1,ndepth);
                ll = [reshape(LAT,[1,nlon*nlat]);reshape(LON,[1,nlon*nlat])];
                if strcmp(mstruct.mapprojection,'xy2latlon')
                    lat0 = mstruct.origin(1);
                    lon0 = mstruct.origin(2);
                    xy = llgrid.ll2xy(ll,lat0,lon0,0);
                else
                    try
                        mstruct = defaultm(mstruct);
                    catch
                        error(['Unknown map projection: ',mstruct.mapprojection]);
                    end
                    switch strtrim(mstruct.mapprojection)
                        case 'utm'
                            [E,N] = mfwdtran(mstruct,LAT,LON);
                            xy(1,:) = reshape(N,1,nlon*nlat);
                            xy(2,:) = reshape(E,1,nlon*nlat);
                        otherwise % e.g., 'eqdcylin','eqacylin','lambertstd','eqaazim'
                            [E,N] = projfwd(mstruct,LAT,LON);
                            xy(1,:) = reshape(N,1,nlon*nlat);
                            xy(2,:) = reshape(E,1,nlon*nlat);
                    end
                end
                XI = reshape(xy(1,:)'*ones(1,ndepth),nTot,1);
                YI = reshape(xy(2,:)'*ones(1,ndepth),nTot,1);
                ZI = reshape(ones(nlat*nlon,1)*temp,nTot,1);
                
                % set up a regular x/y mesh
                [Y,X,Z] = meshgrid(inGrid.yctr,inGrid.xctr,inGrid.zctr);
                
                % interpolate from regular to irregular x/y mesh
                llobj.v = interp3(Y,X,Z,obj.v,YI,XI,ZI,method);
                llobj.v = reshape(llobj.v,[nlat,nlon,ndepth]);
                
                % put together the complete llobj
                llobj.grid  = latgrid;
                llobj.location = 'CELL';
                llobj.modelType  = obj.modelType;
                llobj.modelUnits = obj.modelUnits;
                llobj.paramType  = obj.paramType;
                llobj.displayName = obj.displayName;
                llobj.AirCond  = obj.AirCond;
                llobj.SeaWaterCond  = obj.SeaWaterCond;
                llobj = setOrigin(llobj,lat0,lon0);
                
                % NOW convert llobj back to log if needed
                if strcmp(paramType,'LOG10')
                    llobj = log10(llobj);
                elseif strcmp(paramType,'LOGE')
                    llobj = loge(llobj);
                end
                
                % AND convert vertical distances back to km if needed
                if contains(inUnits,'km')
                    llobj.grid = m2km(llobj.grid);
                end
            
            end
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [lims,LAT,LON] = latlonlims(obj,mstruct)
            % Usage:
            %    [lims,LAT,LON] = latlonlims(obj,mstruct)
            
            inGrid = km2m(obj.grid);
            [Y,X] = meshgrid(inGrid.yctr,inGrid.xctr);
            
            switch strtrim(mstruct.mapprojection)
                case 'xy2latlon'
                    lat0 = mstruct.origin(1);
                    lon0 = mstruct.origin(2);
                    xy = [reshape(X,[1,nlon*nlat]);reshape(Y,[1,nlon*nlat])];
                    ll = llgrid.xy2ll(xy,lat0,lon0);
                    LAT = reshape(ll(1,:),nlon,nlat);
                    LON = reshape(ll(2,:),nlon,nlat);
                case 'utm'
                    [LAT,LON] = minvtran(mstruct,Y,X);
                otherwise
                    [LAT,LON] = projinv(mstruct,Y,X);
            end
            
            lims.latmin = min(min(LAT));
            lims.latmax = max(max(LAT));
            lims.lonmin = min(min(LON));
            lims.lonmax = max(max(LON));
            lims.depthmin = min(obj.grid.z);
            lims.depthmax = max(obj.grid.z);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = setModelVariables(obj,originType,origin,rotation)
            
            % obj = setModelVariables(obj,originType,origin,rotation)
            %
            % sets the model variables needed for ModEM to function
            
            if nargin<2
                originType = 'model_center_ground_level';
            elseif isstruct(originType)
                obj.modelVariables = originType;
                return
            end
            
            if nargin<3
                origin = obj.grid.origin;
            end
            
            if nargin<4
                rotation = 0.0;
            end
            
            % add global variable primary_coords [latlon|xy]
            obj.modelVariables.primary_coords = 'xy';
            
            % origin always defined in meters from upper southwest corner
            obj.modelVariables.origin_description = 'defined in meters from upper southwest corner';
            obj.modelVariables.origin_location = originType;
            obj.modelVariables.origin_x = origin(1);
            obj.modelVariables.origin_y = origin(2);
            obj.modelVariables.origin_z = origin(3);
            
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
                        
            [newgrid,ii,jj,kk] = trim(obj.grid,nTrim,direction);
            
            obj.v = obj.v(ii,jj,kk);
            obj.grid = newgrid;
            
            % update model variables
            obj = obj.setModelVariables;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [h,X] = plot(obj,type,value,P)
            % [h,X] = plot(obj,type,value,P)
            %
            % Sample usage: h = obj.plot('depth',16)
            % Optional structure P contains plotting parameters.
            % By default does 1D interpolation which can be switched off
            % by setting P.nointerpbyvalue = 1.
            
            if strcmp(obj.location,'CELL')
                y = obj.grid.yctr;
                x = obj.grid.xctr;
                depth = [0; cumsum(obj.grid.dz(1:end-1))] + obj.grid.dz/2;
            elseif strcmp(obj.location,'NODE')
                y = obj.grid.y;
                x = obj.grid.x;
                depth = [0; cumsum(obj.grid.dz(1:end))];
            else
                error('Set model location to CELL or NODE before plotting.')
            end
                    
            if strcmp(type,'depth')
                k2 = find(depth>value,1);
                k = max(k2-1,1);
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
                    disp(['Interpolated to depth ' mydepth ' km from depths '...
                        num2str(depth(k)) ' and ' num2str(depth(k2)) ' km']);
                else
                    mydepth = num2str(depth(k));
                    disp(['Plotting the value for depth ' mydepth ' km']);
                end
           elseif strcmp(type,'layer')
                k = value;
                X = obj.v(:,:,k);
                depth1 = 0;
                if k>1
                    depth1 = sum(obj.grid.dz(1:k-1));
                end
                depth2 = sum(obj.grid.dz(1:k));
                disp(['Plotting the value for layer between depth ' num2str(depth1) ' and ' num2str(depth2) ' km']);
            end

            delta = 1e-6;
            clims = [nanmin(nanmin(X))-delta nanmax(nanmax(X))+delta];
            cmap = 'jet';
            nointerp = 0;
            showmesh = 0;
            if nargin >= 5
                if isfield(P,'nointerp')
                    nointerp = P.nointerp;
                end
                if isfield(P,'showmesh')
                    showmesh = P.showmesh;
                end
            end
        
            posn = [1,1,12,7];
            h = figure('Position',100*posn,...
                'PaperPosition',posn,...
                'PaperOrientation','Portrait',...
                'Color',[1 1 1]);

            imagesc(y,x,X); hold on; grid on; 
            h.CurrentAxes.FontWeight = 'demi';
            h.CurrentAxes.FontSize = 20;

            caxis(clims); colormap(cmap); colorbar
            if ~nointerp; shading flat; end
            if showmesh; shading faceted; end


        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = read(cfile,format)
            %
            % Usage:  [cond] = readCond_3D(cfile,format)
            %
            %  Reads in xymodel object as electrical conductivity
            %
            %  format = 'RM' means Randie Mackie's format
            %  format = 'WS' means Weerachai Siripunvaraporn's format
            
            if ~exist(cfile,'file')
                error(['Model file ' cfile ' not found']);
            end
            
            switch format
                case 'WS'
                    [dx,dy,dz,rho,nzAir,type,origin,rotation] = read_WS3d_model(cfile);
                    % convert from resistivity to conductivity
                    switch type
                        case {'LOGE','LOG10'}
                            cond = - rho;
                        case 'LINEAR'
                            cond = 1./(rho);
                        otherwise
                            error(['Unknown conductivity model type ' type]);
                    end
                case 'RM'
                    [dx,dy,dz,rho,nzAir,type,origin,rotation] = read_mackie3d_model(cfile);
                    % convert from resistivity to conductivity
                    switch type
                        case {'LOGE','LOG10'}
                            cond = - rho;
                        case 'LINEAR'
                            cond = 1./(rho);
                        otherwise
                            error(['Unknown conductivity model type ' type]);
                    end
                case 'netcdf'
                    % open NetCDF file and initialize the object with
                    % grid limits
                    ncid = netcdfxyz.open(cfile);
                    try
                        lims = netcdfxyz.getLimits(ncid);
                    catch exception
                        if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                            lims = [];
                        end
                    end
                    try
                        [header,geospatial] = netcdfxyz.getHeader(ncid);
                    catch exception
                        if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                            header = [];
                            geospatial = [];
                        end
                    end
                    % now fill in the grid details from file - the values
                    % are specified at points which is, for our purposes,
                    % equivalent to cell centers - and that's the way we
                    % want it for any future use (plotting or modeling).
                    [yctr,xctr,zctr] = netcdfxyz.getPoints(ncid);
                    try
                        modelvar = netcdfxyz.getModelVariables(ncid);
                        
                        % now use the appropriate origin to accurately compute
                        % grid distances
                        origin = [- modelvar.origin_x; - modelvar.origin_y; - modelvar.origin_z];
                        
                        x = origin(1);
                        dx = zeros(size(xctr));
                        for i = 1:length(xctr)
                            dx(i) = 2*(xctr(i)-x);
                            x = x + dx(i);
                        end
                        y = origin(2);
                        dy = zeros(size(yctr));
                        for i = 1:length(yctr)
                            dy(i) = 2*(yctr(i)-y);
                            y = y + dy(i);
                        end
                        z = origin(3);
                        dz = zeros(size(zctr));
                        for i = 1:length(zctr)
                            dz(i) = 2*(zctr(i)-z);
                            z = z + dz(i);
                        end
                    catch exception
                        if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                            % we make an educated guess at the edges assuming
                            % equal spacing in the center.
                            dx = latlontools.delta(xctr);
                            dy = latlontools.delta(yctr);
                            dz = latlontools.deltaz(zctr);
                            origin = [- sum(dx)/2; - sum(dy)/2; 0.0];
                            
                            % add global variable primary_coords [latlon|xy]
                            modelvar.primary_coords = 'xy';
                            
                            % origin always defined in meters from upper southwest corner
                            modelvar.origin_description = 'defined in meters from upper southwest corner';
                            modelvar.origin_location = 'model_center_ground_level';
                            modelvar.origin_x = - origin(1);
                            modelvar.origin_y = - origin(2);
                            modelvar.origin_z = - origin(3);
                            
                            % model rotation complements data rotation information from data file
                            modelvar.rotation_units = 'degrees';
                            modelvar.rotation_angle = 0.0;
                        end
                    end
                    % now keep reading the variables until 
                    [vars,myind] = netcdfxyz.getValue(ncid);
                    netcdf.close(ncid);

                    % and the model values as X(latitude,longitude,depth)
                    disp(['Reading the variable ' vars(myind).short_name ' from NetCDF file into a model object']);
                    if contains(vars(myind).display_name,'log(10)') || contains(vars(myind).display_name,'log10')
                        type = 'LOG10';
                    elseif contains(vars(myind).display_name,'ln')
                        type = 'LOGE';
                    else
                        type = 'LINEAR';
                    end
                    rotation = 0.0;
                    nzAir = 0;
                    cond = permute(vars(myind).value,[2 1 3]);
                    disp('Output permuted to conform to X(latitude,longitude,depth)');
                otherwise
                    error(['Can''t read ' cfile ': unknown file format ' format]);
            end
            
            grid = xygrid(dx/1000,dy/1000,dz/1000,origin/1000,rotation,'km',nzAir);
                
            obj = xymodel(grid,cond,type);
            
            if exist('modelvar','var')
                obj = obj.setModelVariables(modelvar);
            else
                obj = obj.setModelVariables('data_zero',origin,rotation);
            end
            
            % finally, convert from LOGE to LOG10 at least until LOG10
            % finally gets implemented in ModEM
            if strfind(obj.paramType,'LOGE')
                obj = obj.log10;
            end
            
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
            
            obj = xymodel(grid,rho(1),'LINEAR','resistivity','Ohm*m');
            
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
                    %for i = 1:length(SRPYgrid.depth)-1
                    %    sigma(:,:,i) = log10(1/layerrho);
                    %    layerrho = layerrho/1.05;
                    %end
                otherwise
                    error('Prior types may be uniform or radial');
            end
            
            % convert to linear conductivity
            %obj.v = 1./(obj.v);
            %obj.modelType = 'conductivity';
            %obj.modelUnits = 'S/m';
            
            if nargin > 3
                obj.v(mask == grid.AIR) = 1/obj.AirCond;
                obj.v(mask == grid.OCEAN) = 1/obj.SeaWaterCond;
            end
      
        end
    end

end

    
                                
   
                
                