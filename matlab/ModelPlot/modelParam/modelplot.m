classdef modelplot < latlontools & geoplotaddons
    %   slice and plot xy and ll models
    %   (c) Anna Kelbert, 2014-2016
    
    properties
        grid       % xy or ll grid
        limits     % used for plotting, to trim models etc overrides grid.limits
        isglobal = 0 % logical
        location   = 'CELL'; % vs 'NODE'
        modelName  = '\sigma';
        modelType  = 'electrical conductivity'; % conductivity or resistivity
        modelUnits = 'S/m'; % S/m or Ohm*m
        displayName % for plotting
        paramType  % LOGE, LOG10 or LINEAR
        v          % 3D array
        v_ref      % the prior or reference useful for plotting variations
        AirCond = 1e-10; % this is ALWAYS conductivity; could be log/log10
        SeaWaterCond = 4.8; % this is ALWAYS conductivity; could be log/log10
    end
    
    properties (SetAccess = protected)
        % defines the lat/lon origin and limits
        lat0
        lon0
    end
    
    methods
        function [obj] = modelplot(varargin)
            %   class constructor
            %   [obj] = modelplot(grid,value,paramType,modelType,modelUnits)
            %     OR
            %   [obj] = modelplot(grid,location,paramType,modelType,modelUnits)
            %         
            %   if no value, initialized with NaN's.
            %   example call to initialize without a value:
            %       obj = modelplot(grid,'CELL');
            %   grid always defined at cell centers
            if nargin == 0
                return
            end
            
            if nargin>0
                obj.grid = varargin{1};
                if isa(obj.grid,'xygrid')
                    nx = length(obj.grid.dx);
                    ny = length(obj.grid.dy);
                    nz = length(obj.grid.dz);
                elseif isa(obj.grid,'llgrid')
                    nx = length(obj.grid.dlat);
                    ny = length(obj.grid.dlon);
                    nz = length(obj.grid.dz);
                else
                    error('Unknown grid object');
                end
            end
            
            valueSupplied = 0;
            if nargin == 1
                obj.location = 'CELL';
            elseif ischar(varargin{2})
                obj.location = varargin{2};
            else
                % get location based on the value size
                value = varargin{2};
                valueSupplied = 1;
                if size(value,1) == nx+1
                    obj.location = 'NODE';
                elseif size(value,1) == nx
                    obj.location = 'CELL';
                else
                    error('Grid and value size mismatch on input');
                end
            end
            
            if strcmp(obj.location,'CELL')
                obj.v(1:nx,1:ny,1:nz) = NaN;
                obj.paramType = '';
            elseif strcmp(obj.location,'NODE')
                obj.v(1:nx+1,1:ny+1,1:nz+1) = NaN;
                obj.paramType = '';
            end
            
            if valueSupplied
                obj.v = value;
            end
            
            if nargin == 2
                if ~ischar(varargin{2})
                    error('Missing model parameter type LINEAR, LOG10 or LOGE');
                else
                    return
                end
            end
            
            %to convert to log(10), if needed, use log10(obj)
            if nargin >= 3
                obj.paramType = varargin{3}; % type
            else
                obj.paramType = 'LOG10';
            end
            
            if nargin >= 5
                obj.modelType = varargin{4};
                obj.modelUnits = varargin{5};
            end
            if strfind(obj.modelType,'conductivity')
                obj.modelName = '\sigma';
            elseif strfind(obj.modelType,'resistivity')
                obj.modelName = '\rho';
            end
            
            if strcmp(obj.paramType,'LOG10')
                obj.modelName = ['log(10) (' obj.modelName ')'];
                obj.displayName = ['log(10) ' obj.modelType];
                obj.AirCond = log10(obj.AirCond);
                obj.SeaWaterCond = log10(obj.SeaWaterCond);
            elseif strcmp(obj.paramType,'LOGE')
                obj.modelName = ['ln (' obj.modelName ')'];
                obj.displayName = ['ln ' obj.modelType];
                obj.AirCond = log(obj.AirCond);
                obj.SeaWaterCond = log(obj.SeaWaterCond);
            else
                obj.displayName = [obj.modelType];
            end
            
            disp(['Initializing a model of type: ' obj.displayName ...
                ', units ' obj.modelUnits]);
       end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = fillNaNs(obj,method,bg)
            % obj = fillNaNs(obj,method,bg)
            %
            % replace all NaNs with the depth-averaged value, nearest neighbor
            % value, or the specified depth-specific bg value; comes handy
            % after an llmodel to xymodel conversion, or vice versa
            %
            % if optional background bg is specified, it is linear (not log) 
            % and is either scalar, or has the same number of values as 
            % as the number of vertical layers in the model. 
            % use bg = 0 for conductivity perturbations  
            %
            % options: 'mean', 'background' or 'nearest'; default 'nearest'
            
            if nargin < 2
                method = 'nearest';
            end
            
            nz = size(obj.v,3);
            if strcmp(method,'background')
                if nargin < 3
                    warning('Background conductivity or resistivity has not been specified. Using depth average.');
                    method = 'mean';
                else
                    if isscalar(bg)
                        bg(1:nz) = bg;
                    end
                    if strcmp(obj.paramType,'LOGE') && (min(bg) > 0)
                        bg = log(bg);
                    end
                    if strcmp(obj.paramType,'LOG10') && (min(bg) > 0)
                        bg = log10(bg);
                    end
                end
            end
            
            nantotal = length(find(isnan(obj.v)));
            if nantotal > 0
                disp(['Using method ''' method ''' to fill in ' num2str(nantotal) ' NaNs.']);
            else
                disp('No NaNs to fill in.');
            end
            for k = 1:size(obj.v,3)
                vk = squeeze(obj.v(:,:,k));
                if strcmp(method,'mean')  
                    vk(isnan(vk)) = nanmean(nanmean(vk));
                    obj.v(:,:,k) = vk;
                elseif strcmp(method,'background')
                    vk(isnan(vk)) = bg(k);
                    obj.v(:,:,k) = vk;
                elseif strcmp(method,'nearest')
                    nlon = size(obj.v,2);
                    nlat = size(obj.v,1);
                    nantotal = length(find(isnan(vk)));
                    while nantotal > 0
                        for i = 1:size(obj.v,2)
                            for j = 1:size(obj.v,1)
                                if isnan(vk(j,i))
                                    imin = max(1,i-3);
                                    imax = min(nlon,i+3);
                                    jmin = max(1,j-3);
                                    jmax = min(nlat,j+3);
                                    nearest = vk(jmin:jmax,imin:imax);
                                    obj.v(j,i,k) = nanmean(nanmean(nearest));
                                end
                            end
                        end
                        vk = squeeze(obj.v(:,:,k));
                        nantotal = length(find(isnan(vk)));
                    end
                end   
                
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = loge(obj)            
            switch obj.paramType
                case 'LOG10'
                    obj.v = log(10) * obj.v;
                    obj.AirCond = log(10) * obj.AirCond;
                    obj.SeaWaterCond = log(10)*obj.SeaWaterCond;
                case 'LINEAR'
                    obj.v = log(obj.v);
                    obj.AirCond = log(obj.AirCond);
                    obj.SeaWaterCond = log(obj.SeaWaterCond);
                otherwise
                    if ~strcmp(obj.paramType,'LOGE')
                        error(['Unknown model parameter type ' obj.paramType]);
                    end
            end
            obj.paramType = 'LOGE';
            
            if strfind(obj.modelType,'conductivity')
                obj.modelName = 'ln(\sigma)';
                obj.displayName = 'natural log of electrical conductivity';
            elseif strfind(obj.modelType,'resistivity')
                obj.modelName = 'ln(\rho)';
                obj.displayName = 'natural log of electrical resistivity';
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = log10(obj)            
            switch obj.paramType
                case 'LOGE'
                    obj.v = obj.v/log(10);
                    obj.AirCond = obj.AirCond/log(10);
                    obj.SeaWaterCond = obj.SeaWaterCond/log(10);
                case 'LINEAR'
                    obj.v = log10(obj.v);
                    obj.AirCond = log10(obj.AirCond);
                    obj.SeaWaterCond = log10(obj.SeaWaterCond);
                otherwise
                    if ~strcmp(obj.paramType,'LOG10')
                        error(['Unknown model parameter type ' obj.paramType]);
                    end
            end
            obj.paramType = 'LOG10';            

            if strfind(obj.modelType,'conductivity')
                obj.modelName = 'log10(\sigma)';
                obj.displayName = 'log(10) of electrical conductivity';
            elseif strfind(obj.modelType,'resistivity')
                obj.modelName = 'log10(\rho)';
                obj.displayName = 'log(10) of electrical resistivity';
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = linear(obj)            
            switch obj.paramType
                case 'LOGE'
                    obj.v = exp(obj.v);
                    obj.AirCond = exp(obj.AirCond);
                    obj.SeaWaterCond = exp(obj.SeaWaterCond);
                case 'LOG10'
                    obj.v = 10.^(obj.v);
                    obj.AirCond = 10.^(obj.AirCond);
                    obj.SeaWaterCond = 10.^(obj.SeaWaterCond);
                otherwise
                    if ~strcmp(obj.paramType,'LINEAR')
                        error(['Unknown model parameter type ' obj.paramType]);
                    end
            end
            obj.paramType = 'LINEAR';            

            if strfind(obj.modelType,'conductivity')
                obj.modelName = '\sigma';
                obj.displayName = 'electrical conductivity';
            elseif strfind(obj.modelType,'resistivity')
                obj.modelName = '\rho';
                obj.displayName = 'electrical resistivity';
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = lon360(obj)
            % convert the model parameter value to longitude [0,360)
            % convention
            
            if isnumeric(obj)
                obj = llgrid.lon360(obj);
            elseif strcmp(obj.location,'CELL') 
                midlon = llgrid.lon360(obj.grid.lon(1:end-1)+obj.grid.dlon/2);
                [midlon,ii] = sort(midlon);
                obj.v = obj.v(:,ii,:);
                obj.grid.lon = sort(llgrid.lon360(obj.grid.lon));
                obj.grid.dlon = diff(obj.grid.lon);
                % unreasonable inputs with even a single outlier produce a
                % completely broken set of values - do not use this
                %[obj.grid.dlon,obj.grid.lon] = obj.delta(midlon); 
                if obj.isglobal
                    % deal separately with global model zero longitude problem
                    % which otherwise leaves us with repeated values
                    jj = find(obj.grid.dlon > 0);
                    obj.grid.dlon = obj.grid.dlon(jj);
                    obj.grid.lon = obj.grid.lon(1) + [0; cumsum(obj.grid.dlon)];
                    obj.v = obj.v(:,jj,:);
                    % now pad with a new zero longitude
                    obj.grid.lon(end+1) = obj.grid.lon(1)+360;
                    obj.grid.dlon(end+1) = obj.grid.dlon(1);
                    obj.v(:,end+1,:) = obj.v(:,1,:);
                    % now fix the grid and model limits
                    obj.grid.limits.lonmin = 0;
                    obj.grid.limits.lonmax = 360;                
                    obj.limits.lonmin = 0;
                    obj.limits.lonmax = 360;   
                else         	
                    % now fix the grid and model limits
                    obj.grid.limits.lonmin = llgrid.lon360(obj.grid.limits.lonmin);
                    obj.grid.limits.lonmax = llgrid.lon360(obj.grid.limits.lonmax);                
                    obj.limits.lonmin = llgrid.lon360(obj.limits.lonmin);
                    obj.limits.lonmax = llgrid.lon360(obj.limits.lonmax); 
                end
            else
                warning('Unable to convert this model type to lon360');
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = lon180(obj)
            % convert the model parameter value to longitude (-180,180]
            % convention
           
            if isnumeric(obj)
                obj = llgrid.lon180(obj);
            elseif strcmp(obj.location,'CELL') 
                midlon = llgrid.lon180(obj.grid.lon(1:end-1)+obj.grid.dlon/2);
                [midlon,ii] = sort(midlon);
                obj.v = obj.v(:,ii,:);
                obj.grid.lon = sort(llgrid.lon180(obj.grid.lon));
                obj.grid.dlon = diff(obj.grid.lon);
                % unreasonable inputs with even a single outlier produce a
                % completely broken set of values - do not use this
                %[obj.grid.dlon,obj.grid.lon] = obj.delta(midlon);
                if obj.isglobal
                    % deal separately with global model zero longitude problem
                    % which otherwise leaves us with repeated values
                    jj = find(obj.grid.dlon > 0);
                    obj.grid.dlon = obj.grid.dlon(jj);
                    obj.grid.lon = obj.grid.lon(1) + [0; cumsum(obj.grid.dlon)];
                    obj.v = obj.v(:,jj,:);
                    % now pad with a new zero longitude - pad before
                    obj.grid.lon = [obj.grid.lon(1)-obj.grid.dlon(end); obj.grid.lon];
                    obj.grid.dlon = [obj.grid.dlon(end); obj.grid.dlon];
                    tmp = obj.v;
                    obj.v(:,1,:) = obj.v(:,end,:);
                    obj.v(:,2:end+1,:) = tmp;
                    % ... and pad after
                    %obj.grid.lon(end+1) = obj.grid.lon(end)+obj.grid.dlon(1);
                    %obj.grid.dlon(end+1) = obj.grid.dlon(1);
                    %obj.v(:,end+1,:) = obj.v(:,1,:);
                    % now fix the grid and model limits
                    obj.grid.limits.lonmin = -180;
                    obj.grid.limits.lonmax = 180;                
                    obj.limits.lonmin = -180;
                    obj.limits.lonmax = 180;   
                else
                    obj.grid.limits.lonmin = llgrid.lon180(obj.grid.limits.lonmin);
                    obj.grid.limits.lonmax = llgrid.lon180(obj.grid.limits.lonmax);                
                    obj.limits.lonmin = llgrid.lon180(obj.limits.lonmin);
                    obj.limits.lonmax = llgrid.lon180(obj.limits.lonmax);   
                end
            else
                warning('Unable to convert this model type to lon180');
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setOrigin(obj,lat0,lon0)
            % define the geographic origin and corresponding irregular 
            % spherical grid coordinates LAT & LON
            
            if nargin == 3
                obj.lat0 = lat0;
                obj.lon0 = lon0;
            elseif nargin == 2 && ischar(lat0)
                if strcmp(lat0,'lowerleft')
                    obj.lat0 = obj.grid.lat(1);
                    obj.lon0 = obj.grid.lon(1);
                elseif strcmp(lat0,'center')
                    obj.lat0 = mean(obj.grid.lat);
                    obj.lon0 = mean(obj.grid.lon);
                end
            end     
            
        end     
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setLimits(obj,lims)
            %   obj = setLimits(obj,lims)
            %   Set the lat, lon, depth limits structure to be used for
            %   plotting, trimming models, and other utility routines.
            %
            %   This used to be a protected property, obj = set.limits(obj,lims)
            %   Wasn't such a good idea and caused multiple bugs. Instead,
            %   intended use is now different: create the limits structure
            %   exactly to your liking, then save it within the object.
            %   This is merely a shorthand for this operation.
            %   lims = obj.getLimits 
            %   OR
            %   lims = latlontools.getLimits(region)
            if nargin < 2
                obj.limits = obj.getLimits;
            elseif ~isstruct(lims)
                error('Limits must be a structure that contains min and max values in lon, lat and depth')
            else
                obj.limits.latmin = lims.latmin;
                obj.limits.latmax = lims.latmax;
                obj.limits.lonmin = lims.lonmin;
                obj.limits.lonmax = lims.lonmax;
                if isfield(obj.limits,'depthmin') && isfield(obj.limits,'depthmax')
                    obj.limits.depthmin = max(lims.depthmin,obj.limits.depthmin);
                    obj.limits.depthmax = min(lims.depthmax,obj.limits.depthmax);
                else
                    obj.limits.depthmin = lims.depthmin;
                    obj.limits.depthmax = lims.depthmax;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function lims = getLimits(obj)
            %   lims = getLimits(obj)
            %   Find the lat, lon, depth limits to be used for
            %   plotting, trimming models, and other utility routines.
            %   If needed, save the output: obj.limits = lims
            %   Don't confuse with the static
            %   lims = latlontools.getLimits(region)                
            if isa(obj,'xymodel')
                if isempty(obj.lat0) || isempty(obj.lon0)
                    error('To obtain lat/lon limits, need valid lat0 and lon0');
                end
                [LAT,LON] = latlon(obj.grid,obj.lat0,obj.lon0);
                lims = struct(...
                    'latmin',min(min(LAT(:,:))),...
                    'latmax',max(max(LAT(:,:))),...
                    'lonmin',min(min(LON(:,:))),...
                    'lonmax',max(max(LON(:,:))),...
                    'depthmin',min(obj.grid.z),...
                    'depthmax',max(obj.grid.z));
            elseif isa(obj,'llmodel')
                lims = struct(...
                    'latmin',min(obj.grid.lat),...
                    'latmax',max(obj.grid.lat),...
                    'lonmin',min(obj.grid.lon),...
                    'lonmax',max(obj.grid.lon),...
                    'depthmin',min(obj.grid.depth),...
                    'depthmax',max(obj.grid.depth));
            else
                error('Unknown model object');
            end
  
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function result = zSlice(obj,depths,longitudes,latitudes)
            %   make corners for slice at a given depth; if optional
            %   longitudes/latitudes arguments are missing use all
            lims = obj.limits;
            if nargin < 2
                depths = [lims.depthmin; obj.z(obj.z > lims.depthmin & obj.z <= lims.depthmax)];
            end
            nsections = length(depths);
            if nargin < 3
                longitudes = llgrid.lon180([ lims.lonmin, lims.lonmax]);
            end
            if nargin < 4
                latitudes = [ lims.latmin, lims.latmax];
            end
            if isa(obj,'xymodel')
                latgrid = llgrid(obj.grid,obj.lat0,obj.lon0);
            elseif isa(obj,'llmodel')
                latgrid = obj.grid;
            else
                error('Unknown model object');
            end
            Corners = zeros(3,2,2);
            %  latitude
            Corners(1,1,:) = latitudes(1);
            Corners(1,2,:) = latitudes(2);
            NM(1) = floor(abs(latitudes(2)-latitudes(1))/min(latgrid.dlat));
            %  longitude
            Corners(2,:,1) = longitudes(1);
            Corners(2,:,2) = longitudes(2);
            NM(2) = floor(abs(longitudes(2)-longitudes(1))/min(latgrid.dlon));
            if nsections >= 1
                result = cell(nsections,1);
                for k = 1:nsections
                    Corners(3,:,:) = depths(k);
                    result{k} = struct('Corners',Corners,'NM',NM);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function result = latSlice(obj,latitudes,depths)
            %   make corners for slice at given latitude; if optional
            %   depths argument is missing use all depths
            nsections = length(latitudes);
            lims = obj.limits;
            if nargin == 2
                depths = [ lims.depthmin, lims.depthmax];
            end
            Corners = zeros(3,2,2);
            %   depths
            Corners(3,1,:) = depths(1);
            Corners(3,2,:) = depths(2);
            %  longitude
            Corners(2,:,1) = lims.lonmin;
            Corners(2,:,2) = lims.lonmax;
            %  set default resolution
            NM(1) = obj.grid.nlon;
            NM(2) = obj.grid.nzEarth;
            if nsections > 1
                result = cell(nsections,1);
                for k = 1:nsections
                    Corners(1,:,:) = latitudes(k);
                    result{k} = struct('Corners',Corners,'NM',NM);
                end
            else
                Corners(1,:,:) = latitudes;
                result = Corners;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Corners = lonSlice(obj,longitude,depths)
            %   make corners for slice at given latitude; if optional
            %   depths argument is missing use all depths
            lims = obj.limits;
            if nargin == 2
                depths = [ lims.depthmin, lims.depthmax];
            end
            Corners = zeros(3,2,2);
            %   depths
            Corners(3,1,:) = depths(1);
            Corners(3,2,:) = depths(2);
            %  latitude
            Corners(1,:,1) = lims.latmin;
            Corners(1,:,2) = lims.latmax;
            %  latitude
            Corners(2,:,:) = longitude;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [cond,xy,ll] = slice(obj,Corners,NM)
            %   Corners(3,2,2) = coordinates of 4 corners: longitude,
            %                                            latitude, depth
            %    nominally top/bottom and left/right
            %   NM(2) = number of divisions: top-bottom, left-right
            w1 = (0:NM(1)-1)/(NM(1)-1);
            w1 = w1'*ones(1,NM(2));
            w2 = ones(NM(1),1)*(0:NM(2)-1)/(NM(2)-1);
            ll = zeros(3,NM(1),NM(2));
            for k = 1:3
                ll(k,:,:) = Corners(k,1,1)*(1-w1).*(1-w2) + ...
                    Corners(k,1,2)*(1-w1).*w2 + ...
                    Corners(k,2,1)*w1.*(1-w2) + ...
                    Corners(k,2,2)*w1.*w2;
            end
            ll = reshape(ll,3,NM(1)*NM(2));
            ll(2,:) = llgrid.lon180(ll(2,:));
            xy = ll;
            if isa(obj,'xymodel')
                [Y,X,Z] = meshgrid(obj.grid.yctr,obj.grid.xctr,obj.grid.zctr);
                xy(1:2,:) = obj.ll2xy(ll(1:2,:),obj.lat0,obj.lon0);
                cond = interp3(Y,X,Z,obj.v,xy(2,:),xy(1,:),xy(3,:));
            elseif isa(obj,'llmodel')
                ctrlon = obj.grid.lon(1:end-1) + obj.grid.dlon/2;
                ctrlat = obj.grid.lat(1:end-1) + obj.grid.dlat/2;
                ctrz = obj.grid.depth(1:end-1) + obj.grid.dz/2;
                [LON,LAT,Z] = meshgrid(ctrlon,ctrlat,ctrz);
                cond = interp3(llgrid.lon180(LON),LAT,Z,obj.v,...
                    llgrid.lon180(ll(2,:)),ll(1,:),ll(3,:));
            else
                error('Unknown model object');
            end
            cond = reshape(cond,[NM(1),NM(2)]);
            xy = reshape(xy,[3,NM(1),NM(2)]);
            ll = reshape(ll,[3,NM(1),NM(2)]);
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fencePlot(obj,sections,varargin)
            
            figure('Position',[200,200,1400,900],...
                    'PaperPosition',[1,1,17.8,10]); %widescreen
            nslice = length(sections);
            for k = 1:nslice
                [cond,xy,ll] = slice(obj,sections{k}.Corners,sections{k}.NM); %#ok<ASGLU>
                surf(squeeze(llgrid.lon180(ll(2,:,:))),squeeze(ll(1,:,:)),squeeze(ll(3,:,:)),cond);
                shading flat
                hold on
            end
            set(gca,'Fontweight','demi','FontSize',14,'zdir','reverse')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    

end

    
                                
   
                
                