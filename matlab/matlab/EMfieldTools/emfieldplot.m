classdef emfieldplot < latlontools & geoplotaddons
    %   do coordinate conversions and plot xy and ll EM fields
    %
    %   can be initialized with ModEMM TVector3D_SG or TVector3D_Global
    
    properties
        grid       % xy or ll grid determines cartesian or latlon
        isglobal = 0 % logical
        location   = 'EDGE'; % vs 'FACE'
        emfieldName  = 'E';
        emfieldType  = 'electric field'; % electric or magnetic field
        emfieldUnits = 'V/m';
        displayName % for plotting
        x          % 3D array
        y          % 3D array
        z          % 3D array
        solver     % EM solver type
        period     % in seconds
        source     % for MT, 'X' or 'Y'
    end
    
    properties (SetAccess = protected)
        % defines the lat/lon origin and limits
        lat0
        lon0
        limits
    end
    
    methods
        function [obj] = emfieldplot(varargin)
            %   class constructor
            %   [obj] = emfieldplot(grid,vector,emfieldName,emfieldUnits)
            %     OR
            %   [obj] = emfieldplot(grid,location,emfieldName,emfieldUnits)
            %         
            %   if no value, initialized with NaN's.
            %   example call to initialize without a value:
            %       obj = emfield(grid,'EDGE');
            %   grid determines cartesian or latlon coordinate type
            %
            %   alternatively, vector may be a ModEMM TVector3D_SG
            %   or TVector3D_Global
            % 
            %   in the future might change this to get the grid from the
            %   TVector3D_SG if that is supplied
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
                obj.location = 'EDGE';
            elseif ischar(varargin{2})
                obj.location = varargin{2};
            elseif isa(varargin{2},'TVector3D_SG') || isa(varargin{2},'TVector3D_Global')
                % get location and value based on the supplied ModEMM.Fwd3D.TVector3D_SG
                vector = varargin{2};
                valueSupplied = 1;
                obj.location = upper(vector.type);
                obj.x = vector.x;
                obj.y = vector.y;
                obj.z = vector.z;
                if vector.grid.Nx == obj.grid.nx && vector.grid.Ny == obj.grid.ny ...
                        && vector.grid.NzEarth == obj.grid.nzEarth && vector.grid.Nza == obj.grid.nzAir
                    % the grids are compatible
                else
                    error('Grid and TVector3D_SG size mismatch on input');
                end
            end
            
            if ~valueSupplied
                switch lower(obj.location)
                    case 'edge'
                        obj.x(1:nx,1:ny+1,1:nz+1) = NaN;
                        obj.y(1:nx+1,1:ny,1:nz+1) = NaN;
                        obj.z(1:nx+1,1:ny+1,1:nz) = NaN;
                    case 'face'
                        obj.x(1:nx+1,1:ny,1:nz) = NaN;
                        obj.y(1:nx,1:ny+1,1:nz) = NaN;
                        obj.z(1:nx,1:ny,1:nz+1) = NaN;
                    otherwise
                        error('only edge or face allowed for EM field location')
                end
            end
            
            if nargin > 2
                obj.emfieldName = varargin{3};
                if strfind(obj.emfieldName,'E')
                    obj.emfieldType = 'electric field';
                elseif strfind(obj.emfieldName,'H') || strfind(obj.emfieldName,'B')
                    obj.emfieldType = 'magnetic field';
                end
            end
            
            if nargin > 3
                obj.emfieldUnits = varargin{4};
            end
           
            disp(['Initialized an EM field of type: ' obj.emfieldType ...
                ', units ' obj.emfieldUnits]);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = xy2ll(obj,lat0,lon0)
            % [obj] = xy2ll(obj,lat0,lon0)
            %
            % convert regular cartesian EM field to a regular lat/lon grid.
            %
            % CURRENTLY ONLY WORKS FOR EDGES
            
            if isa(obj.grid,'xygrid')
                newgrid = llgrid(obj.grid,lat0,lon0);
            elseif isa(obj.grid,'llgrid')
                disp('Nothing to do in xy2ll: EM field already a latlon type');
            else
                error('EM field grid has to be of type ModelPlot xygrid or llgrid');
            end
                        
            Nx = obj.grid.nx;
            Ny = obj.grid.ny;
            Nz = obj.grid.nzEarth + obj.grid.nzAir;
            
            nEx=Nx*(Ny+1)*(Nz+1);
            nEy=(Nx+1)*Ny*(Nz+1);
            nEz=(Nx+1)*(Ny+1)*Nz;
            
            %%%for Ex
            [LON,LAT] = meshgrid(newgrid.lon,newgrid.lat(1:end-1,1)+newgrid.dlat/2);
            ll = [reshape(LAT,[1,Nx*(Ny+1)]);reshape(LON,[1,Nx*(Ny+1)])];
            xy = llgrid.ll2xy(ll,lat0,lon0);
            zz = obj.grid.z';
            
            XI = reshape(xy(1,:)'*ones(1,Nz+1),nEx,1);
            YI = reshape(xy(2,:)'*ones(1,Nz+1),nEx,1);
            ZI = reshape(ones(Nx*(Ny+1),1)*zz,nEx,1);
            
            [Y,X,Z] = meshgrid(obj.grid.y,obj.grid.xctr,obj.grid.z);
            
            Ex = interp3(Y,X,Z,obj.x,YI,XI,ZI,'spline');
            Ex = reshape(Ex,[Nx,Ny+1,Nz+1]);
            obj.x=Ex;
            
            %%%for Ey
            [LON,LAT] = meshgrid(newgrid.lon(1:end-1,1)+newgrid.dlon/2,newgrid.lat);
            
            ll = [reshape(LAT,[1,Ny*(Nx+1)]);reshape(LON,[1,Ny*(Nx+1)])];
            xy = llgrid.ll2xy(ll,lat0,lon0);
            zz = obj.grid.z';            
            
            XI = reshape(xy(1,:)'*ones(1,Nz+1),nEy,1);
            YI = reshape(xy(2,:)'*ones(1,Nz+1),nEy,1);
            ZI = reshape(ones((Nx+1)*Ny,1)*zz,nEy,1);
            
            [Y,X,Z] = meshgrid(obj.grid.yctr,obj.grid.x,obj.grid.z);
            
            Ey = interp3(Y,X,Z,obj.y,YI,XI,ZI,'spline');
            Ey = reshape(Ey,[Nx+1,Ny,Nz+1]);
            obj.y = Ey;
            
            %%%%for Ez
            [LON,LAT] = meshgrid(newgrid.lon,newgrid.lat);
            ll = [reshape(LAT,[1,(Nx+1)*(Ny+1)]);reshape(LON,[1,(Nx+1)*(Ny+1)])];
            xy = llgrid.ll2xy(ll,lat0,lon0);
            zz = obj.grid.zctr;
            
            XI = reshape(xy(1,:)'*ones(1,Nz),nEz,1);
            YI = reshape(xy(2,:)'*ones(1,Nz),nEz,1);
            ZI = reshape(ones((Nx+1)*(Ny+1),1)*(zz'/2),nEz,1);
            
            [Y,X,Z] = meshgrid(obj.grid.y,obj.grid.x,obj.grid.zctr);
            
            Ez = interp3(Y,X,Z,obj.z,YI,XI,ZI,'spline');
            Ez = reshape(Ez,[Nx+1,Ny+1,Nz]);
            obj.z=Ez;
            
            obj.grid = newgrid;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = ll2xy(obj,lat0,lon0)
            % [obj] = ll2xy(obj,lat0,lon0)
            %
            % convert regular lat/lon EM field to a regular grid in km.
            %
            % CURRENTLY ONLY WORKS FOR EDGES

            if isa(obj.grid,'llgrid')
                newgrid = xygrid(obj.grid,lat0,lon0);
            elseif isa(obj.grid,'xygrid')
                disp('Nothing to do in ll2xy: EM field already an xy type');
            else
                error('EM field grid has to be of type ModelPlot xygrid or llgrid');
            end

            Nx = obj.grid.nlat;
            Ny = obj.grid.nlon;
            Nz = obj.grid.nzEarth + obj.grid.nzAir;

            nEx=Nx*(Ny+1)*(Nz+1);
            nEy=(Nx+1)*Ny*(Nz+1);
            nEz=(Nx+1)*(Ny+1)*Nz;

            %%%for Ex
            [Y,X]= meshgrid(newgrid.y,newgrid.xctr);
            xy=[reshape(X,[1,Nx*(Ny+1)]);reshape(Y,[1,Nx*(Ny+1)])];
            ll = xygrid.xy2ll(xy,lat0,lon0);
            zz = obj.grid.z';
            
            LATI = reshape(ll(1,:)'*ones(1,(Nz+1)),nEx,1);
            LONI = reshape(ll(2,:)'*ones(1,(Nz+1)),nEx,1);
            ZI = reshape(ones(Nx*(Ny+1),1)*zz,nEx,1);
            
            [LON,LAT,Z] = meshgrid(obj.grid.lon,obj.grid.lat(1:end-1,:)+obj.grid.dlat/2,zz);

            Ex = interp3(LON,LAT,Z,obj.x,LONI,LATI,ZI,'spline');
            Ex = reshape(Ex,[Nx,Ny+1,Nz+1]);
            obj.x=Ex;
            
            %%% for Ey
            [Y,X]= meshgrid(newgrid.yctr,newgrid.x);
            xy=[reshape(X,[1,(Nx+1)*Ny]);reshape(Y,[1,(Nx+1)*Ny])];
            ll = xygrid.xy2ll(xy,lat0,lon0);
            zz = obj.grid.z';
            
            LATI = reshape(ll(1,:)'*ones(1,(Nz+1)),nEy,1);
            LONI = reshape(ll(2,:)'*ones(1,(Nz+1)),nEy,1);
            ZI = reshape(ones((Nx+1)*Ny,1)*zz,nEy,1);
            
            [LON,LAT,Z] = meshgrid(obj.grid.lon(1:end-1,:)+obj.grid.dlon/2,obj.grid.lat,zz);

            Ey = interp3(LON,LAT,Z,obj.y,LONI,LATI,ZI,'spline');
            Ey = reshape(Ey,[Nx+1,Ny,Nz+1]);
            obj.y=Ey;
            
            %%% for Ez
            [Y,X]= meshgrid(newgrid.y,newgrid.x);
            xy=[reshape(X,[1,(Nx+1)*(Ny+1)]);reshape(Y,[1,(Nx+1)*(Ny+1)])];
            ll = xygrid.xy2ll(xy,lat0,lon0);
            zz = obj.grid.zctr;

            LATI = reshape(ll(1,:)'*ones(1,Nz),nEz,1);
            LONI = reshape(ll(2,:)'*ones(1,Nz),nEz,1);
            ZI = reshape(ones((Nx+1)*(Ny+1),1)*zz',nEz,1);

            [LON,LAT,Z] = meshgrid(obj.grid.lon,obj.grid.lat,zz);

            Ez = interp3(LON,LAT,Z,obj.z,LONI,LATI,ZI,'spline');
            Ez = reshape(Ez,[Nx+1,Ny+1,Nz]);
            obj.z=Ez;
            
            obj.grid = newgrid;

        end
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function obj = fillNaNs(obj,method,bg)
%             % obj = fillNaNs(obj,method,bg)
%             %
%             % replace all NaNs with the depth-averaged value, nearest neighbor
%             % value, or the specified depth-specific bg value; comes handy
%             % after an llmodel to xymodel conversion, or vice versa
%             %
%             % if optional background bg is specified, it is linear (not log) 
%             % and is either scalar, or has the same number of values as 
%             % as the number of vertical layers in the model. 
%             % use bg = 0 for conductivity perturbations  
%             %
%             % options: 'mean', 'background' or 'nearest'; default 'nearest'
%             
%             if nargin < 2
%                 method = 'nearest';
%             end
%             
%             nz = size(obj.v,3);
%             if strcmp(method,'background')
%                 if nargin < 3
%                     warning('Background conductivity or resistivity has not been specified. Using depth average.');
%                     method = 'mean';
%                 else
%                     if isscalar(bg)
%                         bg(1:nz) = bg;
%                     end
%                     if strcmp(obj.paramType,'LOGE') && (min(bg) > 0)
%                         bg = log(bg);
%                     end
%                     if strcmp(obj.paramType,'LOG10') && (min(bg) > 0)
%                         bg = log10(bg);
%                     end
%                 end
%             end
%             
%             nantotal = length(find(isnan(obj.v)));
%             if nantotal > 0
%                 disp(['Using method ''' method ''' to fill in ' num2str(nantotal) ' NaNs.']);
%             else
%                 disp('No NaNs to fill in.');
%             end
%             for k = 1:size(obj.v,3)
%                 vk = squeeze(obj.v(:,:,k));
%                 if strcmp(method,'mean')  
%                     vk(isnan(vk)) = nanmean(nanmean(vk));
%                     obj.v(:,:,k) = vk;
%                 elseif strcmp(method,'background')
%                     vk(isnan(vk)) = bg(k);
%                     obj.v(:,:,k) = vk;
%                 elseif strcmp(method,'nearest')
%                     nlon = size(obj.v,2);
%                     nlat = size(obj.v,1);
%                     nantotal = length(find(isnan(vk)));
%                     while nantotal > 0
%                         for i = 1:size(obj.v,2)
%                             for j = 1:size(obj.v,1)
%                                 if isnan(vk(j,i))
%                                     imin = max(1,i-3);
%                                     imax = min(nlon,i+3);
%                                     jmin = max(1,j-3);
%                                     jmax = min(nlat,j+3);
%                                     nearest = vk(jmin:jmax,imin:imax);
%                                     obj.v(j,i,k) = nanmean(nanmean(nearest));
%                                 end
%                             end
%                         end
%                         vk = squeeze(obj.v(:,:,k));
%                         nantotal = length(find(isnan(vk)));
%                     end
%                 end   
%                 
%             end
%             
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = lon360(obj)
            % convert the model parameter value to longitude [0,360)
            % convention
            
            if isnumeric(obj)
                obj = llgrid.lon360(obj);
            else
                midlon = llgrid.lon360(obj.grid.lon(1:end-1)+obj.grid.dlon/2);
                [midlon,ii] = sort(midlon);
                obj.v = obj.v(:,ii,:);
                obj.grid.lon = sort(llgrid.lon360(obj.grid.lon));
                obj.grid.dlon = diff(obj.grid.lon);
                % unreasonable inputs with even a single outlier produce a
                % completely broken set of values - do not use this
                %[obj.grid.dlon,obj.grid.lon] = obj.delta(midlon);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = lon180(obj)
            % convert the model parameter value to longitude (-180,180]
            % convention
           
            if isnumeric(obj)
                obj = llgrid.lon180(obj);
            else
                midlon = llgrid.lon180(obj.grid.lon(1:end-1)+obj.grid.dlon/2);
                [midlon,ii] = sort(midlon);
                obj.v = obj.v(:,ii,:);
                obj.grid.lon = sort(llgrid.lon180(obj.grid.lon));
                obj.grid.dlon = diff(obj.grid.lon);
                % unreasonable inputs with even a single outlier produce a
                % completely broken set of values - do not use this
                %[obj.grid.dlon,obj.grid.lon] = obj.delta(midlon);
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setOrigin(obj,lat0,lon0)
            % define the geographic origin and corresponding irregular 
            % spherical grid coordinates LAT & LON
            
            obj.lat0 = lat0;
            obj.lon0 = lon0;
            
        end     
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = set.limits(obj,lims)
            %   set lat, lon, depth limits
            if ~isstruct(lims)
                error('Limits must be a structure that contains min and max values in lon, lat and depth')
            end
            obj.limits.latmin = lims.latmin;
            obj.limits.latmax = lims.latmax;
            obj.limits.lonmin = llgrid.lon360(lims.lonmin);
            obj.limits.lonmax = llgrid.lon360(lims.lonmax);
            if isfield(obj.limits,'depthmin') && isfield(obj.limits,'depthmax')
                obj.limits.depthmin = max(lims.depthmin,obj.limits.depthmin);
                obj.limits.depthmax = min(lims.depthmax,obj.limits.depthmax);
            else
                obj.limits.depthmin = lims.depthmin;
                obj.limits.depthmax = lims.depthmax;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function lims = get.limits(obj)
            %   find lat, lon, depth limits
            if isstruct(obj.limits)
                lims = obj.limits;
            else
                if isa(obj,'xymodel')
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
            if obj.isglobal
                lims.lonmin = lims.lonmin;
                lims.lonmax = lims.lonmax;
            else
                lims.lonmin = llgrid.lon360(lims.lonmin);
                lims.lonmax = llgrid.lon360(lims.lonmax);
            end   
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [h,X] = pcolor(obj,type,value,component,realORimag,lims,P)
            % Sample usage: h = obj.pcolor('depth',16)
            % Optional structure P contains plotting parameters.
            % By default does 1D interpolation which can be switched off
            % by setting P.nointerpbyvalue = 1.
            % Full example:
            % obj = llmodel.read('SEMUM_percent.nc','netcdf');
            % [lims,posn] = llmodel.getLimits;
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
            
            if strcmp(obj.location,'EDGE')
                if strcmpi(component,'x')
                    lon = obj.grid.lon;
                    lat = obj.grid.lat(1:end-1) + obj.grid.dlat/2;
                    depth = obj.grid.depth;
                elseif strcmpi(component,'y')
                    lon = obj.grid.lon(1:end-1) + obj.grid.dlon/2;
                    lat = obj.grid.lat;
                    depth = obj.grid.depth;
                elseif strcmpi(component,'z')
                    lon = obj.grid.lon;
                    lat = obj.grid.lat;
                    depth = obj.grid.depth(1:end-1) + obj.grid.dz/2;
                end
            elseif strcmp(obj.location,'FACE')
                if strcmpi(component,'x')
                    lon = obj.grid.lon(1:end-1) + obj.grid.dlon/2;
                    lat = obj.grid.lat;
                    depth = obj.grid.depth(1:end-1) + obj.grid.dz/2;
                elseif strcmpi(component,'y')
                    lon = obj.grid.lon;
                    lat = obj.grid.lat(1:end-1) + obj.grid.dlat/2;
                    depth = obj.grid.depth(1:end-1) + obj.grid.dz/2;
                elseif strcmpi(component,'z')
                    lon = obj.grid.lon(1:end-1) + obj.grid.dlon/2;
                    lat = obj.grid.lat(1:end-1) + obj.grid.dlat/2;
                    depth = obj.grid.depth;
                end
            else
                error('Set EM soln location to EDGE or FACE before plotting.')
            end
                
            nlon = length(lon);
            nlat = length(lat);
            
            depth = [- sort(obj.grid.zAir,1,'descend'); obj.grid.depth(2:end)];
            nz = length(depth);
           
            if strcmp(type,'depth') 
                if issorted(depth)
                    k2 = find(depth>value,1);
                    k = max(k2-1,1);
                else
                    k2 = find(depth>value,1,'last');
                    k = min(k2+1,nz);
                end
                alpha = (depth(k2)-value)/(depth(k2)-depth(k));
                X = obj.(component)(:,:,k);
                nointerp = 0;
                if nargin > 6
                    if isfield(P,'nointerpbyvalue')
                        nointerp = P.nointerpbyvalue;
                    end
                end
                if ~nointerp
                    X2 = obj.(component)(:,:,k2);
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
                X = obj.(component)(:,:,k);
                if isfield(obj,'depth_pr') % i.e., global model object
                    if length(obj.depth_pr) == size(obj.v,3)
                        disp(['Plotting the value for layer down to depth ' num2str(obj.depth_pr(k)) ' km']);
                    else
                        disp(['Plotting the value for layer index ' num2str(k) '; depth unknown']);
                    end
                else
                    disp(['Plotting the EM field solution value for layer index ' num2str(k)]);
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
            
            if nargin < 5
                realORimag = 'abs';
            end
            if strcmp(realORimag,'abs') 
                X = abs(X); 
            elseif strcmp(realORimag,'phase')
                X = phase(X);
            elseif strcmp(realORimag,'real')
                X = real(X);
            elseif strcmp(realORimag,'imag')
                X = imag(X);
            end

            proj = 'miller';
            delta = 1e-6;
            clims = [nanmin(nanmin(real(X)))-delta nanmax(nanmax(real(X)))+delta];
            cmap = 'jet';
            nointerp = 0;
            showmesh = 0;
            regional = 1;
            clong = 180;
            titleOff = 0;
            if nargin > 6
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
            end
            
            if nargin < 6
                lonmin = obj.grid.limits.lonmin;
                lonmax = obj.grid.limits.lonmax;
                latmin = obj.grid.limits.latmin;
                latmax = obj.grid.limits.latmax;
            else
                lonmin = lims.lonmin;
                lonmax = lims.lonmax;
                latmin = lims.latmin;
                latmax = lims.latmax;
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
            if ~nointerp; shading flat; end
            if showmesh; shading faceted; end
            colorbar('FontSize',28);
            hold on;
            if strcmp(obj.coords,'geographic')
                cbndry  = [0.1 0.1 0.1];
                %m_gshhs_i('color',cbndry);              % Coastline... higher resolution
                %m_gshhs_i('speckle','color',cbndry);    % with speckle added
                m_coast('color',cbndry);              % Coastline...
                m_coast('speckle','color',cbndry);    % with speckle added
                %plotStates
            end
            if ~isempty(strfind(lower(proj),'miller')) ...
                    || ~isempty(strfind(lower(proj),'equidistant'))
                m_grid('box','fancy','tickdir','in',... %'box','on',...
                    'xtick',2,'linestyle','none',...
                    'XaxisLocation','top',...
                    'fontsize',28,'fontweight','demi');                
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
                    title([obj.emfieldName upper(component) ' at ' mydepth ' km'],'fontsize',38,'fontweight','demi');
                end
            end

        end
        
    end
    
    methods(Static)
        % empty
    end
end

    
                                
   
                
                