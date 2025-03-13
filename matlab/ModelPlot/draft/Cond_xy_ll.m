classdef Cond_xy_ll < latgrid
    %   addd depth and conductivity to the latgrid class, to store
    %   conductivity on a standard cartesian grid with lat/lon info, and
    %   tools for extracting slices in lat/lon/depth coordinates
    
    properties
        CONDxyz    %   this is used to store the model parameter on a 
                   %   regular (tensor product) cartesian grid, with
                   %   coordinates defined by x, y, z --- i.e., the normal
                   %   structure used for modeling and inversion in
                   %   Cartesian coordinates
        modelType  % conductivity or resistivity
        modelUnits % S/m or Ohm*m
        paramType  % LOGE, LOG10 or LINEAR
        AirCond = 1e-10 
        z
        X
        Y
        Z
        CONDll    %   this is used to interpolate the model parameter onto
                  %   a regular (tensor product) latitude/longitude grid;
                  %   this  will generally be a subset of the xyz grid,
                  %   since the minimum width in longitude has to be used 
                  %   to define the model parameter on a regular grid
                  %   NOT used inside class for slice interpolation, just
                  %   for output and plotting
        llgrid    %   grid (lat, lon, depth) used for interpolated model
        limits    %   depthmin, depthmax, latmin, latmax, lonmin, lonmax
                  %   using -180<lon<=180 for output & [0,360] for interp
    end
    
    properties (Dependent = true)
    end
    
    methods
        function [obj] = Cond_xy_ll(Cond,lat0,lon0,varargin)
            %   class constructor    
            obj = obj@latgrid(Cond.grid,lat0,lon0);
 
            if strcmp(Cond.paramType,'LOGE')
                obj.CONDxyz = Cond.v/log(10);
                obj.paramType = 'LOG10';
            else
                obj.paramType = Cond.paramType;
                obj.CONDxyz = Cond.v; 
            end
            obj.z = [0 ; cumsum(Cond.grid.dz)]+Cond.grid.origin(3);
            obj.z = (obj.z(1:end-1)+obj.z(2:end))/2;
            Nxy = Cond.grid.Nx*Cond.grid.Ny;
            obj.LAT = reshape(obj.LAT,Nxy,1)*ones(1,Cond.grid.NzEarth);
            obj.LAT = reshape(obj.LAT,Cond.grid.Nx,Cond.grid.Ny,Cond.grid.NzEarth);
            obj.LON = reshape(obj.LON,Nxy,1)*ones(1,Cond.grid.NzEarth);
            obj.LON = reshape(obj.LON,Cond.grid.Nx,Cond.grid.Ny,Cond.grid.NzEarth);
            [obj.Y,obj.X,obj.Z] = meshgrid(obj.y,obj.x,obj.z);
            obj = CONDxyz2ll(obj,varargin{:});
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = set.limits(obj,lims)
            %   set lat, lon, depth limits
            if ~isstruct(lims)
                error('Limits must be a structure that contains min and max values in lon, lat and depth')
            end
            obj.limits.latmin = lims.latmin;
            obj.limits.latmax = lims.latmax;
            obj.limits.lonmin = latgrid.lon360(lims.lonmin);
            obj.limits.lonmax = latgrid.lon360(lims.lonmax);          
            obj.limits.depthmin = max(lims.depthmin,obj.limits.depthmin);
            obj.limits.depthmax = min(lims.depthmax,obj.limits.depthmax);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function lims = get.limits(obj)
            %   find lat, lon, depth limits
            if isstruct(obj.limits)
                lims = obj.limits;
            else
                lims = struct(...
                    'latmin',min(min(obj.LAT(:,:,1))),...
                    'latmax',max(max(obj.LAT(:,:,1))),...
                    'lonmin',min(min(obj.LON(:,:,1))),...
                    'lonmax',max(max(obj.LON(:,:,1))),...
                    'depthmin',min(obj.Z(1,1,:)),...
                    'depthmax',max(obj.Z(1,1,:)));
            end
            lims.lonmin = latgrid.lon360(lims.lonmin);
            lims.lonmax = latgrid.lon360(lims.lonmax);
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
                longitudes = [ lims.lonmin, lims.lonmax];
            end
            if nargin < 4
                latitudes = [ lims.latmin, lims.latmax];
            end
            Corners = zeros(3,2,2);
            %  latitude
            Corners(1,1,:) = latitudes(1);
            Corners(1,2,:) = latitudes(2);
            NM(1) = sum(obj.llgrid.lat>=latitudes(1) & obj.llgrid.lat<=latitudes(2));
            %  longitude
            Corners(2,:,1) = longitudes(1);
            Corners(2,:,2) = longitudes(2);
            NM(2) = sum(obj.llgrid.lon>=longitudes(1) & obj.llgrid.lon<=longitudes(2));
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
            if nsections > 1
                result = cell(nsections,1);
                for k = 1:nsections
                    Corners(1,:,:) = latitudes(k);
                    result{k} = struct('Corners',Corners);
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
            xy = ll;
            xy(1:2,:) = ll2xy(obj,ll(1:2,:));
            cond = interp3(obj.Y,obj.X,obj.Z,obj.CONDxyz,...
                xy(2,:),xy(1,:),xy(3,:));
            cond = reshape(cond,[NM(1),NM(2)]);
            xy = reshape(xy,[3,NM(1),NM(2)]);
            ll = reshape(ll,[3,NM(1),NM(2)]);
        end  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = CONDxyz2ll(obj,varargin)
            % convert model parameter to regular lat/lon grid
            %   defaults: largest set of lats and lons possible
            lat = obj.LAT(:,1,1);
            ind = find(obj.LON(:,1,1) == max(obj.LON(:,1,1)), 1);
            lon = obj.LON(ind,:,1);
            %   given vertical fields
            depth = squeeze(obj.Z(1,1,:));
            if nargin > 1
                n = length(varargin);
                if mod(n,2)
                    error('Optional arguments must occur in pairs')
                end
                for k = 1:2:n
                    option = lower(varargin{k});
                    switch option
                        case 'depths'
                            depth = varargin{k+1};
                        case 'latitudes'
                            lat = varargin{k+1};
                        case 'longitudes'
                            lon = varargin{k+1};
                        otherwise
                            error('Optional argument not defined')
                    end
                end
            end
            obj.llgrid = struct('lat',lat,'lon',lon,'depth',depth);
            obj.llgrid.lon = latgrid.lon180(lon);
            [LON,LAT] = meshgrid(lon,lat);
            nlon = length(lon);
            nlat = length(lat);
            ndepth = length(depth);
            nTot = nlon*nlat*ndepth;
            ll = [reshape(LAT,[1,nlon*nlat]);reshape(LON,[1,nlon*nlat])];
            xy = ll2xy(obj,ll);
            XI = reshape(xy(1,:)'*ones(1,ndepth),nTot,1);
            YI = reshape(xy(2,:)'*ones(1,ndepth),nTot,1);
            temp = reshape(obj.llgrid.depth,1,ndepth);
            ZI = reshape(ones(nlat*nlon,1)*temp,nTot,1);
            obj.CONDll = interp3(obj.Y,obj.X,obj.Z,obj.CONDxyz,...
                YI,XI,ZI);
            obj.CONDll = reshape(obj.CONDll,[nlat,nlon,ndepth]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fencePlot(obj,sections,varargin)
            
            figure('Position',[200,200,1400,900],...
                    'PaperPosition',[1,1,17.8,10]); %widescreen
            nslice = length(sections);
            for k = 1:nslice
                [cond,xy,ll] = slice(obj,sections{k}.Corners,sections{k}.NM); %#ok<ASGLU>
                surf(squeeze(ll(2,:,:)-360),squeeze(ll(1,:,:)),squeeze(ll(3,:,:)),cond);
                shading flat
                hold on
            end
            set(gca,'Fontweight','demi','FontSize',14,'zdir','reverse')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [CONDll,llgrid] = select(obj,lims)
            
            if nargin < 2
                lims = obj.limits;
            end
            k1 = find(obj.llgrid.depth >= lims.depthmin, 1, 'first');
            k2 = find(obj.llgrid.depth <= lims.depthmax, 1, 'last');
            i1 = find(obj.llgrid.lat >= lims.latmin, 1, 'first');
            i2 = find(obj.llgrid.lat <= lims.latmax, 1, 'last');
            j1 = find(obj.llgrid.lon >= lims.lonmin, 1, 'first');
            j2 = find(obj.llgrid.lon <= lims.lonmax, 1, 'last');
            CONDll = obj.CONDll(i1:i2,j1:j2,k1:k2);
            
            llgrid = obj.llgrid;
            llgrid.depth = obj.llgrid.depth(k1:k2);
            llgrid.lat = obj.llgrid.lat(i1:i2);
            llgrid.lon = obj.llgrid.lon(j1:j2);
            llgrid.nzEarth = length(k1:k2);
            llgrid.nlat = length(i1:i2);
            llgrid.nlon = length(j1:j2);
            llgrid.limits = lims;

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function write(obj,fname,format,lims)
            % format options: Slices, NetCDF
            
            if nargin < 3
                format = 'Slices';
            end
            if nargin < 4
                lims = obj.limits;
            end
            % [CONDll,llgrid] = select(obj,lims);
            switch format
                case 'Slices'
                    if ~exist(fname,'dir'); mkdir(fname); end
                    k = find(obj.llgrid.depth >= lims.depthmin, 1);
                    while obj.llgrid.depth(k) <= lims.depthmax
                        str = sprintf('%06.2f',obj.llgrid.depth(k));
                        f = fopen([fname '/' fname '_' str 'km.cond'],'w');
                        i = find(obj.llgrid.lat >= lims.latmin, 1);
                        while obj.llgrid.lat(i) <= lims.latmax
                            j = find(obj.llgrid.lon >= lims.lonmin, 1);
                            while obj.llgrid.lon(j) <= lims.lonmax
                                fprintf(f,'%f\t%f\t%g\n',obj.llgrid.lat(i),obj.llgrid.lon(j),obj.CONDll(i,j,k));
                                j = j+1;
                            end
                            i = i+1;
                        end
                        fclose(f);
                        k = k+1;
                    end
                case 'NetCDF'
                    disp('Not implemented yet');
                otherwise
                    error(['Can''t write to lat/lon model format ' format ': method unknown']);
            end

        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     methods(Static)
%         function write_NetCDF(CONDll,llgrid,fname)
%             
%         end
%     end
end

    
                                
   
                
                