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
        paramType
        z
        X
        Y
        Z
        CONDll    %   this is used to interpolate the model parameter onto
                  %   a regular (tensor product) latitude/longitude grid;
                  %   this  will generally be a subset of the xyz grid,
                  %   since the minimum width in longitude has to b used to
                  %   define the model parameter on a regular grid
        llgrid    %   grid (lat, lon, z) used for interpolated model
    end
    
    methods
        function [obj] = Cond_xy_ll(Cond,lat0,lon0)
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
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function result = limits(obj)
            %   find lat, lon, depth limits
            result = struct('latmin',min(min(obj.LAT(:,:,1))),...
                            'latmax',max(max(obj.LAT(:,:,1))),...
                            'lonmin',min(min(obj.LON(:,:,1))),...
                            'lonmax',max(max(obj.LON(:,:,1))),...
                            'depthmin',min(obj.Z(1,1,:)),...
                            'depthmax',max(max(obj.Z(1,1,:))));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function result = latSlice(obj,latitudes,varargin)
            %   Usage: Corners = latSlice(obj,'longitudes',[lonMin,lonMax] ...
            %                                 ,'depths','[depthMin,depthMax])
            %
            %   make corners for slice at given latitude; if optional
            %   "depths" or "longitudes" arguments are provided use these to
            %   define corners of slice; otherwise use limits obtained from
            %   object.   If latitudes is a vector, a cell array of section
            %   corners is returned
            lims = obj.limits;
            depths = [ lims.depthmin, lims.depthmax];
            lons = [lims.lonmin lims.lonmax];
            nsections = length(latitudes);
            if nargin > 2
                n = length(varargin);
                if mod(n,2)
                    error('Optional arguments must occur in pairs')
                end
                for k = 1:2:n
                    option = lower(varargin{k});
                    switch option
                        case 'depths'
                            depths = varargin{k+1};
                        case 'longitudes'
                            lons = varargin{k+1};
                            lons  = mod(lons,360)
                        otherwise
                            error('Optional argument not defined')
                    end
                end
            end
            Corners = zeros(3,2,2);
            %   depths
            Corners(3,1,:) = depths(1);
            Corners(3,2,:) = depths(2);
            %  longitude
            Corners(2,:,1) = lons(1);
            Corners(2,:,2) = lons(2);
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
        function Corners = lonSlice(obj,longitudes,varargin)
            %   Usage: Corners = latSlice(obj,'latitudes',[latMin,latMax] ...
            %                                 ,'depths','[depthMin,depthMax])
            %
            %   make corners for slice at given longitude; if optional
            %   "depths" or "latitudes" arguments are provided use these to
            %   define corners of slice; otherwise use limits obtained from
            %   object.  If longitudes is a vector, a cell array of section
            %   corners is returned
            lims = obj.limits;
            depths = [ lims.depthmin, lims.depthmax];
            lats = [lims.latmin lims.latmax];
            longitudes = mod(longitudes,360);
            nsections = length(longitudes); 
            if nargin > 2
                n = length(varargin);
                if mod(n,2)
                    error('Optional arguments must occur in pairs')
                end
                for k = 1:2:n
                    option = lower(varargin{k});
                    switch option
                        case 'depths'
                            depths = varargin{k+1};
                        case 'latitudes'
                            lats = varargin{k+1};
                        otherwise
                            error('Optional argument not defined')
                    end
                end
            end
            Corners = zeros(3,2,2);
            %   depths
            Corners(3,1,:) = depths(1);
            Corners(3,2,:) = depths(2);
            %  latitude
            Corners(1,:,1) = lats(1);
            Corners(1,:,2) = lats(2);
            if nsections > 1
                result = cell(nsections,1);
                for k = 1:nsections
                    Corners(1,:,:) = longitudes(k);
                    result{k} = struct('Corners',Corners);
                end
            else
                Corners(1,:,:) = longitudes;
                result = Corners;
            end
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
            ind = find(obj.LON(:,1,1) == max(obj.LON(:,1,1)),1);
            lon = obj.LON(ind,:,1);obj
            %   given vertical fields
            z = squeeze(obj.Z(1,1,:));
            if nargin > 1
                n = length(varargin);
                if mod(n,2)
                    error('Optional arguments must occur in pairs')
                end
                for k = 1:2:n
                    option = lower(varargin{k});
                    switch option
                        case 'depths'
                            z = varargin{k+1};
                        case 'latitudes'
                            lat = varargin{k+1};
                        case 'longitudes'
                            lon = varargin{k+1};
                        otherwise
                            error('Optional argument not defined')
                    end
                end
            end
            obj.llgrid = struct('lat',lat,'lon',lon,'z',z);
            [LON,LAT] = meshgrid(lon,lat);
            nlon = length(lon);
            nlat = length(lat);
            nz = length(z);
            nTot = nlon*nlat*nz;
            ll = [reshape(LAT,[1,nlon*nlat]);reshape(LON,[1,nlon*nlat])];
            xy = ll2xy(obj,ll);
            XI = reshape(xy(1,:)'*ones(1,length(z)),nTot,1);
            YI = reshape(xy(2,:)'*ones(1,length(z)),nTot,1);
            temp = reshape(obj.llgrid.z,1,nz);
            ZI = reshape(ones(nlat*nlon,1)*temp,nTot,1);
            obj.CONDll = interp3(obj.Y,obj.X,obj.Z,obj.CONDxyz,...
                YI,XI,ZI);
            obj.CONDll = reshape(obj.CONDll,[nlat,nlon,nz]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fencePlot(obj,sections,varargin)
            
            figure('Position',[100,100,600,600])
            nslice = length(sections);
            for k = 1:nslice
                [cond,xy,ll] = slice(obj,sections{k}.Corners,sections{k}.NM);
                surf(squeeze(ll(2,:,:)),squeeze(ll(1,:,:)),squeeze(ll(3,:,:)),cond);
                shading flat
                hold on
            end
            set(gca,'Fontweight','demi','FontSize',14,'zdir','reverse')
        end
    end
end

    
                                
   
                
                