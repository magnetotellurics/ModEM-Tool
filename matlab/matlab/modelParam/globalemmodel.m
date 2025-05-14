classdef globalemmodel < llmodel
    %   Global EM/electrical conductivity model may be defined in layers or
    %   on grid, both options are supported. It can be in geographic or
    %   geomagnetic coordinates. It always spans the whole Earth and has
    %   padding at zero longitude. Otherwise, inherited from llmodel.
    
    properties
        %grid       % llgrid
        nlayers
        headers
        prm
        prm_pr
        depth_pr    % depth to the bottom of each layer
        %v          % 3D array in layers or on grid
        %v_ref      % 1D array
        v_pr        % 3D array in layers or on grid
    end
    
    properties (Constant)
        % Earth radius in km used for global EM modeling
        EarthRadGlobal = 6371.0d0;
    end    
    
    methods
        function [obj] = globalemmodel(varargin)
            %   class constructor
            obj = obj@llmodel(varargin{:});
            obj.modelName  = '\sigma';
            obj.modelType  = 'electrical conductivity';
            obj.modelUnits = 'S/m';
            obj.displayName = 'log(10) \sigma';
            obj.paramType = 'LOG10';
            obj.coords = 'geomagnetic';
            obj.nlayers = 0;
            obj.isglobal = 1;

            if nargin >= 0
                obj.grid = llgrid;
                return
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = read_grid(obj,fname)
            if nargin < 2
                [fname,fdir] = uigetfile('*.grd', 'Global Grid');
            else
                fdir = './';
            end
            ggrid = read_grid([fdir '/' fname]);
            obj.grid = llgrid;
            obj.grid.lat = 90-ggrid.y;
            obj.grid.lon = ggrid.x;
            obj.grid.depth = obj.EarthRadGlobal - ggrid.z(ggrid.nzAir+1:end);
            obj.grid.dlat = abs(diff(obj.grid.lat));
            obj.grid.dlon = abs(diff(obj.grid.lon));
            obj.grid.dz = abs(diff(obj.grid.depth));
            obj.grid.units = 'km';
            obj.grid.zAir = ggrid.z(1:ggrid.nzAir);
            obj.grid.nzAir = ggrid.nzAir;
            obj.grid.nzCrust = ggrid.nzCrust;
            obj.grid.nzEarth = ggrid.nzEarth;
            obj.grid.nlat = ggrid.ny;
            obj.grid.nlon = ggrid.nx;
            % correction to obtain mid-points of cells
            obj.grid.lat = obj.grid.lat(1:end-1)-obj.grid.dlat/2;
            obj.grid.lon = obj.grid.lon(1:end-1)+obj.grid.dlon/2;
            obj.grid.depth = obj.grid.depth(1:end-1)+obj.grid.dz/2;
            % correction for zero longitude
            %obj.grid.lon(end+1) = 360 + obj.grid.lon(1);
            %obj.grid.dlon(end+1) = obj.grid.dlon(1);
            %obj.grid.nlon = obj.grid.nlon+1;
            % update v_pr if needed
            if isempty(obj.v_pr) && ~isempty(obj.v_ref) && obj.nlayers > 0
                obj.v_pr = zeros(obj.grid.nlat,obj.grid.nlon,obj.nlayers);
                for k=1:obj.nlayers
                    obj.v_pr(:,:,k) = obj.v_ref(k);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = write_grid(obj,fname_grid)
            % routine to write an earth3d *.grd file
            obj = obj.lon360;
            dx = obj.grid.dlon;                        % longitude 0..360
            y = 90 - obj.grid.lat(1:end);              % co-latitude 0..180
            if strcmp(obj.grid.units,'km')
                z = obj.EarthRad - obj.grid.depth(1:end);        % z (km)
            else
                z = obj.EarthRad - obj.grid.depth(1:end)/1000.0; % z (km)
            end
            nx = obj.grid.nlon;
            ny = obj.grid.nlat;
            nzAir = 17;
            zAir = [19113.1 19113.0 14563. 10467. 8419. 7395. 6883. ...
                    6627. 6499. 6435. 6403. 6387. 6379. 6375. 6373. 6372. 6371.5];
            nzCrust = obj.grid.nzCrust;
            nzEarth = obj.grid.nzEarth;
            % write to a simple earth3d *.rho ascii file
            fid=fopen(fname_grid,'w');
            fprintf(fid,'%d %d %d %d %d\n\n',nx,ny,nzAir,nzCrust,nzEarth);
            for i = 1:nx
                fprintf(fid,'%6.3f ',dx(i));
                if mod(i,10)==0; fprintf(fid,'\n'); end
            end
            fprintf(fid,'\n\n');
            for i = ny+1:-1:1 
                fprintf(fid,'%6.3f ',y(i)); 
                if mod(ny+2-i,10)==0; fprintf(fid,'\n'); end
            end
            fprintf(fid,'\n\n');
            for i = 1:nzAir
                fprintf(fid,'%6.3f ',zAir(i));
                if mod(i,10)==0; fprintf(fid,'\n'); end
            end
            fprintf(fid,'\n');
            for i = 1:nzCrust+nzEarth+1
                fprintf(fid,'%6.3f ',z(i)); 
                if mod(i,10)==0; fprintf(fid,'\n'); end
            end
            fclose(fid);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use to read the model parameterization; use pr for prior model 
        function [obj] = read_prm(obj,fname)
            if nargin < 2
                [fname,fdir] = uigetfile('*.prm', 'Global Electrical Resistivity Model Coefficients');
            else
                fdir = './';
            end
            obj.prm = read_prm([fdir '/' fname]);
            obj.prm = cat(1,obj.prm{:});
            obj.depth_pr = cat(1,obj.prm.depth);
            d_up=[0; obj.depth_pr(1:end-1)]; d_lo=obj.depth_pr;
            obj.headers=strcat(num2str(d_up),'-',strjust(num2str(d_lo),'left'));
            obj.nlayers=length(obj.prm);
            obj.fileName = [fdir fname];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use to read the gridded model; use pr for prior model 
        function [obj] = read_rho(obj,fname_rho3d)
            if nargin < 2
                [fname_rho3d,fdir] = uigetfile('*.rho', 'Global Electrical Resistivity Model');
            else
                fdir = './';
            end
            % note: optional output arguments are superceded with *grd file
            obj.v=read_rho([fdir fname_rho3d]);
            obj.v= - log10(obj.v);
            % correction for zero longitude
            %obj.v(:,end+1,:) = obj.v(:,1,:);
            obj.fileName = [fdir fname_rho3d];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use to write the gridded model: want to write linear resistivity.
        % Simplified version - at present will assume a regular grid
        function [obj] = write_rho(obj,fname_rho3d)
            % convert to lon360
            obj = obj.lon360;
            % convert to linear (otherwise, type = 'LOG10')
            if strfind(obj.paramType,'LOG10') && strfind(obj.modelType,'conductivity')
                value = 10.^(- obj.v);
            elseif strfind(obj.modelType,'conductivity')
                value = 1./obj.v;
            elseif strfind(obj.paramType,'LOG10')
                value = 10.^(  obj.v);
            else
                value = obj.v;
            end
            type = '';
            % to ignore duplicate zero longitudes, if present
            ilon  = intersect(find(obj.grid.lon <  obj.grid.limits.lonmax), ...
                              find(obj.grid.lon >= obj.grid.limits.lonmin));
            rho = value(:,ilon,:);
            % generalized write_rho(fname_rho3d,value,depth);
            x = obj.grid.lon(ilon);                        % longitude 0..360
            y = obj.grid.lat(2:end);                       % latitude
            if strcmp(obj.grid.units,'km')
                z = obj.grid.depth(1:end-1);               % depth (km)
            else
                z = obj.grid.depth(1:end-1)/1000.0;        % depth (km)
            end
            nx = length(ilon);
            ny = obj.grid.nlat;
            nz = obj.grid.nzCrust+obj.grid.nzEarth;
            % write to a simple earth3d *.rho ascii file
            fid=fopen(fname_rho3d,'w');
            fprintf(fid,'# lon(i), lat(j), depth(k), rho(ijk) written by Matlab\n');
            fprintf(fid,'%d %d %d %s\n',nx,ny,nz,type);
            for k = 1:nz
                for i = 1:nx
                    for j = ny:-1:1
                        fprintf(fid,'%f\t%f\t%f\t%f\n',x(i),y(j),z(k),rho(j,i,k));
                    end
                end
            end
            fclose(fid);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = pr(obj,fname,expand)
            if nargin < 2
                [fname,fdir] = uigetfile('*.prm;*.rho', 'Prior Electrical Resistivity Model');
            else
                fdir = './';
            end
            if nargin < 3
                expand = 0;
            end
            if strfind(fname,'.prm')
                obj.prm_pr = read_prm([fdir '/' fname]);
                obj.prm_pr = cat(1,obj.prm_pr{:});
                obj.v_ref = - cat(1,obj.prm_pr.logrho);
                obj.depth_pr = cat(1,obj.prm_pr.depth);
                d_up=[0; obj.depth_pr(1:end-1)]; d_lo=obj.depth_pr;
                obj.headers=strcat(num2str(d_up),'-',strjust(num2str(d_lo),'left'));
                obj.nlayers=length(obj.prm_pr);
                if obj.grid.nlat>0 && obj.grid.nlon>0
                    obj.v_pr = zeros(obj.grid.nlat,obj.grid.nlon,obj.nlayers);
                    for k=1:obj.nlayers
                        obj.v_pr(:,:,k) = obj.v_ref(k);
                    end
                end
                if expand
                    obj.v_pr = obj.expand(obj.v_pr);
                end
            elseif strfind(fname,'.rho')
                error(['Unable to read the file ' fname ' into prior model']);
            elseif strfind(fname,'.nc')
                error(['Unable to read the file ' fname ' into prior model']);
            else
                error(['Unable to read the file ' fname ' into prior model']);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use to quickly plot an integer layer; otherwise use obj.pcolor
        function plot(obj,k,cx)
            figure; 
            pcolor(obj.grid.lon,obj.grid.lat,obj.v(:,:,k)); shading interp;
            fsize = 16;
            set(gca,'fontsize',fsize)
            if nargin>2
                caxis(cx);
            end
            cbar=colorbar;
            set(cbar,'fontsize',fsize)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use to squeeze values into layers; squeezes ALL so that a
        % variable crust will also be squeezed; fix if needed.
        function newvalue = squeeze(obj,value)
            nz = obj.grid.nzEarth + obj.grid.nzCrust;
            if obj.nlayers == 0 || isempty(obj.nlayers)
                error('Unable to squeeze into layers: no layers defined');
            elseif nz == 0
                error('Unable to squeeze into layers: no grid defined');
            elseif nz ~= size(value,3)
                error('Value size doesn''t match the grid');
            end
            newvalue = zeros(size(value,1),size(value,2),obj.nlayers);
            i1 = 1;
            for k = 1:obj.nlayers
                ctrdepth = obj.grid.depth(1:end-1) + obj.grid.dz/2;
                i2 = find(ctrdepth<=obj.depth_pr(k),1,'last');
                newvalue(:,:,k) = squeeze(value(:,:,i2));
                i1 = i2+1;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use to expand values from layers onto grid
        function newvalue = expand(obj,value)
            nz = obj.grid.nzEarth + obj.grid.nzCrust;
            if obj.nlayers == 0 || isempty(obj.nlayers)
                error('Unable to expand from layers onto grid: no layers defined');
            elseif obj.nlayers ~= size(value,3)
                error('Value size doesn''t match the number of layers');
            elseif nz == 0
                error('Unable to expand from layers onto grid: no grid defined');
            end
            newvalue = zeros(size(value,1),size(value,2),nz);
            i1 = 1;
            for k = 1:obj.nlayers
                i2 = find(obj.grid.depth<=obj.depth_pr(k),1,'last');
                for ii = i1:i2
                    newvalue(:,:,ii) = squeeze(value(:,:,k));
                end
                i1 = i2+1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use to map from geomagnetic to geographic
        function newobj = gm2gg(obj,clong,res)

            if nargin < 2
                clong = 180;
            end            
            if nargin < 3
                res(2) = min(obj.grid.dlon);
                res(1) = min(obj.grid.dlat);
            end            
            newobj = obj;
            if strcmp(obj.coords,'geographic')
                disp('Nothing to convert');
                return
            end    
            warning off;
            if ~isempty(obj.v)
                [newobj.v,lon_gg,lat_gg] = ...
                        InterpGrid(obj.v,obj.grid.lon,obj.grid.lat,clong,res,0,0);
                disp('3D model converted to geographic');
            end
            if ~isempty(obj.v_pr)
                [newobj.v_pr,lon_gg,lat_gg] = InterpGrid(obj.v_pr,obj.grid.lon,obj.grid.lat,clong,res,0,0);
                disp('Prior model converted to geographic');
            end
            newobj.coords = 'geographic';
            newobj.grid.lon = lon_gg;
            newobj.grid.lat = lat_gg;
            newobj.grid.nlon = length(lon_gg);
            newobj.grid.nlat = length(lat_gg);
            newobj.grid.dlon = res(2)*ones(1,newobj.grid.nlon);
            newobj.grid.dlat = res(1)*ones(1,newobj.grid.nlat);
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use to map from geographic to geomagnetic
        % Set the poles to layer average if undefined
        function newobj = gg2gm(obj,clong,res)

            if nargin < 2
                clong = 180;
            end            
            if nargin < 3
                res(2) = min(obj.grid.dlon);
                res(1) = min(obj.grid.dlat);
            end            
            newobj = obj;
            if strfind(obj.coords,'geomagnetic')
                disp('Nothing to convert');
                return
            end 
            lat = obj.grid.lat; % latitude always sorted
            lon = obj.grid.lon;
            value = obj.v;
            poles = repmat(squeeze(mean(mean(obj.v))).',[length(lon) 1]);
            if min(lat) > -90 
                lat = [-90 lat];
                value(2:end+1,:,:) = value;
                value(1,:,:) = poles;
            end
            if max(lat) < 90
                lat = [lat 90];
                value(end+1,:,:) = poles;
            end
            warning off;
            if ~isempty(obj.v)
                [newobj.v,lon_gm,lat_gm] = ...
                        InterpGrid(value,lon,lat,clong,res,0,1);
                disp('3D model converted to geomagnetic');
            end
            lon = newobj.grid.lon; % fix the remaining NaNs (South Pole)
            poles = repmat(squeeze(nanmean(nanmean(newobj.v))).',[length(lon) 1]);
            [i,j,k]=ind2sub(size(newobj.v),find(isnan(newobj.v)));
            for ind=1:length(i)
                newobj.v(i(ind),j(ind),k(ind)) = poles(j(ind),k(ind));
            end
            if ~isempty(obj.v_pr)
                [newobj.v_pr,lon_gm,lat_gm] = InterpGrid(obj.v_pr,obj.grid.lon,obj.grid.lat,clong,res,0,1);
                disp('Prior model converted to geomagnetic');
            end
            %newobj = fix_nans(newobj); % fix the remaining NaNs (South Pole)
            newobj.coords = 'geomagnetic';
            newobj.grid.lon = lon_gm;
            newobj.grid.lat = lat_gm;
            newobj.grid.nlon = length(lon_gm);
            newobj.grid.nlat = length(lat_gm);
            newobj.grid.dlon = res(2)*ones(1,newobj.grid.nlon);
            newobj.grid.dlat = res(1)*ones(1,newobj.grid.nlat);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fix NaN values; replace them with the nanmeans of the layers;
        % typically used for the South Pole which always becomes NaN
        % in conversion from GG back to GM. This is so general so that
        % in the future we could easily change it to something else
        % that varies with latitude, but never with longitude.
        function [obj] = fix_nans(obj)
            lon = obj.grid.lon;
            if ~isempty(obj.v)
                vmean = repmat(squeeze(nanmean(nanmean(obj.v))).',[length(lon) 1]);
                [i,j,k]=ind2sub(size(obj.v),find(isnan(obj.v)));
                for ind=1:length(i)
                    obj.v(i(ind),j(ind),k(ind)) = vmean(j(ind),k(ind));
                end
            end
            if ~isempty(obj.v_pr)
                vmean = repmat(squeeze(nanmean(nanmean(obj.v))).',[length(lon) 1]);
                [i,j,k]=ind2sub(size(obj.v),find(isnan(obj.v)));
                for ind=1:length(i)
                    obj.v(i(ind),j(ind),k(ind)) = vmean(j(ind),k(ind));
                end
            end
        end
    end
    
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = read(cfile,format)
            %
            % Usage:  [obj] = read(cfile,format)
            %
            % Reads in any globalemmodel object in formats NetCDF or global
            % If global, will also ask for the grid and the prior files.
                        
            if ~exist(cfile,'file')
                error(['Model file ' cfile ' not found']);
            end
            
            switch lower(format)
                case 'global'            
                    obj_geomagnetic = globalemmodel;
                    obj_geomagnetic = obj_geomagnetic.read_rho(cfile);
                    obj_geomagnetic = obj_geomagnetic.read_grid;
                    obj_geomagnetic = obj_geomagnetic.pr;
                    obj_geomagnetic.v = obj_geomagnetic.squeeze(obj_geomagnetic.v);
                    %obj.plot(4);
                    %obj.v_pr = obj.expand(obj.v_pr);
                    obj_geomagnetic.coords = 'geomagnetic';
                    obj = obj_geomagnetic.gm2gg;
                case 'netcdf'
                    obj = globalemmodel;
                    obj_geographic = llmodel.read(cfile,'netcdf',1);
                    obj = obj_geographic.copy(obj);
                otherwise
                    error(['Can''t read lat/lon model format ' format ': method unknown']);
            end

        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = llmodel2global(llobj)
            % just copy everything from the general llmodel to global
            obj = globalemmodel;
            obj.fileName = llobj.fileName;
            obj.fileHeader = llobj.fileHeader;
            obj.geospatialInfo = llobj.geospatialInfo;
            obj.grid = llobj.grid;
            obj.limits = llobj.limits;
            obj.location = llobj.location;
            obj.modelName = llobj.modelName;
            obj.modelType = llobj.modelType;
            obj.modelUnits = llobj.modelUnits;
            obj.displayName = llobj.displayName;
            obj.paramType = llobj.paramType;
            obj.v = llobj.v;
            obj.v_ref = llobj.v_ref;
            obj.AirCond = llobj.AirCond;
            obj.SeaWaterCond = llobj.SeaWaterCond;
            obj.coords = llobj.coords;
        end
        
    end
end