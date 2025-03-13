classdef mtperiod < latlontools
    %   read, write and redefine errors on MT data
    %   uses lat/lon conversions and topography
    %   (c) Anna Kelbert, April 2014
    % 
    %   by design, this should always contain a single TF type
    %   (for multiple types, use an array of these elements)
    %
    %   ModEM currently supports the following types:
    %
    %   Full_Impedance
    %   Off_Diagonal_Impedance
    %   Full_Vertical_Components
    %   Full_Interstation_TF
    %   Off_Diagonal_Rho_Phase
    %   Phase_Tensor
    %
    %   ideally, we would want to match ModEM types to those allowed
    %   by the XML definition (and never mix units together into one)...
    %   but that's something to be discussed.
    
    properties
        T           % period in seconds
        Cmplx = 1   % real or complex
        units = '[mV/km]/[nT]'  % use ImpUnits to convert
        signConvention = 1  % 1 or -1
        nComp       % number of *real* transfer function components
        nSite       % number of data sites for this period and data type
        siteLoc     % [nSite x 3] x,y,z site locations in meters
        siteChar    % site names
        TF          % [nSite x nComp/2] if complex; else [nSite x nComp]
        TFerr       % [nSite x nComp/2] if complex; else [nSite x nComp]
                    % TFerr is *standard error* not the variance!
                    % TFerr is sqrt(TFVar/2) for complex values
        orient = 0  % orientation relative to the geographic North
        lat         % nSite vector of latitudes
        lon         % nSite vector of longitudes
        compChar    % [nComp/2 x 3] TF component names
        primaryCoords = 'latlon' % latlon or xy
        txType = 'MT' % defines the transmitter type for writing to file
        type        % determines which transfer function this is!!!
    end
    
    properties (Constant)
        nanvalue = 999999;
    end
    
    properties (SetAccess = protected)
        % x=0,y=0,z=0 in lat/lon coordinates
        % when updated, automatically recompute siteLoc
        % (or lat/lon, depending on which are primary)
        origin
        % mstruct is as used in Matlab Mapping Toolbox; defines projection
        mstruct
        % if we choose to remove topography, store all z values here;
        % to reinsert, use this variable
        topography
        % if we choose to impose error floors, store the true errors here;
        % to reinsert, use this variable
        dataErrors
    end
    
    methods
        function [obj] = mtperiod(varargin)
            %   class constructor
            %   [obj] =
            %   mtperiod(type,nSite,primaryCoords,signConvention)
            obj.origin = zeros(1,3);
            obj.topography = [];
            obj.dataErrors = [];
            
            if nargin == 0
                return
            end 
            
            if nargin>0
                obj.type = varargin{1};
                if strcmp(obj.type,'Full_Impedance')
                    obj.nComp = 8;
                    obj.Cmplx = 1;
                    obj.units = '[mV/km]/[nT]';
                    obj.compChar = ['ZXX';'ZXY';'ZYX';'ZYY'];
                elseif strcmp(obj.type,'Off_Diagonal_Impedance')
                    obj.nComp = 4;
                    obj.Cmplx = 1;
                    obj.units = '[mV/km]/[nT]';
                    obj.compChar = ['ZXY';'ZYX'];
                elseif strcmp(obj.type,'Full_Vertical_Components')
                    obj.nComp = 4;
                    obj.Cmplx = 1;
                    obj.units = '[]';
                    obj.compChar = ['TX ';'TY '];
                elseif strcmp(obj.type,'Full_Interstation_TF')
                    obj.nComp = 8;
                    obj.Cmplx = 1;
                    obj.units = '[]';
                    obj.compChar = ['MXX';'MXY';'MYX';'MYY'];
                elseif strcmp(obj.type,'Off_Diagonal_Rho_Phase')
                    obj.nComp = 4;
                    obj.Cmplx = 0;
                    obj.units = '[]';
                    obj.compChar = ['RHOXY';'PHSXY';'RHOYX';'PHSYX'];
                elseif strcmp(obj.type,'Phase_Tensor')
                    obj.nComp = 4;
                    obj.Cmplx = 0;
                    obj.units = '[]';
                    obj.compChar = ['PTXX';'PTXY';'PTYX';'PTYY'];
                else
                    error('Unable to initialize MT period: unknown data type');
                end
            end
          
            if nargin == 1
                obj.nSite = 0;
                return
            end
           
            if isstruct(varargin{2})
                % simply converting a structure to an object
                % - needed for backward compatibility with old codes
                datastruct = varargin{2};
                obj.T = datastruct.T;
                obj.Cmplx = datastruct.Cmplx;
                obj.units = datastruct.units;
                obj.signConvention = datastruct.signConvention;
                obj.TF = datastruct.Z;
                obj.TFerr = datastruct.Zerr;
                obj.siteLoc = datastruct.siteLoc;
                obj.siteChar = datastruct.siteChar;
                obj.nSite = length(datastruct.siteChar);
                obj.lat = datastruct.lat;
                obj.lon = datastruct.lon;
                obj.origin = datastruct.origin;
                obj.orient = datastruct.orient;
            else
                obj.nSite = varargin{2};
                if obj.nSite < 0
                    error('Number of data sites must be positive or zero')
                end
                
                % if data already initialized, skip this step
                if isempty(obj.TF)
                    if obj.Cmplx == 1
                        obj.TF(1:obj.nSite,1:obj.nComp/2) = NaN + 1i*NaN;
                        obj.TFerr(1:obj.nSite,1:obj.nComp/2) = NaN;
                    else
                        obj.TF(1:obj.nSite,1:obj.nComp) = NaN;
                        obj.TFerr(1:obj.nSite,1:obj.nComp) = NaN;
                    end
                    
                    obj.siteLoc(1:obj.nSite,1:3) = 0;
                    
                    obj.lat(1:obj.nSite) = 0;
                    obj.lon(1:obj.nSite) = 0;
                end
            end
            
            if nargin > 2
                obj.primaryCoords = varargin{3};
            end
            
            if nargin > 3
                obj.signConvention = varargin{4};
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = initData(obj,nSite,dataValue,dataError)
            % obj = initData(obj,nSite,dataValue,dataError)
            %
            % Initializes the data arrays TF and TFerr
            
            if nargin<2
                nSite = obj.nSite;
            end
            
            if nargin<3
                if obj.Cmplx == 1
                    dataValue = 0.0 + 1i*0.0;
                else
                    dataValue = 0.0;
                end
            end
            
            if nargin<4
                dataError = obj.nanvalue;
            end
            
            if isempty(obj.TF)
                if obj.Cmplx == 1
                    obj.TF(1:nSite,1:obj.nComp/2) = dataValue;
                    obj.TFerr(1:nSite,1:obj.nComp/2) = dataError;
                else
                    obj.TF(1:nSite,1:obj.nComp) = dataValue;
                    obj.TFerr(1:nSite,1:obj.nComp) = dataError;
                end
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = merge(obj1,obj2)
            % obj = merge(obj1,obj2)
            %
            % Merges the sites and data arrays Z and Zerr
            %
            % Currently does not locate and remove potential duplicates!
            
            if (abs(obj1.T - obj2.T)/obj1.T) > 1e-6
                warning('Merging two mtperiod objects that have different periods');
            end

            if ~strcmp(obj1.type,obj2.type)
                warning('Merging two mtperiod objects whose data types that are not consistent');
            end

            if ~strcmp(obj1.units,obj2.units)
                warning('Merging two mtperiod objects whose units are not consistent');
            end

            if obj1.signConvention ~= obj2.signConvention
                warning('Merging two mtperiod objects whose sign conventions are not consistent');
            end

            if obj1.origin ~= obj2.origin
                warning('Merging two mtperiod objects whose origins are not consistent!!!');
            end
            
            names1 = obj1.siteChar;
            names2 = obj2.siteChar;
            nchar = max(length(names1(1,:)),length(names2(1,:)));
            nchar1 = length(names1(1,:));
            nchar2 = length(names2(1,:));
            nsites = length(names1) + length(names2);
            nsites1 = length(names1);
            
            obj = obj1;
            obj.nSite = nsites;
            obj.siteChar = repmat(' ',[nsites nchar]);
            obj.lat = [obj1.lat; obj2.lat];
            obj.lon = [obj1.lon; obj2.lon];
            obj.siteChar(1:nsites1,1:nchar1) = char(names1);
            obj.siteChar(nsites1+1:nsites,1:nchar2) = char(names2);
            obj.siteLoc = [obj1.siteLoc; obj2.siteLoc];
            obj.TF = [obj1.Z; obj2.Z];
            obj.TFerr = [obj1.Zerr; obj2.Zerr];

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = gridCells(obj,grid,landORsea,mstruct)
            % obj = gridCells(obj,grid,landORsea,mstruct)
            %
            % A utility to create a synthetic site at the center of every
            % grid cell at the surface of the Earth, no elevation.
            % Use addTopography function to create non-zero elevation.
            % Set landORsea = 'land' to only get continental locations.
            % Set landORsea = 'sea' for seafloor sites.
            % Set landORsea = 'all' for all sites.
            %
            % If grid is of xygrid type, mstruct is required
            % to compute the lats & lons; otherwise not required.
            % See setOrigin() for details.
            
            if isa(grid,'llgrid')
                mygrid = grid;
                if nargin < 4
                    lat0 = mean(mygrid.lat);
                    lon0 = mean(mygrid.lon);
                end
                mstruct.origin = [lat0,lon0,0];
            elseif isa(grid,'xygrid')
                if nargin < 4
                    error('Usage: dat = dat.gridCells(grid,''land'',mstruct)');
                end
                lat0 = mstruct.origin(1);
                lon0 = mstruct.origin(2);
                mygrid = llgrid(grid,lat0,lon0);
            else
                % if this isn't a grid object, can't do much here
                error('Please specify an xygrid or llgrid on which to define the MT data');
            end
 
            if nargin < 3
                landORsea = 'all';
            elseif contains(landORsea,'land') || contains(landORsea,'sea')
                [melev,mlon,mlat] = m_etopo2([grid.limits.lonmin grid.limits.lonmax grid.limits.latmin grid.limits.latmax]);
                [LAT,LON] = meshgrid(mygrid.lat,mygrid.lon);
                H = interp2(mlon,mlat,melev,LON,LAT);
                if sum(sum(isnan(H))) == size(H,1)*size(H,2)
                    warning('No land/sea found. Check your longitude conventions in mtperiod.gridCells');
                end
            end
            
            % Count cells that matter, then initialize
            nCells = 0;
            for i = 1+mygrid.nPad:mygrid.nlat-mygrid.nPad
                for j = 1+mygrid.nPad:mygrid.nlon-mygrid.nPad
                    if contains(landORsea,'land')
                        % skip seafloor sites
                        if H(j,i)<0
                            continue
                        end
                    elseif contains(landORsea,'sea')
                        % skip land sites
                        if H(j,i)>=0
                            continue
                        end
                    end
                    nCells = nCells+1;
                end
            end
            
            % Initialize for this number of sites
            if isempty(obj.TF) || obj.nSite~=nCells
                obj.nSite = nCells;
                
                if obj.Cmplx == 1
                    obj.TF(1:obj.nSite,1:obj.nComp/2) = NaN;
                    obj.TFerr(1:obj.nSite,1:obj.nComp/2) = NaN;
                else
                    obj.TF(1:obj.nSite,1:obj.nComp) = NaN;
                    obj.TFerr(1:obj.nSite,1:obj.nComp) = NaN;
                end
                
                obj.siteChar = cell(1,obj.nSite);
                obj.siteLoc(1:obj.nSite,1:3) = 0;
                
                obj.lat(1:obj.nSite) = 0;
                obj.lon(1:obj.nSite) = 0;
            end
               
            % Now define these sites at cell centers
            k = 0;
            for i = 1+mygrid.nPad:mygrid.nlat-mygrid.nPad
                for j = 1+mygrid.nPad:mygrid.nlon-mygrid.nPad
                    if contains(landORsea,'land')
                        % skip seafloor sites
                        if H(j,i)<0
                            continue
                        end
                    elseif contains(landORsea,'sea')
                        % skip land sites
                        if H(j,i)>=0
                            continue
                        end
                    end
                    k = k+1;
                    sitename = [sprintf('%03d',i-mygrid.nPad) '-' sprintf('%03d',j-mygrid.nPad)];
                    obj.lat(k) = mygrid.lat(i);
                    obj.lon(k) = mygrid.lon(j);
                    if isa(grid,'llgrid')
                        obj.siteLoc(k,:) = [grid.lat(i) grid.lon(j) 0];
                        obj.primaryCoords = 'latlon';
                    elseif isa(grid,'xygrid')
                        xctr = grid.xctr;
                        yctr = grid.yctr;
                        if strcmp(grid.units,'km')
                            xctr = 1e3 * xctr;
                            yctr = 1e3 * yctr;
                        end
                        obj.siteLoc(k,:) = [xctr(i) yctr(j) 0];
                        obj.primaryCoords = 'xy';
                    end
                    obj.siteChar{k} = char(sitename);
                end
            end
            
            obj.lat = obj.lat';
            obj.lon = obj.lon';
            
            if nargin > 3
                obj = obj.setOrigin(mstruct);
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = gridNodes(obj,grid,landORsea,mstruct)
            % obj = gridNodes(obj,grid,landORsea,mstruct)
            %
            % A utility to create a synthetic site at every
            % grid node at the surface of the Earth, no elevation.
            % Use addTopography function to create non-zero elevation.
            % Set landORsea = 'land' to only get continental locations.
            % Set landORsea = 'sea' for seafloor sites.
            % Set landORsea = 'all' for all sites.
            %
            % If you are using topography, it is ill-advised to place
            % synthetic sites at cell nodes. Use gridCells, instead.
            %
            % If grid is of xygrid type, mstruct is required
            % to compute the lats & lons; otherwise not required.
            % See setOrigin() for details.
          
           
            if isa(grid,'llgrid')
                mygrid = grid;
                if nargin < 4
                    lat0 = mean(mygrid.lat);
                    lon0 = mean(mygrid.lon);
                end
                mstruct.origin = [lat0,lon0,0];
            elseif isa(grid,'xygrid')
                if nargin < 4
                    warning('Usage: dat = dat.gridNodes(grid,''land'',mstruct). Cannot compute lats & lons without mstruct');
                else
                    mygrid = llgrid(grid,mstruct);
                end
            else
                % if this isn't a grid object, can't do much here
                error('Please specify an xygrid or llgrid on which to define the MT data');
            end
            
            if nargin < 3
                landORsea = 'all';
            elseif contains(landORsea,'land') || contains(landORsea,'sea')
                [melev,mlon,mlat] = m_etopo2([grid.limits.lonmin grid.limits.lonmax grid.limits.latmin grid.limits.latmax]);
                [LAT,LON] = meshgrid(mygrid.lat,mygrid.lon);
                H = interp2(mlon,mlat,melev,LON,LAT);
                if sum(sum(isnan(H))) == size(H,1)*size(H,2)
                    warning('No land/sea found. Check your longitude conventions in mtperiod.gridNodes');
                end
            end
            
            % Count nodes that matter, then initialize
            nNodes = 0;
            for i = 1+mygrid.nPad:mygrid.nlat+1-mygrid.nPad
                for j = 1+mygrid.nPad:mygrid.nlon+1-mygrid.nPad
                    if contains(landORsea,'land')
                        % skip seafloor sites
                        if H(j,i)<0
                            continue
                        end
                    elseif contains(landORsea,'sea')
                        % skip land sites
                        if H(j,i)>=0
                            continue
                        end
                    end
                    nNodes = nNodes+1;
                end
            end
            
            % Initialize for this number of sites
            if isempty(obj.TF) || obj.nSite~=nNodes
                obj.nSite = nNodes;
                
                if obj.Cmplx == 1
                    obj.TF(1:obj.nSite,1:obj.nComp/2) = NaN;
                    obj.TFerr(1:obj.nSite,1:obj.nComp/2) = NaN;
                else
                    obj.TF(1:obj.nSite,1:obj.nComp) = NaN;
                    obj.TFerr(1:obj.nSite,1:obj.nComp) = NaN;
                end
                
                obj.siteChar = cell(1,obj.nSite);
                obj.siteLoc(1:obj.nSite,1:3) = 0;
                
                obj.lat(1:obj.nSite) = 0;
                obj.lon(1:obj.nSite) = 0;
            end
               
            % Now define these sites at cell nodes
            k = 0;
            for i = 1+mygrid.nPad:mygrid.nlat+1-mygrid.nPad
                for j = 1+mygrid.nPad:mygrid.nlon+1-mygrid.nPad
                    if contains(landORsea,'land')
                        % skip seafloor sites
                        if H(j,i)<0
                            continue
                        end
                    elseif contains(landORsea,'sea')
                        % skip land sites
                        if H(j,i)>=0
                            continue
                        end
                    end
                    k = k+1;
                    sitename = [sprintf('%03d',i-mygrid.nPad) '-' sprintf('%03d',j-mygrid.nPad)];
                    obj.lat(k) = mygrid.lat(i);
                    obj.lon(k) = mygrid.lon(j);
                    if isa(grid,'llgrid')
                        obj.siteLoc(k,:) = [grid.lat(i) grid.lon(j) 0];
                        obj.primaryCoords = 'latlon';
                    elseif isa(grid,'xygrid')
                        x = grid.x;
                        y = grid.y;
                        if strcmp(grid.units,'km')
                            x = 1e3 * x;
                            y = 1e3 * y;
                        end
                        obj.siteLoc(k,:) = [x(i) y(j) 0];
                        obj.primaryCoords = 'xy';
                    end
                    obj.siteChar{k} = char(sitename);
                end
            end

            obj.lat = obj.lat';
            obj.lon = obj.lon';

            if nargin > 3
                obj = obj.setOrigin(mstruct);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setOrigin(obj,mstruct)
            % (re)define the geographic origin and correspondingly
            % (re)compute either lat/lon or x/y site coordinates
            %
            % NOTE: z pointing down! In data files, topography should be
            % +ve above the "origin". But currently I think z is depth
            % (since z is pointing down like in the model file), so in
            % current configuration, -ve values should define topography,
            % and z0 should be a *negative* value at the sea level to
            % indicate the distance from the "top of the mountains" to
            % ground (and the same value should be used in the model file).
            % All this is unnecessarily confusing and needs to be changed
            % in the program, probably to reverse the direction of z in the
            % data file?..
            % Need to figure it all out. What's currently implemented here
            % with respect to z will not work correctly. 
            % This is just a stub to be fixed.
            % In Matlab, topography in the data file should help create the
            % model file, anyway.
            %
            % Soooo easy to mix up km and meters! If you specify the geoid
            % in meters, then no need to multiply by a thousand to compute
            % site coordinates in meters. Otherwise, get huge distances
            % that make no sense! A crude fix is to check the size of the
            % geoid. Doing that here.
            %
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
            % If using UTM, note that the X/Y coordinates need to include
            % the binding Northing and Easting for things to work.  
            
            lat0 = mstruct.origin(1);
            lon0 = mstruct.origin(2);
            if length(mstruct.origin) > 2
                z0 = mstruct.origin(3);
            else
                z0 = 0;
            end
            
            if isfield(mstruct,'geoid')
                if mstruct.geoid(1) > 6300 && mstruct.geoid(1) < 6400
                    disp('Input map projection is in km. Converting to meters!');
                    km2m = 1000;
                elseif mstruct.geoid(1) > 6300000 && mstruct.geoid(1) < 6400000
                    disp('Input map projection is in meters. No conversion needed.');
                    km2m = 1;
                else
                    warning('Unable to resolve the km / meters ambiguity in setOrigin. Check all distances!!!')
                    km2m = 1;
                end
            else
                km2m = 1000;
            end
                
            if strcmp(mstruct.mapprojection,'latlon2xy') ...
                    || strcmp(mstruct.mapprojection,'xy2latlon')
                % using the homemade projection
                obj.origin = [lat0,lon0,z0];
            else
                % using the matlab toolbox
                try
                    mstruct = defaultm(mstruct);
                catch
                    try
                        mstruct = defaultm(mstruct.mapprojection);
                        mstruct.origin = [lat0,lon0,z0];
                        mstruct = defaultm(mstruct);
                    catch
                        error(['Unknown map projection: ',mstruct.mapprojection]);
                    end
                end
                obj.origin = mstruct.origin;
            end
            
            if strcmp(obj.primaryCoords,'latlon')
                % recompute x/y - be careful with array shapes
                sitelatlon = [obj.lat'; obj.lon'];
                if size(sitelatlon,1) ~= 2
                    sitelatlon = reshape(sitelatlon,length(obj.lat),2);
                    sitelatlon = sitelatlon';
                end
                switch strtrim(mstruct.mapprojection)
                    case 'latlon2xy'
                        sitexy = latlontools.ll2xy(sitelatlon,lat0,lon0,z0);
                    case 'utm'
                        [sitexy(2,:),sitexy(1,:)] = mfwdtran(mstruct,sitelatlon(1,:),sitelatlon(2,:));
                    otherwise
                    %case {'eqdcylin','eqacylin','lambertstd','eqaazim'}
                        [sitexy(2,:),sitexy(1,:)] = projfwd(mstruct,sitelatlon(1,:),sitelatlon(2,:));
                end
                obj.siteLoc(:,1) = km2m*sitexy(1,:)';
                obj.siteLoc(:,2) = km2m*sitexy(2,:)';
                sitedepth = obj.siteLoc(:,3)';
                obj.siteLoc(:,3) = z0 + sitedepth';
            elseif strcmp(obj.primaryCoords,'xy')
                % recompute lat/lon
                switch strtrim(mstruct.mapprojection)
                    case 'latlon2xy'
                        sitelatlon = latlontools.xy2ll(obj.siteLoc'/km2m,lat0,lon0,z0);
                    case  'utm'
                        [sitelatlon(1,:),sitelatlon(2,:)]= minvtran(mstruct,obj.siteLoc(:,2)'/km2m,obj.siteLoc(:,1)'/km2m);
                    otherwise
                    %case {'eqdcylin','edacylin','lambertstd','eqaazim'}
                        [sitelatlon(1,:),sitelatlon(2,:)]= projinv(mstruct,obj.siteLoc(:,2)'/km2m,obj.siteLoc(:,1)'/km2m); 
                end
                obj.lat = sitelatlon(1,:)';
                obj.lon = sitelatlon(2,:)';
                sitedepth = obj.siteLoc(:,3)';
                obj.siteLoc(:,3) = z0 + sitedepth';
            end
            
            % store projection details for completeness
            obj.mstruct = mstruct;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = changeDataType(obj,newdatatype)
            % allows to change to new data type; this will work to
            % translate between Full_Impedance, Off_Diagonal_Impedance, and
            % the Off_Diagonal_Rho_Phase. Other options may be added.
            
            if strcmp(obj.type,newdatatype)
                %disp('No type conversion is needed. Command ignored.');
                return;
            end
                        
            newobj = mtperiod(newdatatype,obj.nSite,obj.primaryCoords,obj.signConvention);
            if ~strcmp(obj.units,newobj.units)
                newobj = newobj.changeUnits(obj.units);
            end
            newobj.T = obj.T;
            newobj.nSite = obj.nSite;
            newobj.siteLoc  = obj.siteLoc;
            newobj.siteChar  = obj.siteChar;
            newobj.orient = obj.orient;
            newobj.lat = obj.lat;
            newobj.lon = obj.lon;   
            
            Z = NaN*zeros(2,2,obj.nSite) + 1i*NaN*zeros(2,2,obj.nSite);
            Zstd = NaN*zeros(2,2,obj.nSite);
            
            if strcmp(obj.type,'Full_Impedance')
                Z(1,1,:) = obj.TF(:,1);         
                Z(1,2,:) = obj.TF(:,2); 
                Z(2,1,:) = obj.TF(:,3); 
                Z(2,2,:) = obj.TF(:,4); 
                Zstd(1,1,:) = obj.TFerr(:,1);         
                Zstd(1,2,:) = obj.TFerr(:,2); 
                Zstd(2,1,:) = obj.TFerr(:,3); 
                Zstd(2,2,:) = obj.TFerr(:,4); 
            elseif strcmp(obj.type,'Off_Diagonal_Impedance')
                Z(1,1,:) = NaN + 1i*NaN;         
                Z(1,2,:) = obj.TF(:,1); 
                Z(2,1,:) = obj.TF(:,2); 
                Z(2,2,:) = NaN + 1i*NaN; 
                Zstd(1,1,:) = NaN;        
                Zstd(1,2,:) = obj.TFerr(:,1); 
                Zstd(2,1,:) = obj.TFerr(:,2); 
                Zstd(2,2,:) = NaN;
            elseif strcmp(obj.type,'Off_Diagonal_Rho_Phase')
                for i = 1:obj.nSite
                    apres.xy = obj.TF(i,1);
                    apres.xy_se = obj.TFerr(i,1);
                    phase.xy = obj.TF(i,2);
                    phase.xy_se = obj.TFerr(i,2);
                    apres.yx = obj.TF(i,3);
                    apres.yx_se = obj.TFerr(i,3);
                    phase.yx = obj.TF(i,4);
                    phase.yx_se = obj.TFerr(i,4);
                    [Z1,Z1std] = mttf.apres2imp(obj.T,apres,phase);
                    Z(:,:,i) = Z1;
                    Zstd(:,:,i) = Z1std;
                end
            end

            if strcmp(newdatatype,'Full_Impedance')
                newobj.TF(:,1) = Z(1,1,:);         
                newobj.TF(:,2) = Z(1,2,:); 
                newobj.TF(:,3) = Z(2,1,:); 
                newobj.TF(:,4) = Z(2,2,:); 
                newobj.TFerr(:,1) = Zstd(1,1,:);         
                newobj.TFerr(:,2) = Zstd(1,2,:); 
                newobj.TFerr(:,3) = Zstd(2,1,:); 
                newobj.TFerr(:,4) = Zstd(2,2,:); 
            elseif strcmp(newdatatype,'Off_Diagonal_Impedance')
                newobj.TF(:,1) = Z(1,2,:);         
                newobj.TF(:,2) = Z(2,1,:); 
                newobj.TFerr(:,1) = Zstd(1,2,:);         
                newobj.TFerr(:,2) = Zstd(2,1,:); 
            elseif strcmp(newdatatype,'Off_Diagonal_Rho_Phase')
                for i = 1:obj.nSite
                    Z1 = squeeze(Z(:,:,i));
                    Z1std = squeeze(Zstd(:,:,i));
                    [apres,phase] = mttf.imp2apres(obj.T,Z1,Z1std);
                    newobj.TF(i,1) = apres.xy;
                    newobj.TFerr(i,1) = apres.xy_se;
                    newobj.TF(i,2) = phase.xy;
                    newobj.TFerr(i,2) = phase.xy_se;
                    newobj.TF(i,3) = apres.yx;
                    newobj.TFerr(i,3) = apres.yx_se;
                    newobj.TF(i,4) = phase.yx;
                    newobj.TFerr(i,4) = phase.yx_se;
                end
            end
            
            obj = newobj;

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = changeUnits(obj,newunits)
            % allows to change to new units; this only works correctly if
            % this vector contains one data type
            % unit options for impedances:
            % 1) SI units for E/B: [V/m]/[T] (used internally in ModEM code)
            % 2) practical units for E/B: [mV/km]/[nT]
            % 3) SI units for E/H: [V/m]/[A/m] = Ohm
                        
            SI_factor = ImpUnits(obj.units,newunits);
            obj.TF = SI_factor * obj.TF;
            obj.TFerr = SI_factor * obj.TFerr;
            obj.units = newunits;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = changeSignConvention(obj)
            % allows to flip the sign convention from -ve to +ve and back
                        
            if obj.Cmplx
                obj.TF = conj(obj.TF);
                obj.TFerr = conj(obj.TFerr);
                obj.signConvention = - obj.signConvention;
            else
                disp('Cannot change the sign convention for real data values. Command ignored.')
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = fillNaNs(obj,value)
            % obj = fillNaNs(obj,value)
            %
            % replace all NaNs with obj.nanvalue for output
            % for a forward modeling template, use zeros for TFs
            % then, call fillNaNs(obj,0.0)
            % always use obj.nanvalue for NaN error bars
            
            if nargin<=1
                value = obj.nanvalue;
            end
            
            ii = isnan(obj.TF);
            if obj.Cmplx
                obj.TF(ii) = value + 1i*value;
            else
                obj.TF(ii) = value;
            end           
            
            ii = isnan(obj.TFerr);
            obj.TFerr(ii) = obj.nanvalue;            
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = addTopography(obj)
            % replace existing z coordinates with the stored "topography"
            % (in the future, could extract it from the model, instead)
            obj.siteLoc(:,3) = obj.topography;
           
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = removeTopography(obj)
            % set all z coordinates to zero and store them in "topography"
            obj.topography = obj.siteLoc(:,3);
            obj.siteLoc(:,3) = 0;
                       
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setErrorFloor(obj,relErr,site)
            %  Usage : obj = setErrorFloor(obj,relErr,site);
            %
            %  Based on the external function setErrorFloor3D. Remember, by
            %  design mtperiod objects only include one data type.
            
            if isempty(obj.dataErrors); obj.dataErrors = obj.TFerr; end

            data = obj.TF;
            err = obj.TFerr;
            sites = obj.siteChar;
            nsite  = size(data,1);
            if nargin > 2
                [~,i] = intersect(sites,site,'rows');
            else
                i = 1:nsite;
            end
            
            if isempty(i)
                error(['Cannot set the error floor for site ' site ': not found']);
            end
            
            if strcmp(obj.type,'Full_Impedance')
                errFloor = relErr * sqrt(abs(data(i,2).*data(i,3)));
                errZ = max(err(i,1:4),repmat(errFloor,1,4));
                obj.TFerr(i,1:4) = errZ;
            elseif strcmp(obj.type,'Off_Diagonal_Impedance')
                errFloor = relErr * sqrt(abs(data(i,1).*data(i,2)));
                errZ = max(err(i,1:2),repmat(errFloor,1,2));
                obj.TFerr(i,1:2) = errZ;
            elseif strcmp(obj.type,'Full_Vertical_Components')
                errFloor = relErr;
                errT = max(err(i,1:2),repmat(errFloor,1,2));
                obj.TFerr(i,1:2) = errT;
            elseif strcmp(obj.type,'Off_Diagonal_Rho_Phase')
                errFloor = relErr(1);
                errRho = max(err(i,1),errFloor);
                obj.TFerr(i,1) = errRho;
                errFloor = relErr(2);
                errPhs = max(err(i,2),errFloor);
                obj.TFerr(i,2) = errPhs;
            else
                error(['Cannot set the error floor for data type ' obj.type ': method not defined yet']);
                
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = set(obj,varargin)
            %  Sets the values of individual or multiple TF, TFerr from Z, Zstd
            %
            %  options:
            %   [obj] = obj.set(iSite,Z);
            %   [obj] = obj.set(iSite,Z,Zstd);
            %
            %  Can be extended to set the values for other transfer functions as
            %  long as the format of the array matches the object "type".
            %  Note that nargin also counts the object itself.
            
            if nargin <= 1
                return
            end 
            
            iSite = varargin{1};
            if iSite > obj.nSite
                error('Cannot set the impedance values - site doesn''t exist');
            end
            
            Z = varargin{2};
            for j = 1:length(iSite)
                for iComp = 1:length(obj.compChar)
                    if strcmp('ZXX',obj.compChar(iComp,:))
                        obj.TF(iSite(j),iComp) = squeeze(Z(1,1,j));
                    elseif strcmp('ZXY',obj.compChar(iComp,:))
                        obj.TF(iSite(j),iComp) = squeeze(Z(1,2,j));
                    elseif strcmp('ZYX',obj.compChar(iComp,:))
                        obj.TF(iSite(j),iComp) = squeeze(Z(2,1,j));
                    elseif strcmp('ZYY',obj.compChar(iComp,:))
                        obj.TF(iSite(j),iComp) = squeeze(Z(2,2,j));
                    end
                end
            end
            
            if nargin>3
                Zstd = varargin{3};
                for j = 1:length(iSite)
                    for iComp = 1:length(obj.compChar)
                        if strcmp('ZXX',obj.compChar(iComp,:))
                            obj.TFerr(iSite(j),iComp) = squeeze(Zstd(1,1,j));
                        elseif strcmp('ZXY',obj.compChar(iComp,:))
                            obj.TFerr(iSite(j),iComp) = squeeze(Zstd(1,2,j));
                        elseif strcmp('ZYX',obj.compChar(iComp,:))
                            obj.TFerr(iSite(j),iComp) = squeeze(Zstd(2,1,j));
                        elseif strcmp('ZYY',obj.compChar(iComp,:))
                            obj.TFerr(iSite(j),iComp) = squeeze(Zstd(2,2,j));
                        end
                    end
                end
            end

        end
            
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function apresobj = apres(obj)
%             % create and apparent resistivity and phase object from an
%             % impedance object; can write the new object to file, or plot
%             
%             apresobj = mtperiod('Off_Diagonal_Rho_Phase',obj.nSite,obj.primaryCoords,obj.signConvention);
%             
%             rad_deg = 57.2958;
%             apres.xy = abs(info.data(:,:,2)).^2;
%             apres.xy_se = info.err(:,:,2).^2; % variance of Z as in EDI file
%             apres.yx = abs(info.data(:,:,3)).^2;
%             apres.yx_se = info.err(:,:,3).^2; % variance of Z as in EDI file
%             %phase.xy = rad_deg*atan((imag(info.data(:,:,2)))./real(info.data(:,:,2)));
%             % for now, using atan2 which takes values (-pi,pi] as in ModEM
%             phase.xy = rad_deg*atan2(imag(info.data(:,:,2)),real(info.data(:,:,2)));
%             phase.xy_se = rad_deg*sqrt(apres.xy_se./apres.xy);
%             %phase.yx = rad_deg*atan((imag(info.data(:,:,3)))./real(info.data(:,:,3)));
%             % for now, using atan2 which takes values (-pi,pi] as in ModEM
%             phase.yx = rad_deg*atan2(imag(info.data(:,:,3)),real(info.data(:,:,3))) + 180.;
%             phase.yx_se = rad_deg*sqrt(apres.yx_se./apres.yx);
%             % rescale apparent resistivity by period
%             for l = 1:length(info.per)
%                 apres.yx(:,l) = apres.yx(:,l)*info.per(l)/5. ;
%                 apres.xy(:,l) = apres.xy(:,l)*info.per(l)/5. ;
%                 apres.yx_se(:,l) = sqrt(apres.yx_se(:,l).*apres.yx(:,l)*info.per(l)*4/5.);
%                 apres.xy_se(:,l) = sqrt(apres.xy_se(:,l).*apres.xy(:,l)*info.per(l)*4/5.);
%             end
%             
%         end
        
    end    
end