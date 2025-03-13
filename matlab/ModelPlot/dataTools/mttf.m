classdef mttf
    %   works with a single MT TF to read, write, combine and plot
    %   (c) Anna Kelbert, Feb 2016 - Mar 2017
    % 
    %   reads in Z-files, XML, EDI, BIRRP
    %
    % Description borrowed from Gary Egbert's EMTF/matlab/ZPLT/Z_in.m :
    % nch = total # of channels ; nche = nch-2 = # of predicted channels ; 
    % nbt = # of bands
    %   (NOTE: First two channels are always the "predictors"
    % z(2,nche*nbt) = complex TFs
    %   NOTE: Z(1,1:nche) corresponds to response to Hx sources for first band,
    %         Z(2,1:nche) is Hy for for first band,
    %         Z(1,nche+1:2*nche) corresponds to response to Hx sources for second band,
    %         Z(2,nche+1:2*nche) is Hy for for second band,  ETC. ...
    % sig_s(2,2,nbt) = complex inverse signal covariance
    % sig_e(nche,nche,nbt) = complex residual error covariance
    % stdec(3) = station coordinates, declination
    % periods(nbt) = periods in seconds
    % orient(2,nch) = orientation (degrees E of geomagnetic N) for each channel

    
    properties
        tfdir       % primary file directory
        tfname      % primary file name (no extensions)
        tfext       % primary file extension
        primaryFormat % options are Z, XML, EDI, BIRRP
        Cmplx = 1   % real or complex
        units = '[mV/km]/[nT]'  % impedance ONLY; use ImpUnits to convert
        signConvention = 1  % 1 or -1
        nComp       % number of *real* transfer function components
        nSite = 1   % number of data sites for this period and data type
        siteLoc = [0 0 0] % [1 x 3] x,y,z site locations in meters (rarely needed)
        siteChar    % site name
        TF          % [2,nche,nper] full transfer functions
        TFVar       % [2,nche,nper] variance of complex TFs; divide by 2 for real/imag
        SIG_S       % sig_s(2,2,nper) = complex inverse signal covariance
        SIG_E       % sig_e(nche,nche,nper) = complex residual error covariance
        T           % periods in seconds
        orthogonal  % true or false
        orient      % orientation of all channels relative to theta0
        theta0      % orient relative to geographic North (= declination in the file)
        lat         % site latitude
        lon         % site longitude
        elev        % site elevation
        compChar    % [nComp/2 x 3] TF component names
        primaryCoords = 'latlon' % latlon or xy
        types = {'Full_Vertical_Components','Full_Impedance'}
        chead       % defines the TF type, e.g. 'Robust Remote Reference'
        predicted = 0 % set to 1 for a different treatment of errorbars
        provenance  % a structure as read from file
        copyright   % a structure as read from file or created in Matlab
        metadata    % a structure as read from file or created in Matlab
    end

    properties (SetAccess = protected)
        % we can only get these from a file... don't change them.
        % needed to write back to a file without losing information.
        ndf        % number of data points per band (NOT frequencies)
        stcor      % station coordinates [lat lon]
        stname     % station name, same as siteChar
        ibandlim   % [2 x nper] freq. band
        level      % [1 x nper] decimation level
        sampRate   % [1 x nper] sampling rate
        Nch        % number of channels
        Nche       % number of output channels
        nbt        % same as nper = length(T)
        chid = {'Hx'  'Hy'  'Hz'  'Ex'  'Ey'}
        sta        % like chid but station code for each channel
    end
    
    methods
        function [obj] = mttf(varargin)
            %   class constructor
            %
            %  Note that the TF is TRANSPOSED relative to the way an
            %  impedance is typically recorded. Done for historical reasons
            %  (compatibility with old Gary Egbert's codes). To avoid 
            %  misleading results, it is necessary to access the impedance 
            %  through the mttf.impedance function rather than directly.
            %  The same TF array can hold tipper or others. 
            %
            %  This is intended to be used with one data type at a time, 
            %  but it allows to mix different types both for backward
            %  compatibility and to simplify future work with mixed
            %  covariance matrices. Defaults to tipper+impedance.
            %  Certainly by all means never mix complex and real data types
            %  together; or data types with different units.
            %
            %  options:
            %   [obj] = mttf;
            %   [obj] = mttf(dataType);
            %   [obj] = mttf(dataType,periods);
            %   [obj] = mttf(Zstruct);
            %   [obj] = mttf.read;
            %   [obj] = mttf.read(filename);
            %  follow with
            %   [obj] = obj.set(T,Z);
            %   [obj] = obj.set(T,Z,Zstd);
            %   [obj] = obj.set(T,Z,Zstd,latlon,theta,siteid);
            %   [obj] = obj.set(Zstruct);
            
            if nargin > 0
                % recursive option that inits from old structures
                if isstruct(varargin{1})
                    Zstruct = varargin{1};
                    obj = mttf;
                    obj = obj.struct2mttf(Zstruct);
                    return;
                elseif isnumeric(varargin{1})
                    obj = mttf('Full_Impedance');
                    obj = obj.set(varargin{:});
                    return;
                else
                    obj.types = varargin{1};
                    if ischar(obj.types)
                        obj.types = {obj.types};
                    end
                end
            end
            
            obj.nComp = 0;
            obj.compChar = [];
            obj.units = '[]';
            obj.chid = {'Hx' 'Hy'};
            for j = 1:length(obj.types)
                [nC,comp,typeunits,ochid] = mttf.components(obj.types{j});
                obj.nComp = obj.nComp + nC;
                obj.compChar = [obj.compChar; comp];
                obj.units = typeunits; % use units for last data type
                obj.chid = [obj.chid ochid];
                obj.Cmplx = mttf.isComplex(obj.types{j});
            end
            
            if obj.nComp == 0
                error('Unable to initialize mttf object: unknown data type');
            end
            
            if obj.Cmplx
                obj.Nche = obj.nComp/2 / 2;
            else
                obj.Nche = obj.nComp / 2;
            end
            obj.Nch = obj.Nche + 2;
            
            if nargin > 1
                periods = varargin{2};
                obj.T = periods;
                obj.TF = NaN*(zeros(2,obj.Nche,length(periods)) + 1i*zeros(2,obj.Nche,length(periods)));
                obj.TFVar = NaN*zeros(2,obj.Nche,length(periods));                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = set(obj,varargin)
            %   Sets the values of the TF, TFVar arrays for an object that
            %  is pre-initialized for a particular data type.
            %
            %  Note that the TF is TRANSPOSED relative to the way an
            %  impedance is typically recorded. Done for historical reasons
            %  (compatibility with old Gary Egbert's codes). To avoid 
            %  misleading results, it is necessary to access the impedance 
            %  through the mttf.impedance function rather than directly.
            %  The same TF array can hold tipper or others.
            %
            %  options:
            %   [obj] = obj.set(T,Z);
            %   [obj] = obj.set(T,Z,Zstd);
            %   [obj] = obj.set(T,Z,Zstd,latlon,metadata);
            %   [obj] = obj.set(filename);
            %
            %  Can similarly set the values for other transfer functions as
            %  long as the format of the array matches the object "type".
            %  Note that nargin also counts the object itself, so 7 total.
            
            if nargin <= 1
                return
            end 
            
            if ischar(varargin{1})
                filename = varargin{1};            
                format = mttf.format(filename);
                dataType = obj.types;
                obj.types = {};
                obj.nComp = 0;
                obj.Nche = 0;
                for j = 1:length(dataType)
                    tmp = mttf.read(filename,format,dataType{j});
                    obj = obj.appendType(tmp);
                end
            elseif isnumeric(varargin{1})
                periods = varargin{1};
                if length(periods) ~= length(obj.T)
                    % initialize again with new periods
                    dataType = obj.types;
                    obj = mttf(dataType,periods);
                else
                    obj.T = periods;
                end
                % now read the data inputs
                pred = 1;
                if nargin > 2
                    Z = varargin{2};
                end
                Zstd = [];
                if nargin > 3
                    Zstd = varargin{3};
                    pred = 0;
                end
                % if periods come first, make them last
                if size(Z,1) == length(periods)
                    Z = shiftdim(Z,1);
                end
                if size(Zstd,1) == length(periods)
                    Zstd = shiftdim(Zstd,1);
                end
                % add an additional ghost dimension to, e.g., tippers
                if size(Z,2) == length(periods)
                    Z = reshape(Z,[1 2 length(periods)]);
                end
                if size(Zstd,2) == length(periods)
                    Zstd = reshape(Zstd,[1 2 length(periods)]);
                end
                % set the TRANSPOSED arrays - after some error checking
                if obj.Nche ~= size(Z,1)
                    error('mttf object data types and the size of input array incompatible');
                elseif nargin > 7
                    error('This number of numeric input arguments to mttf.set is not supported');
                end
                obj.TF = NaN*(zeros(2,obj.Nche,length(periods)) + 1i*zeros(2,obj.Nche,length(periods)));
                obj.TFVar = NaN*zeros(2,obj.Nche,length(periods));
                obj.predicted = pred;
                for i = 1:obj.Nche
                    for j = 1:2
                        obj.TF(j,i,:) = Z(i,j,:);
                        if ~obj.predicted
                            obj.TFVar(j,i,:) = 2*(Zstd(i,j,:).^2);
                        end
                    end
                end
                % read metadata values
                latlon = [];
                if nargin > 4
                    latlon = varargin{4};
                end
                if ~isempty(latlon)
                    obj.lat = latlon(1);
                    obj.lon = latlon(2);
                    if length(latlon) > 2
                        obj.elev = latlon(3);
                    end
                end
                if nargin > 5
                    info = varargin{5};
                    orthogonal_coords = info.ORTHOGONAL;
                    th0 = info.THETA0; 
                    th = info.ORIENT;
                    siteid = info.SITEID;
                else
                    info = [];
                    orthogonal_coords = 0;
                    th0 = [];
                    th = [];
                    siteid = '';
                end
                obj.orthogonal = orthogonal_coords;
                obj.theta0 = th0;
                obj.orient = th;
                %if length(theta) > 1
                %    obj.orient = theta - theta(1);
                %end
                if ~isempty(siteid)
                    obj.siteChar = siteid;
                end
                if ~isempty(info)
                    obj.metadata = info;
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = appendType(obj,newobj)
            % [obj] = appendType(obj,newobj)
            %
            % a utility function to append additional data type to obj
            % takes the units, complexity etc of the new type
            % not much error checking at present
            %
            % related functions: set
            
            if isempty(obj.T)
                obj = newobj;
                return
            elseif length(obj.T) ~= length(newobj.T)
                error('Objects not compatible in appendType');
            end
            
            if norm((obj.orient(1,1:2)+obj.theta0) - (newobj.orient(1,1:2)+newobj.theta0))>1e-6
                error('Orientations of input channels do not match. Will not append type')
            end
            
            if norm(obj.theta0 - newobj.theta0)>1e-6
                newobj = newobj.rotate(obj.theta0);
            end
            
            if obj.signConvention ~= newobj.signConvention
                newobj.TF = conj(newobj.TF);
                newobj.SIG_S = conj(newobj.SIG_S);
                newobj.SIG_E = conj(newobj.SIG_E);
            end
            
            oldobj = obj;
            
            obj.TF(:,end+1:end+newobj.Nche,:) = newobj.TF;
            obj.TFVar(:,end+1:end+newobj.Nche,:) = newobj.TFVar;
            obj.types = [obj.types newobj.types];
            obj.nComp = obj.nComp + newobj.nComp;
            obj.Nche = obj.Nche + newobj.Nche;
            obj.Nch = obj.Nche + 2;
            obj.orient = [obj.orient newobj.orient(:,3:end)];
            obj.chid = [obj.chid newobj.chid{3:end}];
            obj.predicted = newobj.predicted;
            obj.units = newobj.units;
            
            if ~isempty(obj.SIG_E)
                obj.SIG_E = zeros(obj.Nche,obj.Nche) + 1i*zeros(obj.Nche,obj.Nche);
                obj.SIG_E(1:oldobj.Nche,1:oldobj.Nche) = oldobj.SIG_E;
                obj.SIG_E(oldobj.Nche+1:obj.Nche,oldobj.Nche+1:obj.Nche,:) = newobj.SIG_E;
            end                        
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [impulse,A4real] = impulse(obj,delta_t,nbounds,option,lambda,errorbars)
            %
            % Usage:  [impulse,A4real] = impulse(obj,delta_t,nbounds,lambda1,lambda2,errorbars)
            % 
            % Computes the discrete time impulse response according to
            % Gary Egbert's 1992 method published in Geophysics
            % 
            % NOT USING ERROR BARS AT PRESENT: DTIR CODED BUT NOT DEBUGGED
            % WITH ERROR BARS; JUST ADD A LOGICAL SWITCH TO USE THESE            
            %
            % Inputs:
            % obj is an object of class mttf
            % delta_t is sampling of the time series in secs
            % nbounds = [nmin nmax] for DTIR
            % option  = 1 MODELSPACE & minimize in linear frequency
            %         = 2 MODELSPACE & minimize in log frequency
            %         = 3 DATASPACE & minimize in linear frequency
            %         = 4 DATASPACE & minimize in log frequency - Q inverse not coded!
            % lambda is the damping parameter for DTIR
            % errorbars = 1 means extract and use standard errors
            %           = 2 means extract and use full covariance matrices
            %
            % Outputs:
            %  impulse - a structure that contains the DTIR for xx,xy,yx,yy
            %  A4real - 4-component real representer matrix; Zpred = A4real*zn
            
            if nargin < 3
                nmin = -800;
                nmax = 2000;
            else
                nmin = nbounds(1);
                nmax = nbounds(2);
            end
            
            if nargin < 4
                option = 3;
            end
            
            if nargin < 5
                lambda = 1e-12;
            end
            
            if nargin < 6
                errorbars = 0;
            end
            
            % Extract the MT impedance tensor from mttf
            periods = obj.T;
            [Z,Zstd] = obj.impedance;
            
            % Get rid of NaNs (crudely - remove all such periods)
            ind = find(~isnan(Z(1,1,:)) & ~isnan(Zstd(1,1,:))); ind1 = ind;
            ind = find(~isnan(Z(1,2,:)) & ~isnan(Zstd(1,2,:))); ind1 = intersect(ind,ind1);
            ind = find(~isnan(Z(2,1,:)) & ~isnan(Zstd(2,1,:))); ind1 = intersect(ind,ind1);
            ind = find(~isnan(Z(2,2,:)) & ~isnan(Zstd(2,2,:))); ind1 = intersect(ind,ind1);
            periods = periods(ind1); Z = Z(:,:,ind1); Zstd = Zstd(:,:,ind1);
            
            % Compute the discrete time impulse response
            if errorbars == 0
                [Zt,A4real] = DTIR(Z,periods,delta_t,[nmin nmax],option,lambda);
            elseif errorbars == 1
                [Zt,A4real] = DTIR(Z,periods,delta_t,[nmin nmax],option,lambda,Zstd);
            end
            impulse.xx = squeeze(Zt(1,1,:));
            impulse.xy = squeeze(Zt(1,2,:));
            impulse.yx = squeeze(Zt(2,1,:));
            impulse.yy = squeeze(Zt(2,2,:));
            
            % Save some TF metadata
            impulse.tfname = obj.tfname;
            impulse.lat = obj.lat;
            impulse.lon = obj.lon;
            impulse.elev = obj.elev;
            impulse.orient = obj.orient + obj.theta0;
            
            % Save some TS metadata
            impulse.delta_t = delta_t;
            impulse.nbounds = nbounds;
            impulse.option = option;
            impulse.lambda = lambda;
            impulse.errorbars = errorbars;
            
        end  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [newobj] = interp(obj,periods,method,extrap)
            % interpolates all values, including errors, to the new
            % periods; set periods outside of the domain to NaNs unless
            % required otherwise by using extrap = 1
            
            if nargin < 2
                error('Please specify periods in seconds');
            end
            
            if nargin < 3
                method = 'linear';
            end
            
            if nargin < 4
                extrap = 0;
            end
             
            newobj = obj;
            newobj.T = periods;
            
            nche = obj.Nche;
            
            % dealing with repeated period values
            [Ti,ii] = unique(obj.T);

            newobj.TF = zeros(2,nche,length(periods)) + 1i*zeros(2,nche,length(periods));
            for j = 1:nche
                newobj.TF(1,j,:) = interp1((log10(Ti)'),squeeze(obj.TF(1,j,ii)).*sqrt(Ti)',log10(periods),method,'extrap');
                newobj.TF(1,j,:) = squeeze(newobj.TF(1,j,:))./sqrt(periods)';
                newobj.TF(2,j,:) = interp1((log10(Ti)'),squeeze(obj.TF(2,j,ii)).*sqrt(Ti)',log10(periods),method,'extrap');
                newobj.TF(2,j,:) = squeeze(newobj.TF(2,j,:))./sqrt(periods)';
            end
            
            if ~obj.predicted
                newobj.TFVar = zeros(2,nche,length(periods));
                for j = 1:nche
                    newobj.TFVar(1,j,:) = interp1((log10(Ti)'),squeeze(obj.TFVar(1,j,ii)).*sqrt(Ti)',log10(periods),method,'extrap');
                    newobj.TFVar(1,j,:) = squeeze(newobj.TFVar(1,j,:))./sqrt(periods)';
                    newobj.TFVar(2,j,:) = interp1((log10(Ti)'),squeeze(obj.TFVar(2,j,ii)).*sqrt(Ti)',log10(periods),method,'extrap');
                    newobj.TFVar(2,j,:) = squeeze(newobj.TFVar(2,j,:))./sqrt(periods)';
                end
            end
            
            if ~isempty(obj.SIG_S)
                newobj.SIG_S = zeros(2,2,length(periods)) + 1i*zeros(2,2,length(periods));
                for j = 1:2
                    for i = 1:2
                        newobj.SIG_S(i,j,:) = interp1((log10(Ti)'),squeeze(obj.SIG_S(i,j,ii)).*sqrt(Ti)',log10(periods),method,'extrap');
                        newobj.SIG_S(i,j,:) = squeeze(newobj.SIG_S(i,j,:))./sqrt(periods)';
                    end
                end
            end
            
            if ~isempty(obj.SIG_E)
                newobj.SIG_E = zeros(nche,nche,length(periods)) + 1i*zeros(nche,nche,length(periods));
                for j = 1:nche
                    for i = 1:nche
                        newobj.SIG_E(i,j,:) = interp1((log10(Ti)'),squeeze(obj.SIG_E(i,j,ii)).*sqrt(Ti)',log10(periods),method,'extrap');
                        newobj.SIG_E(i,j,:) = squeeze(newobj.SIG_E(i,j,:))./sqrt(periods)';
                    end
                end
            end

            periodMin = min(obj.T);
            periodMax = max(obj.T);
           
            % dealing with periods outside of the defined range
            iper = find(periods < periodMin | periods > periodMax);
            if ~extrap
                newobj.TF(:,:,iper) = NaN + 1i*NaN;
                if ~isempty(newobj.TFVar)
                    newobj.TFVar(:,:,iper) = NaN;
                end
                if ~isempty(newobj.SIG_S)
                    newobj.SIG_S(:,:,iper) = NaN + 1i*NaN;
                end
                if ~isempty(newobj.SIG_E)
                    newobj.SIG_E(:,:,iper) = NaN + 1i*NaN;
                end
            end
            
            % dealing with additional Z-file info for output
            newobj.ndf = zeros(1,length(periods));
            newobj.ibandlim = zeros(2,length(periods));
            newobj.level = zeros(1,length(periods));
            newobj.sampRate = zeros(1,length(periods));
            newobj.nbt = length(periods);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = rotate(obj,newtheta,neworient)
            % [obj] = rotate(obj,newtheta,neworient)
            %
            % rotates any transfer function, orthogonal or otherwise,
            % to orthogonal coordinates with specified azimuth newtheta
            % relative to geographic North. Defaults to orthogonal
            % geographic. 
            %
            % Advanced functionality: if neworient is specified, FIRST
            % rotates to orthogonal geographic coordinates, THEN rotates
            % to channel orientations defined by newtheta + neworient.
            % Thus, this really supports rotation between any two
            % coordinate systems.
            %
            % Uses the algorithm from EMTF-FCU, (a) A. Kelbert, 2017
            
            fwdORinv = 'FWD';
            if nargin < 2
                disp(['Rotating ',obj.tfname,' to orthogonal geographic...']);
                newtheta = 0.0;
            elseif length(newtheta) == 1
                disp(['Rotating ',obj.tfname,' to ',num2str(newtheta),' angle to geographic North.']);
            else
                error('New angle to geographic North should be a scalar');
            end
            
            % update the rotation metadata to new orthogonal coordinates
            if nargin < 3
                for i=1:length(obj.chid)
                    if ~isempty(strfind(lower(obj.chid{i}),'x'))
                        neworient(i) = 0.0;
                    elseif strfind(lower(obj.chid{i}),'y')
                        neworient(i) = 90.0;
                    elseif ~isempty(strfind(lower(obj.chid{i}),'z'))
                        neworient(i) = 0.0;
                    else
                        error('Unknown channel ID - cannot determine orientation');
                    end
                end
                obj.orthogonal = 1;
            else
                % instead of checking for a bunch of things here, first
                % rotate to orthogonal geographic!!! - from there, setup
                % an "inverse" rotation to arbitrary orientations.
                if length(neworient) ~= obj.Nch
                    error('Orientations aren''t specified for all channels, cannot rotate');
                end
                obj = obj.rotate(0.0);
                fwdORinv = 'INV';
                obj.orthogonal = 0;
            end         
            
            if strcmpi(fwdORinv,'FWD')
                [U,V] = mttf.setup_rotation(obj.orient(1,:),obj.theta0,newtheta);
            else
                [U,V] = mttf.setup_rotation(neworient,newtheta,obj.theta0,'INV');
            end
            
            if isempty(obj.SIG_S) || isempty(obj.SIG_E)
                warning('Attempting a rotation without the full covariance matrix.')
                warning('The rotated error bars cannot be robustly interpreted.')
                fullcov = 0;
            else
                fullcov = 1;
            end
            
            for k=1:length(obj.T)
                
                data = squeeze(obj.TF(:,:,k)).';
                data = V * data * U';
                obj.TF(:,:,k) = data.';
                
                if ~fullcov
                    % In the absence of covariance matrices, rotate the variances matrix just like we rotate the TF
                    % and let's be clear that this is not the right way to do it!
                    var = squeeze(obj.TFVar(:,:,k)).';
                    var = V * var * U';
                    obj.TFVar(:,:,k) = var.';
                    
                else
                    S = squeeze(obj.SIG_S(:,:,k));
                    S = U * S * U';
                    obj.SIG_S(:,:,k) = S;
                    N = squeeze(obj.SIG_E(:,:,k));
                    N = V * N * V';
                    obj.SIG_E(:,:,k) = N;
                    
                    % finally, update the variances
                    % A.K. NOTE AS OF 8/14/2017: the variance of the real or imaginary component
                    % would be this divided by 2. However, p 52 of the EDI manual explains that e.g. ZXY.VAR
                    % refers to the complex variance (while ZXYR.VAR would mean the variance of the real part).
                    % To be consistent with this definition, we are no longer dividing by two.
                    % This was fixed in the reading routine, but not here, so are now fixing the XML files
                    % in the database.
                    for i=1:obj.Nche
                        for j=1:2
                            var(i,j) = (N(i,i)*S(j,j));
                        end
                    end
                    obj.TFVar(:,:,k) = var.';
                    
                end
            end
            
            % update the rotation metadata to new coordinates
            obj.orient = neworient;
            obj.theta0 = newtheta;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = truncate(obj,periodRange)
            % truncates to the period range [periodMin periodMax]
            
            if length(periodRange) ~= 2
                error('Please specify [periodMin periodMax] in seconds');
            elseif periodRange(2) <= periodRange(1)
                error('Maximum period needs to be greater than minimum period for truncation');
            end
            
            periodMin = periodRange(1);
            periodMax = periodRange(2);
           
            iper = find(obj.T >= periodMin & obj.T <= periodMax);
            obj.TF = obj.TF(:,:,iper);
            if ~isempty(obj.TFVar)
                obj.TFVar = obj.TFVar(:,:,iper);
            end
            if ~isempty(obj.SIG_S)
                obj.SIG_S = obj.SIG_S(:,:,iper);
                obj.SIG_E = obj.SIG_E(:,:,iper);
            end
            obj.T = obj.T(iper);
            if ~isempty(obj.ndf)
                obj.ndf = obj.ndf(iper);
            end
            if ~isempty(obj.ibandlim)
                obj.ibandlim = obj.ibandlim(:,iper);
            end
            if ~isempty(obj.level)
                obj.level = obj.level(iper);
            end
            if ~isempty(obj.sampRate)
                obj.sampRate = obj.sampRate(iper);
            end
            obj.nbt = length(iper);
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = insert(obj,otherobj,param)
            %
            % [obj] = insert(obj,otherobj,param)
            %
            % insert certain period or a certain mode from another mttf object 
            % into this object to replace the bad ones - will only work if this 
            % is the same site but different period bands or modes
            %
            % options are:
            % 1) insert a period range: e.g. 
            %    param = [3000 10000]
            % will replace these periods with those from otherobj mttf
            % 2) insert a channel: e.g.
            %    param = 'EY'
            % will replace all modes that use this channel with otherobj
            % currently mutually exclusive but that can be changed.
            %
            % use this carefully! we're not requiring that the TF name is
            % the same because it might be different coming from different
            % files... will check the latitude and longitude but that's
            % error prone so make sure you know what you're going when you
            % use this.
            
            if abs(obj.lat - otherobj.lat) > 1e-6 || abs(obj.lon - otherobj.lon) > 1e-6
                warning(['Merging MT TFs for sites ' ...
                    strtrim(obj.siteChar) ' and ' strtrim(otherobj.siteChar) ...
                    ' with different locations ' ...
                    '... make sure that these are the same site!']);
            end
            
            % Make sure that channel orientations match, too. On rare
            % occasions, might want to merge even when some channel
            % orientations don't match - will deal with these if they
            % arise.
            if max(abs((obj.theta0+obj.orient) - (otherobj.theta0+otherobj.orient))) > 1e-6
                % rotate before merging.
            end
            
            % By default, insert the full TF
            if nargin <= 2
                param = [min(otherobj.T) max(otherobj.T)];
            end
                
            if isnumeric(param) % INSERT A PERIOD RANGE
                if length(param) ~= 2
                    error('To merge by period, please specify [periodMin periodMax] in seconds');
                elseif param(2) <= param(1)
                    error('Maximum period needs to be greater than minimum period in a period range');
                end
                if param(1) > max(obj.T)
                    disp('Period ranges in MT TFs do not intersect; merging...');
                    obj1 = obj;
                    obj2 = otherobj.truncate(param);
                elseif param(2) < min(obj.T)
                    disp('Period ranges in MT TFs do not intersect; merging...');
                    obj1 = otherobj.truncate(param);
                    obj2 = obj;
                elseif min(obj.T) < param(1) && max(obj.T) > param(2)
                    disp(['Replacing MT TF periods longer than '...
                        num2str(param(1)) ' secs and shorter than '...
                        num2str(param(2)) ' secs with data from ' strtrim(otherobj.siteChar)]);
                    obj1 = obj.truncate([min(obj.T) param(1)]);
                    obj2 = otherobj.truncate([param(1) param(2)]);
                    obj3 = obj.truncate([param(2) max(obj.T)]);
                elseif min(obj.T) <= param(1) && max(obj.T) <= param(2)
                    disp('Period ranges in MT TFs partially intersect; merging...');
                    obj1 = obj.truncate([min(obj.T) param(1)]);
                    obj2 = otherobj.truncate([param(1) param(2)]);
                elseif min(obj.T) >= param(1) && max(obj.T) >= param(2)
                    disp('Period ranges in MT TFs partially intersect; merging...');
                    obj1 = otherobj.truncate([param(1) param(2)]);
                    obj2 = obj.truncate([param(2) max(obj.T)]);
                elseif max(obj.T) > param(1) 
                    disp(['Replacing MT TF periods longer than '...
                        num2str(param(1)) ' secs with data from ' strtrim(otherobj.siteChar)]);
                    obj1 = obj.truncate([min(obj.T) param(1)]);
                    obj2 = otherobj.truncate(param);
                elseif min(obj.T) < param(2)
                    disp(['Replacing MT TF periods shorter than '...
                        num2str(param(2)) ' secs with data from ' strtrim(otherobj.siteChar)]);
                    obj1 = otherobj.truncate(param);
                    obj2 = obj.truncate([param(2) max(obj.T)]);
                else
                    error('Can''t merge by period: period ranges not compatible');
                end
                obj = mttf.merge(obj1,obj2);
                if exist('obj3','var')
                    obj = mttf.merge(obj,obj3);
                end
            end
            
            if ischar(param) % REPLACE OR INSERT AN OUTPUT CHANNEL
                [ch1,chtype1] = obj.findChannel(param);
                [ch2,chtype2] = otherobj.findChannel(param);
                if isempty(ch1) && ~isempty(ch2)
                    disp(['Channel ' param ' is not found in MT TF ' obj.siteChar '. Inserting at the end.']);
                    ch1 = obj.Nch+1;
                elseif isempty(ch2)
                    error(['Channel ' param ' is not found in MT TF ' otherobj.siteChar]);
                elseif strcmp(chtype1,'input') || strcmp(chtype2,'input')
                    error('Can only merge output channels.');
                end
                otherobj = otherobj.truncate([min(obj.T) max(obj.T)]);
                obj = obj.truncate([min(otherobj.T) max(otherobj.T)]);
                obj.tfname = [obj.tfname '_' otherobj.tfname param];
                obj.TF(:,ch1-2,:) = otherobj.TF(:,ch2-2,:);
                if ~isempty(otherobj.TFVar)
                    obj.TFVar(:,ch1-2,:) = otherobj.TFVar(:,ch2-2,:);
                end
                if ~isempty(otherobj.SIG_E)
                    obj.SIG_E(ch1-2,ch1-2,:) = otherobj.SIG_E(ch2-2,ch2-2,:);
                end
                obj.chid{ch1} = otherobj.chid{ch2};
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function d = norm(obj,dim)
            % d = norm(obj,dim)
            % 
            % Computes the Frobenois norm of the impedance only
            % as a frequency-dependent array.
            % If dimensionality dim=1,2,3 extracts the relevant component
            % and outputs the norm of that.
            
            if nargin < 2
                Z = obj.impedance;
            elseif dim==1
                Z = obj.extract1d;
            elseif dim==2
                [ZA,ZD] = obj.extract2d;
                Z = ZA + ZD;
            elseif dim==3
                Z = obj.extract3d;
            else % default
                Z = obj.impedance;
            end
            
            d = zeros(1,length(obj.T));
            for k = 1:length(obj.T)
                d(k) = norm(squeeze(Z(:,:,k)),'fro');
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Z1,p1,obj1d] = extract1d(obj)
            % [z1,p1,obj1d] = extract1d(obj)
            % 
            % Get the closest approximation to 1D impedance using a crude
            % rule: set diagonal elements to zero and average the
            % amplitudes of the off-diagonals.
            %
            % To make things formal, use Jeff Love's formulation here:
            % 
            % Z = z1*[ 0 1] + z2A*[0 1] + z2D*[1  0] + z3*[1 0]
            %        [-1 0]       [1 0]       [0 -1]      [0 1]
            %
            % Doing nothing to the vertical field transfer functions.
            %
            % note that in Gary's TF matrix, index 1 is input, index 2 is
            % output - meaning that e.g. (ignoring Hz and period)
            % TF(1,2) is Ey/Hx = Zyx NOT Zxy as one might assume.
            % here we transpose to make it simple to use impedance matrix.
            % Trivial to extend the code for general number of output channels.
            % If more than two channels, assume that the first is Hz, the
            % others are the electrics. Can add additional checks.

            Z = impedance(obj);
            z1 = (Z(1,2,:) - Z(2,1,:))/2;
            
            Z1 = Z;
            Z1(1,1,:) = 0;
            Z1(1,2,:) = z1;
            Z1(2,1,:) = - z1;
            Z1(2,2,:) = 0;
            
            p1 = zeros(1,size(Z,3));
            for iper = 1:size(Z,3)
                Z1Dp = squeeze(Z1(:,:,iper));
                Zp = squeeze(Z(:,:,iper));
                p1(iper) = norm(Z1Dp,'fro')/norm(Zp,'fro');
            end

            if nargout>2
                if obj.Nche > 2
                    k = 1;
                else
                    k = 0;
                end
                obj1d = obj;
                for i = 1:2
                    for j = 1:2
                        obj1d.TF(j,k+i,:) = Z1(i,j,:);
                    end
                end
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Z2A,Z2D,p2A,p2D,obj2d] = extract2d(obj)
            % [z2A,z2D,p2A,p2D,obj2d] = extract2d(obj)
            % 
            % Really just getting the parts of obj that are not
            % rotationally invariant.
            %
            % To make things formal, use Jeff Love's formulation here:
            % 
            % Z = z1*[ 0 1] + z2A*[0 1] + z2D*[1  0] + z3*[1 0]
            %        [-1 0]       [1 0]       [0 -1]      [0 1]
            %
            % Doing nothing to the vertical field transfer functions.

            Z = impedance(obj);
            z2A = (Z(1,2,:) + Z(2,1,:))/2;
            z2D = (Z(1,1,:) - Z(2,2,:))/2;

            Z2A = Z;
            Z2A(1,1,:) = 0;
            Z2A(1,2,:) = z2A;
            Z2A(2,1,:) = z2A;
            Z2A(2,2,:) = 0;

            Z2D = Z;
            Z2D(1,1,:) = z2D;
            Z2D(1,2,:) = 0;
            Z2D(2,1,:) = 0;
            Z2D(2,2,:) = - z2D;
            
            p2A = zeros(1,size(Z,3));
            for iper = 1:size(Z,3)
                Z2Ap = squeeze(Z2A(:,:,iper));
                Zp = squeeze(Z(:,:,iper));
                p2A(iper) = norm(Z2Ap,'fro')/norm(Zp,'fro');
            end

            p2D = zeros(1,size(Z,3));
            for iper = 1:size(Z,3)
                Z2Dp = squeeze(Z2D(:,:,iper));
                Zp = squeeze(Z(:,:,iper));
                p2D(iper) = norm(Z2Dp,'fro')/norm(Zp,'fro');
            end

%             p2 = zeros(1,size(Z,3));
%             for iper = 1:size(Z,3)
%                 Z2p = squeeze(Z2A(:,:,iper))+squeeze(Z2D(:,:,iper));
%                 Zp = squeeze(Z(:,:,iper));
%                 p2(iper) = norm(Z2p,'fro')/norm(Zp,'fro');
%             end
            
            if nargout>4
                if obj.Nche > 2
                    k = 1;
                else
                    k = 0;
                end
                obj2d = obj;
                for i = 1:2
                    for j = 1:2
                        obj2d.TF(j,k+i,:) = Z2A(i,j,:) + Z2D(i,j,:);
                    end
                end
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Z3,p3,obj3d] = extract3d(obj)
            % [z3,p3,obj3d] = extract3d(obj)
            % 
            % Getting the parts of obj that are not 1D but
            % rotationally invariant.
            %
            % To make things formal, use Jeff Love's formulation here:
            % 
            % Z = z1*[ 0 1] + z2A*[0 1] + z2D*[1  0] + z3*[1 0]
            %        [-1 0]       [1 0]       [0 -1]      [0 1]
            %
            % Doing nothing to the vertical field transfer functions.

            Z = impedance(obj);
            z3 = (Z(1,1,:) + Z(2,2,:))/2;

            Z3 = Z;
            Z3(1,1,:) = z3;
            Z3(1,2,:) = 0;
            Z3(2,1,:) = 0;
            Z3(2,2,:) = z3;
            
            p3 = zeros(1,size(Z,3));
            for iper = 1:size(Z,3)
                Z3p = squeeze(Z3(:,:,iper));
                Zp = squeeze(Z(:,:,iper));
                p3(iper) = norm(Z3p,'fro')/norm(Zp,'fro');
            end
            
            if nargout>2
                if obj.Nche > 2
                    k = 1;
                else
                    k = 0;
                end
                obj3d = obj;
                for i = 1:2
                    for j = 1:2
                        obj3d.TF(j,k+i,:) = Z3(i,j,:);
                    end
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Z,Zstd,Sigma,SigmaDiag] = impedance(obj)
            %Z     means exactly what it should: complex impedance tensor
            %Zstd  standard deviation for impedance
            %Sigma we can make a complete covariance matrix if you want it
            %SigmaDiag usual diagonal covariance matrix that scales by error
            %
            % assume that HZ is the first output channel followed by
            % electrics - this function ignores HZ
            %
            % FULL AND DIAGONAL COVARIANCE NOT TESTED YET
            
            % note that in Gary's TF matrix, index 1 is input, index 2 is
            % output - meaning that e.g. (ignoring Hz and period) 
            % TF(1,2) is Ey/Hx = Zyx NOT Zxy as one might assume.
            % here we transpose to make a simple to use impedance matrix.
            % Trivial to extend the code for general number of output channels.
            % If more than two channels, usually the first is Hz, the
            % others are the electrics. But this should handle other
            % configurations. Can easily generalize for off diagonals.

            jj = findNche(obj,'Full_Impedance');

            if isempty(jj)
                disp(['No full impedance available for ' obj.tfname]);
                Z = [];
                Zstd = [];
                return
            elseif length(jj) ~= 2
                error('Something went wrong while searching for full impedance indices');
            end
            
            Z(1,1,:) = obj.TF(1,jj(1),:); % Zxx
            Z(1,2,:) = obj.TF(2,jj(1),:); % Zxy
            Z(2,1,:) = obj.TF(1,jj(2),:); % Zyx
            Z(2,2,:) = obj.TF(2,jj(2),:); % Zyy
            
            if isempty(obj.TFVar)
                for i=1:2
                    for j=jj
                        obj.TFVar(i,j,:) = real(obj.SIG_S(i,i,:).*obj.SIG_E(j,j,:));
                    end
                end
            end                
            
            % Zstd is standard deviation of the real/imag component
            Zstd(1,1,:) = sqrt(obj.TFVar(1,jj(1),:)/2); % Zxx
            Zstd(1,2,:) = sqrt(obj.TFVar(2,jj(1),:)/2); % Zxy
            Zstd(2,1,:) = sqrt(obj.TFVar(1,jj(2),:)/2); % Zyx
            Zstd(2,2,:) = sqrt(obj.TFVar(2,jj(2),:)/2); % Zyy
                       
            if nargout > 2
                nper = length(obj.T);
                Sigma = zeros(4*nper,4*nper);
                for j = 1:nper
                    Sigma(1*j,1*j) = N(1,1,j)*S(1,1,j);
                    Sigma(1*j,2*j) = N(1,1,j)*S(1,2,j);
                    Sigma(1*j,3*j) = N(1,1,j)*S(2,1,j);
                    Sigma(1*j,4*j) = N(1,1,j)*S(2,2,j);
                    Sigma(2*j,1*j) = N(1,2,j)*S(1,1,j);
                    Sigma(2*j,2*j) = N(1,2,j)*S(1,2,j);
                    Sigma(2*j,3*j) = N(1,2,j)*S(2,1,j);
                    Sigma(2*j,4*j) = N(1,2,j)*S(2,2,j);
                    Sigma(3*j,1*j) = N(2,1,j)*S(1,1,j);
                    Sigma(3*j,2*j) = N(2,1,j)*S(1,2,j);
                    Sigma(3*j,3*j) = N(2,1,j)*S(2,1,j);
                    Sigma(3*j,4*j) = N(2,1,j)*S(2,2,j);
                    Sigma(4*j,1*j) = N(2,2,j)*S(1,1,j);
                    Sigma(4*j,2*j) = N(2,2,j)*S(1,2,j);
                    Sigma(4*j,3*j) = N(2,2,j)*S(2,1,j);
                    Sigma(4*j,4*j) = N(2,2,j)*S(2,2,j);
                end
            end
            
            if nargout > 3
                SigmaDiag = zeros(4*nper,4*nper);
                for j = 1:(4*nper)
                    SigmaDiag(j,j) = 1/(abs(Sigma(j,j))^2);
                end
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [T,Tstd] = verticalFieldTFs(obj)
            %T     vertical magnetic field TFs
            %Tstd  standard deviation for vertical magnetic field TFs
            %
            % N = obj.SIG_E(1,1,:);
            % S = obj.SIG_S;
            % Tstd(1,:) = sqrt(real(N(1,1,:).*S(1,1,:))/2);
            % Tstd(2,:) = sqrt(real(N(1,1,:).*S(2,2,:))/2);
            %
            % Should do more checking based on obj.types strings
            % but this is not a priority right now
            %
            % FULL AND DIAGONAL COVARIANCE NOT IMPLEMENTED YET
            
            jj = findNche(obj,'Full_Vertical_Components');

            if isempty(jj)
                disp(['No vertical field TFs available for ' obj.tfname]);
                T = [];
                Tstd = [];
                return
            elseif length(jj) ~= 1
                error('Something went wrong while searching for tipper indices');                
            end
                
            T = obj.TF(:,jj,:); 
            
            if isempty(obj.TFVar)
                for i=1:2
                    for j=jj
                        obj.TFVar(i,j,:) = real(obj.SIG_S(i,i,:).*obj.SIG_E(j,j,:));
                    end
                end
            end                
            
            % Zstd is standard deviation of the real/imag component
            Tstd(1,:) = sqrt(obj.TFVar(1,jj,:)/2); % Tx
            Tstd(2,:) = sqrt(obj.TFVar(2,jj,:)/2); % Ty
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = struct2mttf(obj,Zstruct)
            %   create an mttf object from a Zstruct as used in
            %   readZfile.m and Zplot.m in EMTF/matlab
            
            % create a new mttf object only if it doesn't exist
            %if ~exist('obj','mttf')
            %    obj = mttf;
            %end
            
            % copy all structure fields to the mttf object
            for fn = fieldnames(Zstruct)'
                obj.(fn{1}) = Zstruct.(fn{1});
            end
            
            % compute the variance (divide by two for real/imag parts)
            for i=1:2
                for j=1:obj.Nche
                    obj.TFVar(i,j,:) = real(obj.SIG_S(i,i,:).*obj.SIG_E(j,j,:));
                end
            end
            
            % update all dependent fields
            obj.lat = obj.stcor(1);
            obj.lon = obj.stcor(2);
            obj.siteChar = obj.stname;
            obj.nComp = 2*obj.Nche;
            
            if obj.nComp == 12
                obj.compChar = ['TX ';'TY ';'ZXX';'ZXY';'ZYX';'ZYY'];
                obj.types = {'Full_Vertical_Components','Full_Impedance'};
            elseif obj.nComp == 8
                obj.compChar = ['ZXX';'ZXY';'ZYX';'ZYY'];
                obj.types = {'Full_Impedance'};
            end
            
            % always assume not orthogonal unless we know otherwise
            obj.orthogonal = 0;
                
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [index,type] = findChannel(obj,chname)
            %
            % Usage:  [index,type] = findChannel(obj,chname)
            
            ch = zeros(1,obj.Nch);
            for i = 1:obj.Nch
                if cell2mat(strfind(lower(cellstr(obj.chid(i))),lower(chname)));
                    ch(i) = 1;
                end
            end
            index = find(ch);
            
            type = cell(1,length(index));
            for i = 1:length(index)
                if index(i) <= 2
                    type{i} = 'input';
                else
                    type{i} = 'output';
                end
            end
            
        end
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [index] = findNche(obj,type)
            %
            % Usage:  [index] = findNche(obj,type)
            %
            % Utility function. Terse and perhaps too general but this
            % deals with multi-type objects and allows to locate the
            % indices of the output channels for any specific type.
            % Used in impedance and verticalFieldTFs functions.
            % Currently assumes two input channels. Can be changed.
            
            jmax = find(strncmp(type,obj.types,40));
            nche = zeros(1,jmax);
            for j = 1:jmax
                ncomp = mttf.components(obj.types{j});
                if obj.Cmplx
                    nche(j) = ncomp/2 / 2;
                else
                    nche(j) = ncomp / 2;
                end
            end
            index = sum(nche(1:jmax-1))+1:sum(nche(1:jmax));
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = mark(obj,mcolor,msize)
            
            % [] = mark(obj,mcolor,msize)
            %
            % Marks a cross at a specific site location on a pcolor map.
            % The site can be defined by [lat lon] or an mttf.
            %
            % Optionally, can specify the color of the cross (default
            % black).
                         
            if nargin < 1
                disp('Usage: mttf.mark([lat lon],[0.5 0.5 0.5])');
                return
            end
            
            if isobject(obj)
                xlat = obj.lat;
                xlon = obj.lon;
            elseif isnumeric(obj)
                xlat = obj(1);
                xlon = obj(2);
            else
                disp('Cannot figure out the lat and lon in mttf.mark');
                return
            end
            
            if nargin < 2
                mcolor = [0.1 0.1 0.1];
            end
            
            if nargin < 3
                msize = 14;
            end
            
            m_plot(xlon,xlat,'x','color',mcolor,'MarkerSize',msize,'LineWidth',2)

        end
                 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f]= Zplot(obj)
            % cleaned up version of EMTF/matlab/ZPLT/Zplot.m
            % which should also work with these mttf objects
            
            [Z,Zstd] = impedance(obj);
            
            f = figure('Position',[100,100,950,750],...
                'PaperPosition',[1,1,9.5,7.5]);
            if obj.predicted
                sreal = '-';
                simag = '-';
            else
                sreal = '+';
                simag = '*';
            end
            xl = .40;
            yl = .40;
            sp = .05;
            y0 = .5;
            YL = 50;
            XL = [1,100000];
            IMPlabels = {'Zxx','Zxy','Zyx','Zyy'};
            jk = 0;
            for j = 1:2
                y = y0 - (j-1)*(yl+sp);
                x0 = .1;
                for k = 1:2
                    jk = jk+1;
                    x = x0 + (k-1)*(xl+sp);
                    axes('Position',[x,y,xl,yl])
                    ZT = squeeze(Z(j,k,:)).*sqrt(obj.T)';
                    ZTstd = squeeze(Zstd(j,k,:)).*sqrt(obj.T)';
                    errorbar(obj.T,real(ZT),ZTstd,['b' sreal],'markersize',8,'linewidth',2); 
                    hold on
                    errorbar(obj.T,imag(ZT),ZTstd,['r' simag],'markersize',8,'linewidth',2);
                    semilogx(XL,[0,0],'k--')
                    set(gca,'xscale','log');
                    set(gca,'FontWeight','bold','FontSize',14)
                    xlim(XL)
                    ylim([-YL,YL])
                    text(3000,YL*.8,IMPlabels{jk},'FontWeight','bold',...
                        'FontSize',18)
                end
            end
            legend('Real','Imag','Location','SouthWest')
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f,p] = apresplt(obj,xlims,yfactor,modes)
            % [f,p] = apresplt(obj,xlims,yfactor,modes)
            %
            % Usage options: obj.apresplt([1 1e5],100,'xy,yx,xx,yy')
            % Usage options: obj.apresplt([1 1e5],[1e-3 1e1],'xy,yx,xx,yy')
            %                obj.apresplt(gcf)
            %
            % plot the apparent resistivity and phase from an impedance
            % object; designed to imitate EMTF/matlab/ZPLT/pltrhom.m
            % but written from scratch. This one is cleaner, more readable
            % and has additional functionality but no menus
            %
            % can print using
            %   print(f,'-dpng',[obj.tfname '_apres']);
            %
            % if needed, individual axes can be obtained using
            % rho_axes = f.Children(3);
            % rho_label = f.Children(2);
            % phs_axes = f.Children(1);
            
            % First check for alternative usage to plot on top of existing
            % figure
            existing_figure = 0;
            fixed_range = 0;
            if nargin > 1
                if ~isnumeric(xlims)
                    existing_figure = 1;
                    f = xlims;
                end
            end

            if existing_figure
                ax = get(f,'CurrentAxes');
                xlims = ax.XLim;
            else
                if nargin < 2
                    xlims(1) = min(min(obj.T),1);
                    xlims(2) = max(max(obj.T),1e5);
                elseif isempty(xlims) || ~isnumeric(xlims)
                    clear xlims;
                    xlims(1) = min(min(obj.T),1);
                    xlims(2) = max(max(obj.T),1e5);
                end
                xlims(1) = 10^(log10(xlims(1)));
                xlims(2) = 10^(log10(xlims(2)));
            end
            
            if existing_figure
                ax = get(f,'CurrentAxes');
                yfactor = ax.YLim(2)/ax.YLim(1);
            else
                if nargin < 3
                    yfactor = 100;
                elseif length(yfactor) == 2
                    ymin = yfactor(1);
                    ymax = yfactor(2);
                    yfactor = sqrt(ymax/ymin);
                    fixed_range = 1;
                end
            end

            if nargin < 4
                modes = 'xy,yx';
            end
            modes = upper(modes);

            
            if obj.predicted
                sxy = '-';
                syx = '-';
                sxx = '-';
                syy = '-';
            else
                sxy = 'o';
                syx = 'x';
                sxx = '*';
                syy = '^';
            end

            [Z,Zstd] = obj.impedance;
            [apres,phase] = obj.imp2apres(obj.T,Z,Zstd);
            
            % warn about electrode polarity errors
            pol = mttf.test_electrode_polarity(phase);
            if (pol(1) && pol(2))
                warning('Ex and Ey polarity reversed');
            elseif (pol(1))
                warning('Ex polarity reversed');
            elseif (pol(2))
                warning('Ey polarity reversed');
            end
            
            % convert the phases to 1st quadrant for plotting
            phase1 = phase;
            phase1.xy = atand(tand(phase.xy));
            phase1.yx = atand(tand(phase.yx));
            phase1.xx = atand(tand(phase.xx));
            phase1.yy = atand(tand(phase.yy));
            
            % colors extracted from EMTF/matlab/ZPLT/set_fig.m for
            % compatibility with Gary Egbert's plots
            symbol_colors = [ 0 0 .7; .7 0 0; 0 .7 0; .35 0 .35; 0 .35 .35 ; 0 0 0 ; .3 .3 0];

            siteName = obj.tfname;
            strlon = num2str(obj.lon);
            strlat = num2str(obj.lat);
            per = obj.T;
            
            ii = per>xlims(1) & per<xlims(2);
            per = per(ii);
            apres.xx = apres.xx(ii);
            apres.xy = apres.xy(ii);
            apres.yx = apres.yx(ii);
            apres.yy = apres.yy(ii);
            if ~obj.predicted
                apres.xx_se = apres.xx_se(ii);
                apres.xy_se = apres.xy_se(ii);
                apres.yx_se = apres.yx_se(ii);
                apres.yy_se = apres.yy_se(ii);
            end
            phase1.xx = phase1.xx(ii);
            phase1.xy = phase1.xy(ii);
            phase1.yx = phase1.yx(ii);
            phase1.yy = phase1.yy(ii);
            if ~obj.predicted
                phase1.xx_se = phase1.xx_se(ii);
                phase1.xy_se = phase1.xy_se(ii);
                phase1.yx_se = phase1.yx_se(ii);
                phase1.yy_se = phase1.yy_se(ii);
            end
                
            
            if ~existing_figure
                f=figure('Position',[200 200 750 900],'PaperPosition',[0.25 0.5 7.5 9]); clf;
            else
                figure(f); hold on;
            end
            
            % Plot apparent resistivities
            subplot(10,1,1:7); 
            i = 1;
            if strfind(modes,'XY')
                p(i)=loglog(per,apres.xy(:),sxy,'color',symbol_colors(1,:),'linewidth',1.5); hold on
                if ~obj.predicted
                    p(i)=errorbar(per,apres.xy(:),2*apres.xy_se(:),sxy,'color',symbol_colors(1,:),'linewidth',1.5);
                end
                str(i,:) = 'XY'; i = i+1; hold on; 
            end
            if strfind(modes,'YX')
                p(i)=loglog(per,apres.yx(:),syx,'color',symbol_colors(2,:),'linewidth',1.5); hold on
                if ~obj.predicted
                    p(i)=errorbar(per,apres.yx(:),2*apres.yx_se(:),syx,'color',symbol_colors(2,:),'linewidth',1.5);
                end
                str(i,:) = 'YX'; i = i+1; hold on; 
            end
            if strfind(modes,'XX')
                p(i)=loglog(per,apres.xx(:),sxx,'color',symbol_colors(3,:),'linewidth',1.5); hold on
                if ~obj.predicted
                    p(i)=errorbar(per,apres.xx(:),2*apres.xx_se(:),sxx,'color',symbol_colors(3,:),'linewidth',1.5);
                end
                str(i,:) = 'XX'; i = i+1; hold on; 
            end
            if strfind(modes,'YY')
                p(i)=loglog(per,apres.yy(:),syy,'color',symbol_colors(4,:),'linewidth',1.5); hold on
                if ~obj.predicted
                    p(i)=errorbar(per,apres.yy(:),2*apres.yy_se(:),syy,'color',symbol_colors(4,:),'linewidth',1.5);
                end
                str(i,:) = 'YY'; i = i+1; hold on; 
            end
            hold on

            jmin = find(per>xlims(1),1,'first');
            jmax = find(per<xlims(2),1,'last');
            ymean = (nanmedian(apres.xy(jmin:jmax))+nanmedian(apres.yx(jmin:jmax)))/2;
            if isnan(ymean)>0 || ymean==0, ymean=0.01;
               fprintf('BAD PNG: ');
             end
            if ~existing_figure
                if fixed_range
                    set(gca,'xlim',xlims,'ylim',[ymin ymax]);
                else
                    set(gca,'xlim',xlims,'ylim',[ymean/yfactor ymean*yfactor]);
                end
                set(gca,'xminortick','on','yminortick','on','xgrid','on','ygrid','on');
            end
            if ~existing_figure
                xlim1 = get(gca,'xlim'); ylim1 = get(gca,'ylim');
                pbaspect([log10(xlim1(2)/xlim1(1)) log10(ylim1(2)/ylim1(1)) 1]); % axes square
                posn1 = get(gca,'position');
            end
            set(gca,'fontsize',13,'fontweight','bold');
            legend(p,str(1:i-1,:),'location','southeast');
            ylabel('\rho_a (\Omega m)','fontsize',16,'fontweight','demi');
            if ~existing_figure
                %RMS = sqrt(nansum(obj.v.res(i,:))/sum(~isnan(obj.v.err(i,:))));
                %title([char(OBS) ' [ LON = ' lon '; LAT = ' lat '; RMS = ' num2str(RMS)  ' ]'],'fontweight','demi');
                title([char(siteName) ' [ LON = ' strlon '; LAT = ' strlat ' ]'],...
                    'fontsize',16,'fontweight','bold','interpreter','none');
            end
            
            % Plot phases
            subplot(10,1,8:10); 
            if strfind(modes,'XY')
                p(i)=semilogx(per,phase1.xy(:),sxy,'color',symbol_colors(1,:),'linewidth',1.5); hold on
                if ~obj.predicted
                    p(i)=errorbar(per,phase1.xy(:),2*phase1.xy_se(:),sxy,'color',symbol_colors(1,:),'linewidth',1.5);
                end
                i = i+1; hold on
            end
            if strfind(modes,'YX')
                p(i)=semilogx(per,phase1.yx(:),syx,'color',symbol_colors(2,:),'linewidth',1.5); hold on
                if ~obj.predicted
                    p(i)=errorbar(per,phase1.yx(:),2*phase1.yx_se(:),syx,'color',symbol_colors(2,:),'linewidth',1.5);
                end
                i = i+1; hold on;
            end
            if strfind(modes,'XX')
                inphase = find(phase1.xx>0);
                ibreak = find(diff(inphase)>1);
                if isempty(ibreak); ibreak = 0; end
                if ibreak > 0
                    p(i)=semilogx(per(inphase(1:ibreak)),phase1.xx(inphase(1:ibreak)),sxx,'color',symbol_colors(3,:),'linewidth',1.5);
                    i = i+1; hold on;
                end
                p(i)=semilogx(per(inphase(ibreak+1:end)),phase1.xx(inphase(ibreak+1:end)),sxx,'color',symbol_colors(3,:),'linewidth',1.5);
                i = i+1; hold on;
                if ~obj.predicted
                    p(i)=errorbar(per,phase1.xx(:),2*phase1.xx_se(:),sxx,'color',symbol_colors(3,:),'linewidth',1.5);
                end
                hold on
            end
            if strfind(modes,'YY')
                inphase = find(phase1.yy>0);
                ibreak = find(diff(inphase)>1);
                if isempty(ibreak); ibreak = 0; end
                if ibreak > 0
                    p(i)=semilogx(per(inphase(1:ibreak)),phase1.yy(inphase(1:ibreak)),syy,'color',symbol_colors(4,:),'linewidth',1.5);
                    i = i+1; hold on;
                end
                p(i)=semilogx(per(inphase(ibreak+1:end)),phase1.yy(inphase(ibreak+1:end)),syy,'color',symbol_colors(4,:),'linewidth',1.5);
                i = i+1; hold on;
                if ~obj.predicted
                    p(i)=errorbar(per,phase1.yy(:),2*phase1.yy_se(:),syy,'color',symbol_colors(4,:),'linewidth',1.5);
                end
                hold on
            end
            if ~existing_figure
                posn2 = get(gca,'position');
                set(gca,'position',[posn1(1) posn2(2) posn1(3) posn2(4)]);
            end
            hold off
            %ymean = (nanmean(phase1.xy(:))+nanmean(phase1.yx(:)))/2;
            %set(gca,'xlim',xlims,'ylim',[ymean-45 ymean+45]);
            set(gca,'xlim',xlims,'ylim',[0 90]);
            set(gca,'xminortick','on','yminortick','on','xgrid','on','ygrid','on');
            set(gca,'fontsize',13,'fontweight','bold');
            ylabel('\phi (Degrees)','fontsize',16,'fontweight','demi');
            xlabel('Period (secs)','fontsize',16,'fontweight','demi');

            %print(f,'-dpng',[siteName '_apres']);
                        
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f,p] = impplt(obj,dataType,f)
            % Sample usage: [f,p] = obj.impplt(1,gcf)
            % Use this to plot the impedances, or the predicted impedances
            % on top of measured data (or vice versa). Can also use this
            % for vertical field transfer functions.
            % 
            % dataType = 1 for Full_Impedance
            %          = 2 for Full_Vertical_Components
            % by default, plot the impedances
            %
            % can print using
            %   print(f,'-dpng',[obj.tfname '_impedance']);
            %   print(f,'-dpng',[obj.tfname '_verticalTFs']);
                         
            if nargin < 2
                dataType = 1;
            end
            
            existing_figure = 0;
            if nargin >= 3
                existing_figure = 1;
            end
            
            % colors extracted from EMTF/matlab/ZPLT/set_fig.m for
            % compatibility with Gary Egbert's plots and apresplt
            symbol_colors = [ 0 0 .7; .7 0 0; 0 .7 0; .35 0 .35; 0 .35 .35 ; 0 0 0 ; .3 .3 0];
            
            % old colors: [0. 0.5 0.8], [0.8 0.5 0]

            fileName = [obj.tfdir '/' obj.tfname];
            siteName = obj.tfname;
            strlon = num2str(obj.lon);
            strlat = num2str(obj.lat);
            per = obj.T';
            range = 60; %100; %20;
            xlims(1) = 10^floor(log10(min(per)));
            xlims(2) = 10^ceil(log10(max(per)));
            
            if obj.predicted
                sreal = 'x-';
                simag = 'o-';
                symbol_colors = [ 0 0 0; 0 0 0; 0 0 0; 0 0 0];
                W = 2;
            else
                sreal = 'x';
                simag = 'o';
                W = 1;
            end
            
            if ~isempty(find(dataType==1, 1))
                [Z,Zstd] = obj.impedance;
                if isempty(Z)
                    disp(['Unable to plot the impedances for ' siteName ': not found']);
                else
                    if ~existing_figure
                        f=figure('Position',[200,200,1400,900],...
                            'PaperPosition',[1,1,18,10],...
                            'PaperOrientation','Portrait'); clf;
                    else
                        figure(f); hold on;
                    end
                    % Diagonal TFs
                    subplot(2,2,1)
                    data = squeeze(Z(1,1,:)).*sqrt(per);
                    std2 = 2*squeeze(Zstd(1,1,:)).*sqrt(per);
                    TFinfo = [char(siteName) ' (' strlon '; ' strlat ')'];
                    p(1)=semilogx(per,real(data),sreal,'color',symbol_colors(3,:),'linewidth',1.5*W); hold on
                    if ~ obj.predicted
                        p(1)=errorbar(per,real(data),std2,sreal,'color',symbol_colors(3,:),'linewidth',1.5*W); hold on
                    end
                    p(2)=semilogx(per,imag(data),simag,'color',symbol_colors(3,:),'linewidth',W); hold on
                    if ~ obj.predicted
                        p(2)=errorbar(per,imag(data),std2,simag,'color',symbol_colors(3,:),'linewidth',W); hold on
                    end
                    ymean = real(nanmedian(data)); if isnan(ymean); ymean = 0; end
                    if ~existing_figure
                        if range<=1e-6; range = 60;end
                        set(gca,'xlim',xlims,'ylim',[ymean-range/2 ymean+range/2]);
                        set(gca,'xminortick','on','yminortick','on');
                        set(gca,'fontsize',13,'fontweight','bold');
                        legend(p(1:2),'Re(ZXX)','Im(ZXX)');
                        semilogx(xlims,[0,0],'k--'); hold on
                        ylabel('ZXX * sqrt(T)','fontsize',16,'fontweight','demi');
                        xlabel('Period (secs)','fontsize',16,'fontweight','demi');
                        title(TFinfo,'fontsize',16,'fontweight','demi','interpreter','none');
                    end
                    subplot(2,2,4)
                    data = squeeze(Z(2,2,:)).*sqrt(per);
                    std2 = 2*squeeze(Zstd(2,2,:)).*sqrt(per);
                    TFinfo = [siteName ' (' strlon '; ' strlat ')'];
                    p(3)=semilogx(per,real(data),sreal,'color',symbol_colors(4,:),'linewidth',1.5*W); hold on
                    if ~ obj.predicted
                        p(3)=errorbar(per,real(data),std2,sreal,'color',symbol_colors(4,:),'linewidth',1.5*W); hold on
                    end
                    p(4)=semilogx(per,imag(data),simag,'color',symbol_colors(4,:),'linewidth',W); hold on
                    if ~ obj.predicted
                        p(4)=errorbar(per,imag(data),std2,simag,'color',symbol_colors(4,:),'linewidth',W); hold on
                    end
                    ymean = real(nanmedian(data)); if isnan(ymean); ymean = 0; end                   
                    if ~existing_figure
                        if range<=1e-6; range = 60;end
                        set(gca,'xlim',xlims,'ylim',[ymean-range/2 ymean+range/2]);
                        set(gca,'xminortick','on','yminortick','on');
                        set(gca,'fontsize',13,'fontweight','bold');
                        legend(p(3:4),'Re(ZYY)','Im(ZYY)');
                        semilogx(xlims,[0,0],'k--'); hold on
                        ylabel('ZYY * sqrt(T)','fontsize',16,'fontweight','demi');
                        xlabel('Period (secs)','fontsize',16,'fontweight','demi');
                        title(TFinfo,'fontsize',16,'fontweight','demi','interpreter','none');
                    end
                    % Off-diagonal TFs
                    subplot(2,2,2)
                    data = squeeze(Z(1,2,:)).*sqrt(per);
                    std2 = 2*squeeze(Zstd(1,2,:)).*sqrt(per);
                    TFinfo = [siteName ' (' strlon '; ' strlat ')'];
                    p(5)=semilogx(per,real(data),sreal,'color',symbol_colors(1,:),'linewidth',1.5*W); hold on
                    if ~ obj.predicted
                        p(5)=errorbar(per,real(data),std2,sreal,'color',symbol_colors(1,:),'linewidth',1.5*W); hold on
                    end
                    p(6)=semilogx(per,imag(data),simag,'color',symbol_colors(1,:),'linewidth',W); hold on
                    if ~ obj.predicted
                        p(6)=errorbar(per,imag(data),std2,simag,'color',symbol_colors(1,:),'linewidth',W); hold on
                    end
                    ymean = real(nanmedian(data)); if isnan(ymean); ymean = 0; end
                    if ~existing_figure
                        if range<=1e-6; range = 60;end
                        set(gca,'xlim',xlims,'ylim',[ymean-range/2 ymean+range/2]);
                        set(gca,'xminortick','on','yminortick','on');
                        set(gca,'fontsize',13,'fontweight','bold');
                        legend(p(5:6),'Re(ZXY)','Im(ZXY)');
                        semilogx(xlims,[0,0],'k--'); hold on
                        ylabel('ZXY * sqrt(T)','fontsize',16,'fontweight','demi');
                        xlabel('Period (secs)','fontsize',16,'fontweight','demi');
                        title(TFinfo,'fontsize',16,'fontweight','demi','interpreter','none');
                    end
                    subplot(2,2,3)
                    data = squeeze(Z(2,1,:)).*sqrt(per);
                    std2 = 2*squeeze(Zstd(2,1,:)).*sqrt(per);
                    TFinfo = [siteName ' (' strlon '; ' strlat ')'];
                    p(7)=semilogx(per,real(data),sreal,'color',symbol_colors(2,:),'linewidth',1.5*W); hold on
                    if ~ obj.predicted
                        p(7)=errorbar(per,real(data),std2,sreal,'color',symbol_colors(2,:),'linewidth',1.5*W); hold on
                    end
                    p(8)=semilogx(per,imag(data),simag,'color',symbol_colors(2,:),'linewidth',W); hold on
                    if ~ obj.predicted
                        p(8)=errorbar(per,imag(data),std2,simag,'color',symbol_colors(2,:),'linewidth',W); hold on
                    end
                    ymean = real(nanmedian(data)); if isnan(ymean); ymean = 0; end
                    if ~existing_figure
                        if range<=1e-6; range = 60;end
                        set(gca,'xlim',xlims,'ylim',[ymean-range/2 ymean+range/2]);
                        set(gca,'xminortick','on','yminortick','on');
                        set(gca,'fontsize',13,'fontweight','bold');
                        legend(p(7:8),'Re(ZYX)','Im(ZYX)');
                        semilogx(xlims,[0,0],'k--'); hold off
                        ylabel('ZYX * sqrt(T)','fontsize',16,'fontweight','demi');
                        xlabel('Period (secs)','fontsize',16,'fontweight','demi');
                        title(TFinfo,'fontsize',16,'fontweight','demi','interpreter','none');
                    end
                    %print(f,'-dpng',[fileName '_impedance']);
                end
            elseif ~isempty(find(dataType==2, 1))
                % Vertical TFs
                [HT,HTstd] = obj.verticalFieldTFs;
                if isempty(HT)
                    disp(['Unable to plot vertical field TFs for ' siteName ': not found']);
                else
                    if nargin < 3
                        f=figure('Position',[200,200,900,900],...
                            'PaperPosition',[1,1,10,10],...
                            'PaperOrientation','Portrait'); clf;
                    else
                        figure(f); hold on;
                    end
                    subplot(2,1,1)
                    data = HT(1,:);
                    std2 = 2*HTstd(1,:);
                    TFinfo = [siteName ' (' strlon '; ' strlat ')'];
                    p(1)=semilogx(per,real(data),['m' sreal],'linewidth',1.5); hold on
                    if ~ obj.predicted
                        p(1)=errorbar(per,real(data),std2,['m' sreal],'linewidth',1.5); hold on
                    end
                    p(2)=semilogx(per,imag(data),['g' simag],'linewidth',1.5); hold on
                    if ~ obj.predicted
                        p(2)=errorbar(per,imag(data),std2,['g' simag],'linewidth',1.5); hold on
                    end
                    ymean = 0; %real(nanmedian(data));
                    if ~existing_figure
                        set(gca,'xlim',xlims,'ylim',[ymean-0.5 ymean+0.5]);
                        set(gca,'xminortick','on','yminortick','on');
                        set(gca,'fontsize',13,'fontweight','bold');
                        legend(p(1:2),'Re(TX)','Im(TX)');
                        semilogx(xlims,[0,0],'k--'); hold on
                        ylabel('TX','fontsize',16,'fontweight','demi');
                        xlabel('Period (secs)','fontsize',16,'fontweight','demi');
                        title(TFinfo,'fontsize',16,'fontweight','demi','interpreter','none');
                    end
                    subplot(2,1,2)
                    data = HT(2,:);
                    std2 = 2*HTstd(2,:);
                    %TFinfo = [siteName ' (' strlon '; ' strlat ')'];
                    p(3)=semilogx(per,real(data),['m' sreal],'linewidth',1.5); hold on
                    if ~ obj.predicted
                        p(3)=errorbar(per,real(data),std2,['m' sreal],'linewidth',1.5); hold on
                    end
                    p(4)=semilogx(per,imag(data),['g' simag],'linewidth',1.5); hold on
                    if ~ obj.predicted
                        p(4)=errorbar(per,imag(data),std2,['g' simag],'linewidth',1.5); hold on
                    end
                    ymean = 0; %real(nanmedian(data));
                    if ~existing_figure
                        set(gca,'xlim',xlims,'ylim',[ymean-0.5 ymean+0.5]);
                        set(gca,'xminortick','on','yminortick','on');
                        set(gca,'fontsize',13,'fontweight','bold');
                        legend(p(3:4),'Re(TY)','Im(TY)');
                        semilogx(xlims,[0,0],'k--'); hold off
                        ylabel('TY','fontsize',16,'fontweight','demi');
                        xlabel('Period (secs)','fontsize',16,'fontweight','demi');
                        %title(TFinfo,'fontsize',16,'fontweight','demi','interpreter','none');
                    end
                    %print(f,'-dpng',[fileName '_verticalTFs']);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [info] = write_xml(obj,cfile,project,survey,year)
            % [info] = write_xml(obj,cfile,project,survey,year)
            % 
            % Write the minimal EMTF XML format data file. Uses xml_io_tools.
            % Only write the impedance Z and the vertical field TFs T
            % to the file, or whichever of the two are found in obj.
            %
            % If the file name is not specified, this puts together and outputs
            % the info structure for external modification. To write it out
            % later, use xml_write(cfile,info).
            %
            % Project, survey & year are required metadata. If no metadata 
            % information exists, these can optionally be included as 
            % input arguments. Project is an abbreviation. Survey can be
            % descriptive.
            %
            % Note: THIS IS THE MINIMAL FILE. IT IS NOT INTENDED TO BE
            % COMPREHENSIVE AND IS TO BE USED ONLY FOR CONVENIENCE, NOT FOR
            % ARCHIVING. SIGNAL COVARIANCE MATRICES ARE NOT WRITTEN.
            % COPYRIGHT INFORMATION IS NOT WRITTEN. FOR ANY LONG TERM
            % PURPOSES, the Fortran EMTF File Conversion Utilities (FCU)
            % ARE TO BE USED, INSTEAD.
            %
            % @ Anna Kelbert, USGS. May 5, 2016. Last mod. July 6, 2017.
            
            [Z,Zstd] = impedance(obj);
            [T,Tstd] = verticalFieldTFs(obj);
            
            if ~isempty(Z) && ~isempty(T)
                tags = 'impedance,tipper';
            elseif ~isempty(Z)
                tags = 'impedance';
            elseif ~isempty(T)
                tags = 'tipper';
            end
                                               
            if nargin < 5
                if isfield(obj.metadata,'YEAR')
                    year = obj.metadata.YEAR;
                else
                    today = datevec(now); year = today(1);
                end
            end
            if nargin < 4
                if isfield(obj.metadata,'SURVEY')
                    survey = obj.metadata.SURVEY;
                else
                    survey = '';
                end
            end
            if nargin < 3
                if isfield(obj.metadata,'PROJECT')
                    project = obj.metadata.PROJECT;
                else
                    project = '';
                end
            end
            if isnumeric(year); year = num2str(year); end
            
            % ProductId is critical: it needs to be unique in the SPUD EMTF
            % database for possible future archiving, so take care with it
            ProductId = [strtrim(project) '.' ...
                        strtrim(obj.tfname) '.' ...
                        year];
            disp(ProductId);
                                    
            info.Description = 'Magnetotelluric Transfer Functions';
            info.ProductId = ProductId;
            info.SubType = 'MT_TF';
            info.Notes = ['Written from Matlab. Original format: ' obj.primaryFormat];
            info.Tags = tags;
            info.PrimaryData = [];
            info.Attachment = [];
            info.Provenance.CreateTime = datestr(now,'yyyy-mm-ddTHH:MM:SS');
            info.Provenance.CreatingApplication = 'mttf.write_xml in Matlab; author Anna Kelbert';
            info.Copyright.ReleaseStatus = 'Unrestricted Release';
            info.Site.Id = strtrim(obj.tfname);
            if ~isempty(obj.siteChar)
                info.Site.Name = strtrim(obj.siteChar);
            else
                info.Site.Name = '';
            end
            info.Site.Project = project;
            info.Site.Survey = survey;
            info.Site.YearCollected = year;
            info.Site.Location.Latitude = obj.lat;
            info.Site.Location.Longitude = obj.lon;
            if isempty(obj.elev)
                info.Site.Location.Elevation.CONTENT = 0;
            else
                info.Site.Location.Elevation.CONTENT = obj.elev;
            end
            info.Site.Location.Elevation.ATTRIBUTE.units = 'meters';
            if obj.orthogonal
                info.Site.Orientation.CONTENT = 'orthogonal';
                info.Site.Orientation.ATTRIBUTE.angle_to_geographic_north = obj.theta0 + obj.orient(1);
            else
                info.Site.Orientation = 'sitelayout';
            end
            info.Site.Start = [year(1:4) '-01-01T00:00:00'];
            info.Site.End = [year(end-3:end) '-12-31T24:00:00'];
            if obj.signConvention > 0
                info.ProcessingInfo.SignConvention = 'exp(+ i\omega mu)';
            else
                info.ProcessingInfo.SignConvention = 'exp(- i\omega mu)';
            end
            info.ProcessingInfo.ProcessingTag = info.Site.Id;

            if ~isempty(obj.metadata)
                str = fields(obj.metadata);
                for i = 1:length(str)
                    switch str{i}
                        case 'IMAGE'
                            info.PrimaryData.Filename = obj.metadata.(str{i});
                        case 'ORIGINAL'
                            original = obj.metadata.(str{i});
                            if ischar(original)
                                info.Attachment.Filename = original;
                            else
                                for k=1:length(original)
                                    info.Attachment{k}.Filename = original{k};
                                end
                            end
                        case 'ORTHOGONAL'
                            % overwritten by obj.orthogonal
                        %    orthogonal = obj.metadata.(str{i});
                        %    if ~ orthogonal
                        %        info.Site.Orientation.CONTENT = 'SiteLayout';
                        %        clear info.Site.Orientation.ATTRIBUTE;
                        %    end
                        case 'SITEID'
                            info.Site.Id = obj.metadata.(str{i});
                        case 'SITENAME'
                            info.Site.Name = obj.metadata.(str{i});
                        case 'PROJECT'
                            % already there
                        case 'SURVEY'
                            % already there
                        case 'YEAR'
                            % already there
                        case 'COUNTRY'
                            info.Site.Country = obj.metadata.(str{i});
                        case 'ORIENT'
                            % overwritten by obj.theta0 + obj.orient
                        case 'SIGNCONVENTION'
                            % overwritten by obj.signConvention
                        case 'DECLINATION'
                            info.Site.Declination.CONTENT = obj.metadata.(str{i});
                            info.Site.Declination.ATTRIBUTE.epoch = '1995.0';
                        case 'START'
                            info.Site.Start = obj.metadata.(str{i});
                        case 'END'
                            info.Site.End = obj.metadata.(str{i});
                        case 'ACQUIREDBY'
                            info.Site.AcquiredBy = obj.metadata.(str{i});
                        case 'PROCESSEDBY'
                            info.ProcessingInfo.ProcessedBy = obj.metadata.(str{i});
                        case 'PROCESSINGSOFTWARE'
                            info.ProcessingInfo.ProcessingSoftware.Name = obj.metadata.(str{i});
                        case 'PROCESSINGSOFTWARELASTMOD'
                            info.ProcessingInfo.ProcessingSoftware.LastMod = obj.metadata.(str{i});
                        case 'PROCESSINGSOFTWAREAUTHOR'
                            info.ProcessingInfo.ProcessingSoftware.Author = obj.metadata.(str{i});
                        case 'PROCESSINGTAG'
                            info.ProcessingInfo.ProcessingTag = obj.metadata.(str{i});
                        case 'RUNLIST'
                            % won't bother recording these separately
                            info.ProcessingInfo.FieldNotes.RunList = obj.metadata.(str{i});
                        case 'REMOTEREF'
                            info.ProcessingInfo.RemoteRef.ATTRIBUTE.type = obj.metadata.(str{i});
                        case 'REMOTESITE'
                            info.ProcessingInfo.RemoteInfo.Site.Id = obj.metadata.(str{i});
                        case 'QUALITYRATING'
                            info.Site.DataQualityNotes.Rating = obj.metadata.(str{i});
                        case 'GOODFROMPERIOD'
                            info.Site.DataQualityNotes.GoodFromPeriod = obj.metadata.(str{i});
                        case 'GOODTOPERIOD'
                            info.Site.DataQualityNotes.GoodToPeriod = obj.metadata.(str{i});
                        case 'QUALITYCOMMENTS'
                            info.Site.DataQualityNotes.Comments.ATTRIBUTE.author = obj.metadata.PROCESSEDBY;
                            info.Site.DataQualityNotes.Comments.CONTENT = obj.metadata.(str{i});
                        case 'WARNINGFLAG'
                            info.Site.DataQualityWarnings.Flag = obj.metadata.(str{i});
                        case 'WARNINGCOMMENTS'
                            info.Site.DataQualityNotes.Comments.ATTRIBUTE.author = obj.metadata.PROCESSEDBY;
                            info.Site.DataQualityWarnings.Comments.CONTENT = obj.metadata.(str{i});
                        case 'UNITS'
                            % overwritten by obj.units
                        otherwise
                            % ignore for now
                            disp(['Metadata field ' str{i} ' not recorded in xml']);
                    end
                end
            end
            
            if ~isempty(obj.provenance)
                info.Provenance = obj.provenance;
            end
            
            if ~isempty(obj.copyright)
                info.Copyright = obj.copyright;
            end
            
            if isempty(obj.orient) || isempty(obj.theta0)
                theta = [0 90 0 0 90];
            elseif length(obj.orient) == 1
                theta = obj.theta0 + obj.orient + [0 90 0 0 90];
            else
                theta = obj.theta0 + obj.orient;
            end  
            % writing the channels... this can be generalized to death, but
            % currently don't see the point of going through this trouble
            ind = obj.findChannel('Hx');
            if ~isempty(ind)
                info.SiteLayout.InputChannels.Magnetic(1).ATTRIBUTE.name = obj.chid{ind};
                info.SiteLayout.InputChannels.Magnetic(1).ATTRIBUTE.orientation = theta(ind);
            end
            ind = obj.findChannel('Hy');
            if ~isempty(ind)
                info.SiteLayout.InputChannels.Magnetic(2).ATTRIBUTE.name = obj.chid{ind};
                info.SiteLayout.InputChannels.Magnetic(2).ATTRIBUTE.orientation = theta(ind);
            end
            ind = obj.findChannel('Hz');
            if ~isempty(ind)
                info.SiteLayout.OutputChannels.Magnetic(1).ATTRIBUTE.name = obj.chid{ind};
                info.SiteLayout.OutputChannels.Magnetic(1).ATTRIBUTE.orientation = theta(ind);
            end
            ind = obj.findChannel('Ex');
            if ~isempty(ind)
                info.SiteLayout.OutputChannels.Electric(1).ATTRIBUTE.name = obj.chid{ind};
                info.SiteLayout.OutputChannels.Electric(1).ATTRIBUTE.orientation = theta(ind);
            end
            ind = obj.findChannel('Ey');
            if ~isempty(ind)
                info.SiteLayout.OutputChannels.Electric(2).ATTRIBUTE.name = obj.chid{ind};
                info.SiteLayout.OutputChannels.Electric(2).ATTRIBUTE.orientation = theta(ind);
            end
            
            datatype.impedance.Tag = 'impedance';
            datatype.impedance.Intention = 'primary data type';
            datatype.impedance.ExternalUrl = 'http://www.iris.edu/dms/products/emtf/impedance.html';
            datatype.impedance.Description = 'MT impedance';
            datatype.impedance.ATTRIBUTE.name = 'Z';
            datatype.impedance.ATTRIBUTE.type = 'complex';
            datatype.impedance.ATTRIBUTE.input = 'H';
            datatype.impedance.ATTRIBUTE.output = 'E';
            datatype.impedance.ATTRIBUTE.units = '[mV/km]/[nT]';
            
            datatype.tipper.Tag = 'tipper';
            datatype.tipper.Intention = 'primary data type';
            datatype.tipper.ExternalUrl = 'http://www.iris.edu/dms/products/emtf/tipper.html';
            datatype.tipper.Description = 'Vertical Field Transfer Functions (Tipper)';
            datatype.tipper.ATTRIBUTE.name = 'T';
            datatype.tipper.ATTRIBUTE.type = 'complex';
            datatype.tipper.ATTRIBUTE.input = 'H';
            datatype.tipper.ATTRIBUTE.output = 'H';
            datatype.tipper.ATTRIBUTE.units = '[]';
  
            if ~isempty(strfind(tags,'impedance')) && ~isempty(strfind(tags,'tipper'))
                info.DataTypes.DataType(1) = datatype.impedance;
                info.DataTypes.DataType(2) = datatype.tipper;
            elseif ~isempty(strfind(tags,'impedance'))
                info.DataTypes.DataType(1) = datatype.impedance;
            elseif ~isempty(strfind(tags,'tipper'))
                info.DataTypes.DataType(1) = datatype.tipper;
            end

            info.Data.ATTRIBUTE.count = length(obj.T);
            for j = 1:length(obj.T)
                info.Data.Period(j).ATTRIBUTE.units = 'secs';
                info.Data.Period(j).ATTRIBUTE.value = obj.T(j);
            end
            
            if strfind(tags,'impedance')
                type = 'Z';
                nchout = 2;
                size = [2 2];
                type_units = obj.units;
                element = [type '.VAR']; variance = strrep(element,'.','0x2E');
                for j = 1:length(obj.T)
                    % first write the data...
                    info.Data.Period(j).(type).ATTRIBUTE.size = size;
                    info.Data.Period(j).(type).ATTRIBUTE.type = 'complex';
                    info.Data.Period(j).(type).ATTRIBUTE.units = type_units;
                    info.Data.Period(j).(type).value(1).ATTRIBUTE.input = 'Hx';
                    info.Data.Period(j).(type).value(1).ATTRIBUTE.output = 'Ex';
                    info.Data.Period(j).(type).value(1).ATTRIBUTE.name = 'Zxx';
                    Zxx = [num2str(real(Z(1,1,j)),'%15.9e') ' ' num2str(imag(Z(1,1,j)),'%15.9e')];
                    info.Data.Period(j).(type).value(1).CONTENT = Zxx;
                    info.Data.Period(j).(type).value(2).ATTRIBUTE.input = 'Hy';
                    info.Data.Period(j).(type).value(2).ATTRIBUTE.output = 'Ex';
                    info.Data.Period(j).(type).value(2).ATTRIBUTE.name = 'Zxy';
                    Zxy = [num2str(real(Z(1,2,j)),'%15.9e') ' ' num2str(imag(Z(1,2,j)),'%15.9e')];
                    info.Data.Period(j).(type).value(2).CONTENT = Zxy;
                    info.Data.Period(j).(type).value(3).ATTRIBUTE.input = 'Hx';
                    info.Data.Period(j).(type).value(3).ATTRIBUTE.output = 'Ey';
                    info.Data.Period(j).(type).value(3).ATTRIBUTE.name = 'Zyx';
                    Zyx = [num2str(real(Z(2,1,j)),'%15.9e') ' ' num2str(imag(Z(2,1,j)),'%15.9e')];
                    info.Data.Period(j).(type).value(3).CONTENT = Zyx;
                    info.Data.Period(j).(type).value(4).ATTRIBUTE.input = 'Hy';
                    info.Data.Period(j).(type).value(4).ATTRIBUTE.output = 'Ey';
                    info.Data.Period(j).(type).value(4).ATTRIBUTE.name = 'Zyy';
                    Zyy = [num2str(real(Z(2,2,j)),'%15.9e') ' ' num2str(imag(Z(2,2,j)),'%15.9e')];
                    info.Data.Period(j).(type).value(4).CONTENT = Zyy;
                    % then write the variance.
                    if ~obj.predicted
                        info.Data.Period(j).(variance).ATTRIBUTE.size = size;
                        info.Data.Period(j).(variance).ATTRIBUTE.type = 'real';
                        info.Data.Period(j).(variance).value(1).ATTRIBUTE.input = 'Hx';
                        info.Data.Period(j).(variance).value(1).ATTRIBUTE.output = 'Ex';
                        info.Data.Period(j).(variance).value(1).ATTRIBUTE.name = 'Zxx';
                        Zxx = num2str(2*(Zstd(1,1,j)^2),'%15.9e');
                        info.Data.Period(j).(variance).value(1).CONTENT = Zxx;
                        info.Data.Period(j).(variance).value(2).ATTRIBUTE.input = 'Hy';
                        info.Data.Period(j).(variance).value(2).ATTRIBUTE.output = 'Ex';
                        info.Data.Period(j).(variance).value(2).ATTRIBUTE.name = 'Zxy';
                        Zxy = num2str(2*(Zstd(1,2,j)^2),'%15.9e');
                        info.Data.Period(j).(variance).value(2).CONTENT = Zxy;
                        info.Data.Period(j).(variance).value(3).ATTRIBUTE.input = 'Hx';
                        info.Data.Period(j).(variance).value(3).ATTRIBUTE.output = 'Ey';
                        info.Data.Period(j).(variance).value(3).ATTRIBUTE.name = 'Zyx';
                        Zyx = num2str(2*(Zstd(2,1,j)^2),'%15.9e');
                        info.Data.Period(j).(variance).value(3).CONTENT = Zyx;
                        info.Data.Period(j).(variance).value(4).ATTRIBUTE.input = 'Hy';
                        info.Data.Period(j).(variance).value(4).ATTRIBUTE.output = 'Ey';
                        info.Data.Period(j).(variance).value(4).ATTRIBUTE.name = 'Zyy';
                        Zyy = num2str(2*(Zstd(2,2,j)^2),'%15.9e');
                        info.Data.Period(j).(variance).value(4).CONTENT = Zyy;
                    end
                end
            end
            
            if strfind(tags,'tipper')
                type = 'T';
                nchout = 1;
                size = [1 2];
                type_units = '[]';
                element = [type '.VAR']; variance = strrep(element,'.','0x2E');
                for j = 1:length(obj.T)
                    % first write the data...
                    info.Data.Period(j).(type).ATTRIBUTE.size = size;
                    info.Data.Period(j).(type).ATTRIBUTE.type = 'complex';
                    info.Data.Period(j).(type).ATTRIBUTE.units = type_units;
                    info.Data.Period(j).(type).value(1).ATTRIBUTE.input = 'Hx';
                    info.Data.Period(j).(type).value(1).ATTRIBUTE.output = 'Hz';
                    info.Data.Period(j).(type).value(1).ATTRIBUTE.name = 'Tx';
                    Tx = [num2str(real(T(1,j)),'%15.9e') ' ' num2str(imag(T(1,j)),'%15.9e')];
                    info.Data.Period(j).(type).value(1).CONTENT = Tx;
                    info.Data.Period(j).(type).value(2).ATTRIBUTE.input = 'Hy';
                    info.Data.Period(j).(type).value(2).ATTRIBUTE.output = 'Hz';
                    info.Data.Period(j).(type).value(2).ATTRIBUTE.name = 'Ty';
                    Ty = [num2str(real(T(2,j)),'%15.9e') ' ' num2str(imag(T(2,j)),'%15.9e')];
                    info.Data.Period(j).(type).value(2).CONTENT = Ty;
                    % then write the variance.
                    if ~obj.predicted
                        info.Data.Period(j).(variance).ATTRIBUTE.size = size;
                        info.Data.Period(j).(variance).ATTRIBUTE.type = 'real';
                        info.Data.Period(j).(variance).value(1).ATTRIBUTE.input = 'Hx';
                        info.Data.Period(j).(variance).value(1).ATTRIBUTE.output = 'Hz';
                        info.Data.Period(j).(variance).value(1).ATTRIBUTE.name = 'Tx';
                        Tx = num2str(2*(Tstd(1,j)^2),'%15.9e');
                        info.Data.Period(j).(variance).value(1).CONTENT = Tx;
                        info.Data.Period(j).(variance).value(2).ATTRIBUTE.input = 'Hy';
                        info.Data.Period(j).(variance).value(2).ATTRIBUTE.output = 'Hz';
                        info.Data.Period(j).(variance).value(2).ATTRIBUTE.name = 'Ty';
                        Ty = num2str(2*(Tstd(2,j)^2),'%15.9e');
                        info.Data.Period(j).(variance).value(2).CONTENT = Ty;
                    end
                end
            end
            
            % Allows to write out the structure for external modification
            if nargin > 1
                Pref.StructItem = 0;
                Pref.CellItem = 0;
                xml_write(cfile,info,'EM_TF',Pref);
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = write(obj,cfile,format)
            %
            % Usage:  [] = write(obj,cfile,format)
            %
            %  Writes an mttf object to a Z-file, XML, EDI etc.
            
            if nargin < 2
                error(['Output file name ' cfile ' not specified']);
            elseif nargin < 3
                % determine file format from filename extension
                format = mttf.format(cfile);
            end
                        
            switch format
                case 'Z' % borrowed from EMTF/matlab/ZPLT/readZfile.m
                    if isempty(obj.sta)
                       for k=1:obj.Nch
                           obj.sta{k} = obj.tfname;
                       end
                    end
                    if size(obj.orient,1)==1
                        obj.orient(2,obj.Nch) = 0;
                    end
                    writeZfile(obj,cfile);
                case 'XML'
                    write_xml(obj,cfile);
                case 'EDI' % NOT IMPLEMENTED YET
                    error('EDI output is not implemented yet - please contribute that code in your spare time!')
                case 'BIRRP' % NOT IMPLEMENTED YET
                    error('BIRRP output is not implemented yet - please contribute that code in your spare time!')
                otherwise
                    error(['Can''t write MT impedance in the format ' format ': method unknown']);
            end
            
            disp(['MT transfer function ' strtrim(obj.siteChar) ' written to ' format '-file ' cfile]);
            
        end      
        
    end
        
    
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = read(cfile,format,type)
            %
            % Usage:  [obj] = read(cfile,format,type)
            %
            %  Reads an MT transfer function into mttf object
            %
            % Currently, not reading the units or the sign convention;
            % assumes practical units and exp(+ i\omega mu). This can and
            % should be generalized.
            
            obj = mttf;
            
            if nargin == 0
               [mtfile,mtdir] = uigetfile('*');
               cfile = [mtdir mtfile];
            end
            
            if nargin < 2
               format = mttf.format(cfile);
            end
            
            % For XML and EDI formats: loop over types by default
            % for backwards compatibility with old codes, some functionality 
            % currently requires that the tipper comes first, if present
            % so the order in supported_types is important:
            % supported_types = {'Full_Vertical_Components','Full_Impedance'};
            if nargin < 3
                if strcmpi(format,'XML') || strcmpi(format,'EDI')
                    obj = mttf();
                    supported_types = obj.types;
                    for k = 1:length(supported_types)
                        obj1 = mttf.read(cfile,format,supported_types{k});
                        if ~isempty(obj1.TF)
                            obj = obj.appendType(obj1);
                        end
                    end
                    return
                else
                    type = 'ALL';
                end
            end

            if ~exist(cfile,'file')
                error(['Data file ' cfile ' not found']);
            else
                disp(['Reading data type ',char(type),' from file ',cfile,'...']);
            end
                        
            switch upper(format)
                case 'Z' % borrowed from EMTF/matlab/ZPLT/readZfile.m
                    Z = readZfile(cfile);
                    obj = struct2mttf(obj,Z);
                case 'XML'
                    switch type
                        case {'Full_Impedance','Off_Diagonal_Impedance',...
                                'Off_Diagonal_Rho_Phase','Phase_Tensor'}
                            [periods,Z,Zstd,latlon,info,copyrightinfo] = mttf.read_xml(cfile,1);
                            obj = mttf('Full_Impedance');
                            obj = obj.set(periods,Z,Zstd,latlon,info);
                            obj.copyright = copyrightinfo;
                        case {'Full_Vertical_Components'}
                            [periods,Ti,Tstd,latlon,info,copyrightinfo] = mttf.read_xml(cfile,2);
                            obj = mttf('Full_Vertical_Components');
                            obj = obj.set(periods,Ti,Tstd,latlon,info);
                            obj.copyright = copyrightinfo;
                        case {'Full_Interstation_TF'}
                            disp(['Reading of ' type ' not coded yet']);
                        otherwise
                            disp(['Unknown data type ' type]);
                    end
%                    error('XML input is not implemented yet - feel free to adapt ModelPlot/matlab/dataTools/readXML.m')
%                     [obj.d,obj.units,obj.signConvention] ...
%                         = readXML(cfile,newUnits,dataType);
%                     % a patch to convert from structure to mtperiod class
%                     % to avoid modifying the data list reading code
%                     obj.nPeriods = length(obj.d);
%                     for i=1:obj.nPeriods
%                         obj.d{i} = mtperiod(dataType,obj.d{i});
%                     end
                case 'EDI'
                    switch type
                        case {'Full_Impedance','Off_Diagonal_Impedance',...
                                'Off_Diagonal_Rho_Phase','Phase_Tensor'}
                            [periods,Z,Zstd,latlon,info] = mttf.read_edi(cfile,1);
                            obj = mttf('Full_Impedance');
                            obj = obj.set(periods,Z,Zstd,latlon,info);
                        case {'Full_Vertical_Components'}
                            [periods,Ti,Tstd,latlon,info] = mttf.read_edi(cfile,2);
                            obj = mttf('Full_Vertical_Components');
                            obj = obj.set(periods,Ti,Tstd,latlon,info);
                        case {'Full_Interstation_TF'}
                            disp(['Reading of ' type ' not coded yet']);
                        otherwise
                            disp(['Unknown data type ' type]);
                    end
                case 'BIRRP'
                    [periods,apres,phase] = mttf.read_birrp(cfile);
                    [Z,Zstd] = mttf.apres2imp(periods,apres,phase);
                    obj = mttf('Full_Impedance');
                    obj = obj.set(periods,Z,Zstd);
                otherwise
                    error(['Can''t read MT impedance in the format ' format ': method unknown']);
            end
            
            % Now extract all that we can from the file name...
            ipt = strfind(cfile,'.');
            isl = strfind(cfile,'/'); if isempty(isl); isl = 0; end
            if isempty(ipt)
                error(['MT Impedance file name ',cfile,' needs to have a valid extension']);
            end
            % Homemade recipe doesn't work with SPUD bundles which include
            % dots in their directory names... use fileparts.
            %obj.tfdir = '.'; if isl>0; obj.tfdir = cfile(1:isl(end)); end
            %obj.tfname = cfile(isl(end)+1:ipt(1)-1);
            %obj.tfext = cfile(ipt(1)+1:end);
            [obj.tfdir,obj.tfname,obj.tfext] = fileparts(cfile);
            obj.primaryFormat = format;
            
        end      
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function format = format(filename)
            
            ipt = strfind(filename,'.');
            if isempty(ipt)
                error(['MT Impedance file name ',filename,' needs to have a valid extension']);
            end
            
            ext = filename(ipt(end)+1:end);           
            switch ext
                case {'zss','zrr','zmm'}
                    format = 'Z'; % Gary Egbert
                case 'xml'
                    format = 'XML'; % Anna Kelbert
                case 'edi'
                    format = 'EDI'; % SEG 1987
                case 'rp'
                    format = 'BIRRP'; % Alan Chave
                otherwise
                    error(['Can''t determine the format of a file with extension ' ext ': method unknown']);
            end
            
            disp(['MT transfer function file format recognized as: ',format]);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Zrms,Trms] = rms(obj1,obj2)
            %
            % [Zrms,Trms] = rms(obj1,obj2)
            %
            % Computes the component by component Root Mean Square misfit
            % between two mttf objects. Assumes that one of the mttfs has
            % errorbar information (ie, ~obj.predicted). If both have
            % errorbar information, will use the first one for errorbars.
            %
            % Assumes that the periods are compatible. If not compatible,
            % run obj2 = obj2.interp(obj1.T) before calling the rms function.
            
            if length(obj1.T) ~= length(obj2.T)
                error(['Cannot compute the RMS for two mttfs with different ' ...
                    ' period sets. Try ''interp''.']);
            end
            
            periods = obj1.T;
            [Z1,Z1std] = impedance(obj1);
            [Z2,Z2std] = impedance(obj2);
            if ~obj1.predicted
                Zstd = Z1std;
            else
                Zstd = Z2std;
            end
            
            Zrms = zeros(2,2);
            for i = 1:2
                for j = 1:2
                    Zdiff = squeeze(Z1(i,j,:)-Z2(i,j,:));
                    Zrms(i,j) = nansum(squeeze(real(Zdiff)).^2./squeeze(Zstd(i,j,:)).^2) + ...
                        nansum(squeeze(imag(Zdiff)).^2./squeeze(Zstd(i,j,:)).^2);
                    ind = ~isnan(Zdiff) & ~isnan(squeeze(Zstd(i,j,:))); 
                    nper = length(periods(ind));
                    Zrms(i,j) = sqrt(Zrms(i,j)/(2*nper));
                end
            end
            
            if nargout > 1
                [T1,T1std] = verticalFieldTFs(obj1);
                [T2,T2std] = verticalFieldTFs(obj2);
                if ~obj1.predicted
                    Tstd = T1std;
                else
                    Tstd = T2std;
                end
                
                Trms = zeros(1,2);
                for i = 1:2
                    for j = 1:2
                        Tdiff = squeeze(T1(j,:)-T2(j,:));
                        ind = ~isnan(Tdiff) & ~isnan(squeeze(Tstd(j,:))); 
                        nper = length(periods(ind));
                        Trms(j) = nansum(squeeze(real(Tdiff)).^2./squeeze(Tstd(j,:)).^2) + ...
                            nansum(squeeze(imag(Tdiff)).^2./squeeze(Tstd(j,:)).^2);
                        Trms(j) = sqrt(Trms(j))/(2*nper);
                    end
                end
                
            end
            
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = merge(obj1,obj2)
            %
            % [obj] = merge(obj1,obj2)
            %
            % Will only merge two MT TFs if their period ranges do not
            % intersect; obj2 assumed a natural continuation of obj1.
            %
            % Related functions: rotate, truncate, insert.
            
            if max(obj1.T) > min(obj2.T)
                error(['Merging MT TFs ' strtrim(obj1.siteChar) ' and ' strtrim(obj2.siteChar) ...
                    ' by period not possible: period ranges cannot intersect. Try ''insert''.']);
            end
            
            % jump through some hoops to avoid duplicate periods
            delta = 1e-6 * min(obj1.T);
            if abs(min(obj2.T) - max(obj1.T)) < delta
                obj2 = obj2.truncate([obj2.T(1)+delta obj2.T(end)+delta]);
            end
            
            nper1 = length(obj1.T);
            nper2 = length(obj2.T);
            obj = obj1;
            obj.tfname = [obj1.tfname '_' obj2.tfname];
            obj.TF = zeros(2,obj1.Nche,nper1+nper2);
            obj.TF(:,:,1:nper1) = obj1.TF;
            obj.TF(:,:,nper1+1:nper1+nper2) = obj2.TF;
            if ~isempty(obj1.TFVar) && ~isempty(obj2.TFVar)
                obj.TFVar = zeros(2,obj1.Nche,nper1+nper2);
                obj.TFVar(:,:,1:nper1) = obj1.TFVar;
                obj.TFVar(:,:,nper1+1:nper1+nper2) = obj2.TFVar;
            end
            if ~isempty(obj1.SIG_S) && ~isempty(obj2.SIG_S)
                obj.SIG_S = zeros(2,2,nper1+nper2);
                obj.SIG_S(:,:,1:nper1) = obj1.SIG_S;
                obj.SIG_S(:,:,nper1+1:nper1+nper2) = obj2.SIG_S;
            end
            if ~isempty(obj1.SIG_E) && ~isempty(obj2.SIG_E)
                obj.SIG_E = zeros(obj1.Nche,obj1.Nche,nper1+nper2);
                obj.SIG_E(:,:,1:nper1) = obj1.SIG_E;
                obj.SIG_E(:,:,nper1+1:nper1+nper2) = obj2.SIG_E;
            end
            obj.T(1:nper1) = obj1.T;
            obj.T(nper1+1:nper1+nper2) = obj2.T;
            obj.ndf(1:nper1) = obj1.ndf;
            obj.ndf(nper1+1:nper1+nper2) = obj2.ndf;
            obj.ibandlim = zeros(2,nper1+nper2);
            obj.ibandlim(:,1:nper1) = obj1.ibandlim;
            obj.ibandlim(:,nper1+1:nper1+nper2) = obj2.ibandlim;
            obj.level = [obj1.level obj2.level];
            obj.sampRate = [obj1.sampRate obj2.sampRate];
            obj.nbt = nper1+nper2;
            disp(['New period range: ' num2str(min(obj.T)) ' to ' num2str(max(obj.T))]);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Tpred,Zpred,A4] = predictZ(impulse,delta_t,nbounds,Trange)
            
            % [Tpred,Zpred,A4] = predictZ(impulse,delta_t,nbounds,Trange)
            %
            % Complements the non-causal Discrete Time Impulse Response (DTIR) function.
            % Uses the impulse response output; recomputes the representer matrix A4
            % for the periods supplied as input and predicts the complex impedance
            % at these periods.
            %
            % For simplicity of computation, using the COMPLEX version of the
            % representer matrix A here. The REAL version was used in DTIR by
            % necessity. Both are equivalent.
            %
            % Inputs:
            %  impulse - discrete time inpulse response
            %  Trange - original periods in seconds (1xnper)
            %  delta_t - a scalar sampling rate of the time series in seconds
            %  nbounds - defines the length of the impulse response in samples
            %
            % Outputs:
            %  Tpred - periods in seconds (1xnper)
            %  Zpred - predicted complex impedance tensor (2x2xnper)
            %  A4 - 4-component complex representer matrix; Zpred = A4*zn
            %
            % (c) Anna Kelbert, 2015-2016, USGS Geomagnetism Program
            
            if nargin < 4
                Trange = [1 1000000];
            elseif nargin < 3
                error('Need all three input parameters to predict impedance');
            end
            
            Tpred = logspace(log10(Trange(1)*1e-1),log10(Trange(end)*100),300);
            
            % setup
            nmin = nbounds(1);
            nmax = nbounds(2);
            nout = (nmax-nmin+1);
            nper = length(Tpred);
            
            % check for input consistency
            % if length(impulse.xx) ~= nout
            %     disp('Impulse response and time lag sample bounds inconsistent');
            % end
            if length(impulse.xx) ~= nout
                disp('Impulse response and time lag sample bounds inconsistent');
            end
            
            % compute angular frequencies
            omega = zeros(1,nper);
            for j = 1:nper
                omega(j) = 2*pi/Tpred(j);
            end
            
            % set up the COMPLEX representer matrix A
            A = zeros(nper,nout) + 1i*zeros(nper,nout);
            for j = 1:nper % for all frequencies...
                for k = 1:(nmax-nmin+1) % and for n = ..., -1, 0, 1, ...
                    n = nmin+k-1;
                    A(j,k) = exp(-1i*omega(j)*n*delta_t);
                end
            end
            A4 = blkdiag(A,A,A,A);
            
            % compute Zpred = A4*zn
            zn = [impulse.xx; impulse.xy; impulse.yx; impulse.yy];
            %zn = reshape(impulse,1,4*nout); zn = zn.';
            tmp = A4*zn;
            
            % reshape for output
            %Zpred = reshape(tmp,nper,4); Zpred = Zpred.';
            % Zpred = reshape(Zpred,2,2,nper);
            Zpred(1,1,:) = tmp(1:nper); % xx
            Zpred(1,2,:) = tmp(nper+1:2*nper); % xy
            Zpred(2,1,:) = tmp(2*nper+1:3*nper); % yx
            Zpred(2,2,:) = tmp(3*nper+1:4*nper); % yy
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Tout,Zout] = padZ(T,Z,method)
            
            % [Tout,Zout] = padZ(T,Z,method)
            %
            % Used to artificially extend the period range of an impedance.
            %
            % Inputs:
            %  T - periods in seconds (1xnper)
            %  Z - the complex impedance tensor in MEASUREMENT coordinates (2x2xnper)
            %  method - 'linear','nearest','next','previous','spline','pchip','cubic'
            %            (using interp1)
            %            OR 'zeros'
            
            if nargin < 3
                method = 'pchip';
            end
            
            %Tout = logspace(1,6,50); % pad short and long periods
            Tout = logspace(1,6,50); % pad long periods only
            %iper = intersect(find(T>7),find(T<2e4));
            iper = 1:length(T);
            
            Zxx = squeeze(Z(1,1,:));
            Zxy = squeeze(Z(1,2,:));
            Zyx = squeeze(Z(2,1,:));
            Zyy = squeeze(Z(2,2,:));
            
            figure;
            
            if strcmp(method,'zeros')
                f_xx = zeros(1,length(Tout)) + 1i*zeros(1,length(Tout));
            else
                %[p,s,mu] = polyfit((log10(T(iper))'),Zxx(iper).*sqrt(T(iper))',order);
                %f_xx = polyval(p,log10(Tout),[],mu);
                f_xx = interp1((log10(T(iper))'),Zxx(iper).*sqrt(T(iper))',log10(Tout),method,'extrap');
            end
            subplot(2,2,1);
            semilogx(T,real(Zxx).*sqrt(T)','b*'); hold on;
            semilogx(Tout,real(f_xx),'g'); hold on;
            semilogx(T,imag(Zxx).*sqrt(T)','r*'); hold on;
            semilogx(Tout,imag(f_xx),'k'); hold on;
            ylabel('Zxx * sqrt(T)')
            
            if strcmp(method,'zeros')
                f_xy = zeros(1,length(Tout)) + 1i*zeros(1,length(Tout));
            else
                %[p,s,mu] = polyfit((log10(T(iper))'),Zxy(iper).*sqrt(T(iper))',order);
                %f_xy = polyval(p,log10(Tout),[],mu);
                f_xy = interp1((log10(T(iper))'),Zxy(iper).*sqrt(T(iper))',log10(Tout),method,'extrap');
            end
            subplot(2,2,2);
            semilogx(T,real(Zxy).*sqrt(T)','b*'); hold on;
            semilogx(Tout,real(f_xy),'g'); hold on;
            semilogx(T,imag(Zxy).*sqrt(T)','r*'); hold on;
            semilogx(Tout,imag(f_xy),'k'); hold on;
            ylabel('Zxy * sqrt(T)')
            
            if strcmp(method,'zeros')
                f_yx = zeros(1,length(Tout)) + 1i*zeros(1,length(Tout));
            else
                %[p,s,mu] = polyfit((log10(T(iper))'),Zyx(iper).*sqrt(T(iper))',order);
                %f_yx = polyval(p,log10(Tout),[],mu);
                f_yx = interp1((log10(T(iper))'),Zyx(iper).*sqrt(T(iper))',log10(Tout),method,'extrap');
            end
            subplot(2,2,3);
            semilogx(T,real(Zyx).*sqrt(T)','b*'); hold on;
            semilogx(Tout,real(f_yx),'g'); hold on;
            semilogx(T,imag(Zyx).*sqrt(T)','r*'); hold on;
            semilogx(Tout,imag(f_yx),'k'); hold on;
            ylabel('Zyx * sqrt(T)')
            
            if strcmp(method,'zeros')
                f_yy = zeros(1,length(Tout)) + 1i*zeros(1,length(Tout));
            else
                %[p,s,mu] = polyfit((log10(T(iper))'),Zyy(iper).*sqrt(T(iper))',order);
                %f_yy = polyval(p,log10(Tout),[],mu);
                f_yy = interp1((log10(T(iper))'),Zyy(iper).*sqrt(T(iper))',log10(Tout),method,'extrap');
            end
            subplot(2,2,4);
            semilogx(T,real(Zyy).*sqrt(T)','b*'); hold on;
            semilogx(Tout,real(f_yy),'g'); hold on;
            semilogx(T,imag(Zyy).*sqrt(T)','r*'); hold on;
            semilogx(Tout,imag(f_yy),'k'); hold on;
            ylabel('Zyy * sqrt(T)')
            
            Tnew = [Tout(Tout<T(1)) T Tout(Tout>T(end))];
            Zout = zeros(2,2,length(Tnew)) + 1i*zeros(2,2,length(Tnew));
            lnew1 = length(find(Tout<T(1)));
            lnew2 = length(find(Tout>T(end)));
            Zout(1,1,1:lnew1) = f_xx(Tout<T(1))./sqrt(Tout(Tout<T(1)));
            Zout(1,1,lnew1+1:end-lnew2) = Zxx;
            Zout(1,1,end-lnew2+1:end) = f_xx(Tout>T(end))./sqrt(Tout(Tout>T(end)));
            Zout(1,2,1:lnew1) = f_xy(Tout<T(1))./sqrt(Tout(Tout<T(1)));
            Zout(1,2,lnew1+1:end-lnew2) = Zxy;
            Zout(1,2,end-lnew2+1:end) = f_xy(Tout>T(end))./sqrt(Tout(Tout>T(end)));
            Zout(2,1,1:lnew1) = f_yx(Tout<T(1))./sqrt(Tout(Tout<T(1)));
            Zout(2,1,lnew1+1:end-lnew2) = Zyx;
            Zout(2,1,end-lnew2+1:end) = f_yx(Tout>T(end))./sqrt(Tout(Tout>T(end)));
            Zout(2,2,1:lnew1) = f_yy(Tout<T(1))./sqrt(Tout(Tout<T(1)));
            Zout(2,2,lnew1+1:end-lnew2) = Zyy;
            Zout(2,2,end-lnew2+1:end) = f_yy(Tout>T(end))./sqrt(Tout(Tout>T(end)));
            Tout = Tnew;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [apres,phase] = imp2apres(T,Z,Zstd)
            % T are the periods in seconds
            % Z is the impedance (2,2,nper)
            % Zstd is the standard error for real or imag parts (2,2,nper)
            %
            % apres apparent resistivity xx,xx_se,xy,xy_se,yx,yx_se,yy,yy_se
            % phase (-pi,pi] as in ModEM xx,xx_se,xy,xy_se,yx,yx_se,yy,yy_se
            %
            % rho = 0.2*T*abs(z)^2, where T is period is secs
            %     = 0.2*(abs(z)/sqrt(f))^2
            %
            % phase = arctan(imag(z)/real(z))
            %
            % see Vozoff 1972, Geophysics
            %
            % the error estimates are using so-called "delta method"
            % based on Tailor series expansion
            % ref Stuart & Ord 1994 Sec 10.5
            % see also Chave & Lezaeta 2007
            %
            % phases use atan2 which takes values (-pi,pi] as in ModEM
            % but we convert fully negative phases to fully positive
            % by adding 2*pi
            %
            % same computation as used in EMTF/matlab/ZPLT/ap_res.m
            % but decoded for readability!
            %
            % Corollary:
            %
            % Re(Z) = sqrt(5*rho/T) * cos(phase)
            % Im(Z) = sqrt(5*rho/T) * sin(phase)
            %
            
            rad_deg = 57.2958;
            
            % Off-diagonal components: Zxy
            value = squeeze(Z(1,2,:)).';
            error = squeeze(Zstd(1,2,:)).'; % standard error of the real or imag component
            apres.xy = (T/5.).*abs(value).^2; % actual rho, rescaled by period
            apres.xy_se = 2*error.*(T/5.).*abs(value);
            phase.xy = rad_deg*atan2(imag(value),real(value));
            phase.xy_se = rad_deg*abs(error./max(value,1e-20));
            
            % Off-diagonal components: Zyx
            value = squeeze(Z(2,1,:)).';
            error = squeeze(Zstd(2,1,:)).'; % standard error of the real or imag component
            apres.yx = (T/5.).*abs(value).^2; % actual rho, rescaled by period
            apres.yx_se = 2*error.*(T/5.).*abs(value);
            phase.yx = rad_deg*atan2(imag(value),real(value));
            phase.yx_se = rad_deg*abs(error./max(value,1e-20));
            
            % Diagonal components: Zxx
            value = squeeze(Z(1,1,:)).';
            error = squeeze(Zstd(1,1,:)).'; % standard error of the real or imag component
            apres.xx = (T/5.).*abs(value).^2; % actual rho, rescaled by period
            apres.xx_se = 2*error.*(T/5.).*abs(value);
            phase.xx = rad_deg*atan2(imag(value),real(value));
            phase.xx_se = rad_deg*abs(error./max(value,1e-20));
            
            % Diagonal components: Zyy
            value = squeeze(Z(2,2,:)).';
            error = squeeze(Zstd(2,2,:)).'; % standard error of the real or imag component
            apres.yy = (T/5.).*abs(value).^2; % actual rho, rescaled by period
            apres.yy_se = 2*error.*(T/5.).*abs(value);
            phase.yy = rad_deg*atan2(imag(value),real(value));
            phase.yy_se = rad_deg*abs(error./max(value,1e-20));
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Z,Zstd] = apres2imp(T,apres,phase)
            % T are the periods in seconds
            % Z is the impedance (2,2,nper)
            % Zstd is the standard error for real or imag parts (2,2,nper)
            %
            % apres apparent resistivity xx,xx_se,xy,xy_se,yx,yx_se,yy,yy_se
            % phase (-pi,pi] as in ModEM xx,xx_se,xy,xy_se,yx,yx_se,yy,yy_se
            %
            % reverse of imp2apres for use with Chave's BIRRP files
            %
            % rho = 0.2*T*abs(z)^2, where T is period is secs
            %     = 0.2*(abs(z)/sqrt(f))^2
            %
            % phase = arctan(imag(z)/real(z))
            %
            % see Vozoff 1972, Geophysics
            %
            % the error estimates are using so-called "delta method"
            % based on Tailor series expansion
            % ref Stuart & Ord 1994 Sec 10.5
            % see also Chave & Lezaeta 2007
            %
            % phases use atan2 which takes values (-pi,pi] as in ModEM
            %
            % Corollary:
            %
            % Re(Z) = sqrt(5*rho/T) * cos(phase)
            % Im(Z) = sqrt(5*rho/T) * sin(phase)
            %

            Z = NaN*(zeros(2,2,length(T)) + 1i*zeros(2,2,length(T)));
            Zstd = NaN*zeros(2,2,length(T));
            
            % Off-diagonal components: Zxy
            value = sqrt(5*apres.xy./T') .* (cosd(phase.xy) + 1i*sind(phase.xy));
            error =  (5*apres.xy_se./T') ./ (2*abs(value)); 
            Z(1,2,:) = value;
            Zstd(1,2,:) = error; % standard error of the real or imag component

            % Off-diagonal components: Zyx
            value = sqrt(5*apres.yx./T') .* (cosd(phase.yx) + 1i*sind(phase.yx));
            error =  (5*apres.yx_se./T') ./ (2*abs(value)); 
            Z(2,1,:) = value;
            Zstd(2,1,:) = error; % standard error of the real or imag component

            if isfield(apres,'xx')
                
            % Diagonal components: Zxx
            value = sqrt(5*apres.xx./T') .* (cosd(phase.xx) + 1i*sind(phase.xx));
            error =  (5*apres.xx_se./T') ./ (2*abs(value)); 
            Z(1,1,:) = value;
            Zstd(1,1,:) = error; % standard error of the real or imag component

            % Diagonal components: Zyy
            value = sqrt(5*apres.yy./T') .* (cosd(phase.yy) + 1i*sind(phase.yy));
            error =  (5*apres.yy_se./T') ./ (2*abs(value)); 
            Z(2,2,:) = value;
            Zstd(2,2,:) = error; % standard error of the real or imag component
            
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [pol] = test_electrode_polarity(phase)
            % For mid-range and long periods, verifies that XY phase is in 
            % 1st quadrant and YX phase is in 3rd quadrant. Outputs are zero
            % if the phases pass the test; otherwise electrode polarity
            % may be reversed. Assuming atan2 was used in imp2apres.
            
            nper = length(phase.xy);
            if nper > 10
                ip = 10:nper;
            else
                ip = 1:nper;
            end
            pol = [0 0];
            if ~any(phase.xy(ip) > 0 & phase.xy(ip) < 90) % should be in 1st quadrant
                pol(1) = 1;
            end
            if ~any(phase.yx(ip) > -270 & phase.yx(ip) < -90) % should be in 3rd quadrant
                pol(2) = 1;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [periods,data,std,latlon,metadata,copyright,SIG_S,SIG_E] = read_xml(cfile,dataType)
            % [periods,data,std,latlon,metadata,copyright,SIG_S,SIG_E] = read_xml(cfile,dataType)
            %
            % Usage: [periods,Z,Zstd,metadata,copyright,SIG_S,SIG_E] = read_xml(cfile,1)
            %        [periods,T,Tstd,metadata,copyright,SIG_S,SIG_E] = read_xml(cfile,2)
            %
            % Reads the EMTF XML format data file. Uses xml_io_tools.
            % Only reads the impedance Z and the vertical field TFs T
            % from the file, or optionally one of these.
            %
            % dataType = 1 for Full_Impedance
            %          = 2 for Full_Vertical_Components
            % by default, read the impedance.
            %
            % metadata is a structure that includes:
            %   SITEID=AK1A
            %   PROJECT=1D Synthetic response
            %   SURVEY=1D 3D comp
            %   YEAR=2017
            %   ORIENT=[0 90 0 0 90]
            %   SIGNCONVENTION=exp(+ i\omega t)
            %   PROCESSEDBY=Paul Bedrosian
            %   PROCESSINGSOFTWARE=MT1D_pick
            %   PROCESSINGTAG=AK1A
            %   SITENAME=Adirondack Mountains, NY, USA
            %   RUNLIST=synth
            %   REMOTEREF=Robust Remote Reference
            %   REMOTESITE=NA
            %   UNITS='[mV/km]/[nT]'
            %   QUALITYRATING=5
            %   GOODFROMPERIOD=10
            %   GOODTOPERIOD=1e5
            %   QUALITYCOMMENT=''
            %   WARNINGFLAG=0
            %   WARNINGCOMMENT=''
            %   PREDICTED=0
            %
            % Note: for simplicity, this code uses the "name" attribute of
            % the transfer functions in the XML file. In general we only
            % use it for convenience, it could potentially even be missing
            % from file. If generality is needed in the future, will work
            % on that.
            %
            % Could and should also read in and make use of: units, sign
            % convention, data quality, and period range.
                        
            % by default, read metadata only - much faster!
            typeMismatch = 0;
            if nargin < 2
                dataType = 0;
                typeMismatch = 1;
                type = 'Z';
                nchout = 2;
            end
            
            % ingest the XML
            Pref.Str2Num = 'never';
            info = xml_read(cfile,Pref);
           
            if dataType == 1
                if isempty(strfind(info.Tags,'impedance'))
                    typeMismatch = 1;
                end
                type = 'Z';
                nchout = 2;
            elseif dataType == 2
                if isempty(strfind(info.Tags,'tipper'))
                    typeMismatch = 1;
                end
                type = 'T';
                nchout = 1;
            end
            
            element = [type '.VAR']; variance = strrep(element,'.','0x2E');
            element = [type '.INVSIGCOV']; invsigcov = strrep(element,'.','0x2E');
            element = [type '.RESIDCOV']; residcov = strrep(element,'.','0x2E');
            
            nper = length(info.Data.Period);
            periods = zeros(1,nper);
            data = (NaN + 1i*NaN) * zeros(2,nchout,nper);
            std = NaN * zeros(2,nchout,nper);
            latlon(1) = str2num(info.Site.Location.Latitude);
            latlon(2) = str2num(info.Site.Location.Longitude);
            latlon(3) = str2num(info.Site.Location.Elevation.CONTENT);
            if ~ strfind(info.Site.Location.Elevation.ATTRIBUTE.units,'meters')
                disp(['Warning: elevation is given in ' info.Site.Location.Elevation.ATTRIBUTE.units]);
            end
            if nargout > 4
               metadata.SITEID = info.Site.Id;
               metadata.SITENAME = info.Site.Name;
               metadata.PROJECT = info.Site.Project;
               metadata.SURVEY = info.Site.Survey;
               metadata.YEAR = info.Site.YearCollected;
               metadata.PREDICTED = 0;
               metadata.ORTHOGONAL = 0;
               metadata.THETA0 = 0.0;
               if isfield(info,'SiteLayout') % new XML format as of 2017
                   InputChannels = info.SiteLayout.InputChannels;
                   OutputChannels = info.SiteLayout.OutputChannels;
                   if ischar(info.Site.Orientation)
                       Orientation = info.Site.Orientation;
                   else
                       Orientation = info.Site.Orientation.CONTENT;
                       metadata.THETA0 = str2double(info.Site.Orientation.ATTRIBUTE.angle_to_geographic_north);
                   end
                   if strcmpi(Orientation,'Site Layout') || strcmpi(Orientation,'SiteLayout')
                       metadata.ORTHOGONAL = 0;
                       disp('Warning: use original site layout for TF orientation');
                   else
                       metadata.ORTHOGONAL = 1; 
                       disp('Warning: data already rotated to an orthogonal coordinate system. Ignore site layout');
                   end
               else % needed for backwards compatibility with old XML files
                   InputChannels = info.InputChannels;
                   OutputChannels = info.OutputChannels;
               end
               for i = 1:length(InputChannels.Magnetic)
                   if strcmpi(InputChannels.Magnetic(i).ATTRIBUTE.name(1:2),'Hx')
                       metadata.ORIENT(1) = str2num(InputChannels.Magnetic(i).ATTRIBUTE.orientation);
                   end
                   if strcmpi(InputChannels.Magnetic(i).ATTRIBUTE.name(1:2),'Hy')
                       metadata.ORIENT(2) = str2num(InputChannels.Magnetic(i).ATTRIBUTE.orientation);
                   end
               end
               if isfield(OutputChannels,'Magnetic')
                   for i = 1:length(OutputChannels.Magnetic)
                       if strcmpi(OutputChannels.Magnetic(i).ATTRIBUTE.name(1:2),'Hz')
                           metadata.ORIENT(3) = str2num(OutputChannels.Magnetic(i).ATTRIBUTE.orientation);
                       end
                   end
               else
                   metadata.ORIENT(3) = 0;
               end
               for i = 1:length(OutputChannels.Electric)
                   if strcmpi(OutputChannels.Electric(i).ATTRIBUTE.name(1:2),'Ex')
                       metadata.ORIENT(4) = str2num(OutputChannels.Electric(i).ATTRIBUTE.orientation);
                   end
                   if strcmpi(OutputChannels.Electric(i).ATTRIBUTE.name(1:2),'Ey')
                       metadata.ORIENT(5) = str2num(OutputChannels.Electric(i).ATTRIBUTE.orientation);
                   end
               end
               if metadata.ORTHOGONAL % ignore site layout!
                   metadata.ORIENT = [0 90 0 0 90];
               end
               % this can be generalized to death but I won't do it
               if strcmp(type,'T')
                   metadata.ORIENT = metadata.ORIENT(1:3);
               elseif strcmp(type,'Z')
                   metadata.ORIENT = metadata.ORIENT([1:2 4:5]);
               end
               if isfield(info,'ProcessingInfo')
                   metadata.SIGNCONVENTION = info.ProcessingInfo.SignConvention;
                   if isfield(info.ProcessingInfo,'ProcessedBy')
                       metadata.PROCESSEDBY = info.ProcessingInfo.ProcessedBy;
                       metadata.PROCESSINGSOFTWARE = info.ProcessingInfo.ProcessingSoftware.Name;
                       metadata.PROCESSINGTAG = info.ProcessingInfo.ProcessingTag;
                       if isfield(info,'FieldNotes')
                           tmp = info.FieldNotes;
                       else
                           tmp = '';
                       end
                       metadata.RUNLIST = '';
                       for i=1:length(tmp); metadata.RUNLIST = [metadata.RUNLIST ' ' tmp(i).ATTRIBUTE.run]; end
                       metadata.REMOTEREF = info.ProcessingInfo.RemoteRef.ATTRIBUTE.type;
                   end
                   if isfield(info.ProcessingInfo,'RemoteInfo')
                       metadata.REMOTESITE = info.ProcessingInfo.RemoteInfo.Site.Id;
                   end
               else
                   % assume +ve sign convention by default
                   metadata.SIGNCONVENTION = 'exp(+ i\omega mu)';
               end
               if isfield(info.Site,'DataQualityNotes')
                   metadata.QUALITYRATING = info.Site.DataQualityNotes.Rating;
                   if str2num(metadata.QUALITYRATING) == 0
                       metadata.GOODFROMPERIOD = 0;
                       metadata.GOODTOPERIOD = 0;
                       metadata.QUALITYCOMMENTS = 'Unrated';
                   else
                       metadata.GOODFROMPERIOD = str2num(info.Site.DataQualityNotes.GoodFromPeriod);
                       metadata.GOODTOPERIOD = str2num(info.Site.DataQualityNotes.GoodToPeriod);
                       metadata.QUALITYCOMMENTS = '';
                       if isfield(info.Site.DataQualityNotes,'Comments')
                           metadata.QUALITYCOMMENTS = info.Site.DataQualityNotes.Comments.CONTENT;
                       end
                   end
               end
               if isfield(info.Site,'DataQualityWarnings')
                   metadata.WARNINGFLAG = info.Site.DataQualityWarnings.Flag;
                   metadata.WARNINGCOMMENTS = '';
                   if isfield(info.Site.DataQualityWarnings,'Comments')
                       metadata.WARNINGCOMMENTS = info.Site.DataQualityWarnings.Comments.CONTENT;
                   end
               end
               if ~typeMismatch
                metadata.UNITS = info.Data.Period(1).(type).ATTRIBUTE.units;
               end
            end
            if nargout > 5
                if isfield(info,'Copyright')
                    copyright = info.Copyright;
                else
                    copyright = mttf.newcopyright;
                end
            end
            if nargout > 6
                SIG_S = (NaN + 1i*NaN) * zeros(2,2,nper);
                SIG_E = (NaN + 1i*NaN) * zeros(nchout,nchout,nper);
            end
            
            for j = 1:length(info.Data.Period)
                periods(j) = str2num(info.Data.Period(j).ATTRIBUTE.value);
                if typeMismatch; continue; end % skip the rest of the loop if type not present
                element = info.Data.Period(j).(type);
                if ~isfield(element,'value')
                    warning([info.Site.Id ': Empty element ' type ' for period #' num2str(j)]);
                else
                    for i = 1:length(element.value)
                        if isnumeric(element.value(i).CONTENT)
                            value = element.value(i).CONTENT;
                        else
                            value = str2num(element.value(i).CONTENT); %#ok<*ST2NM>
                        end
                        if strcmp(type,'Z')
                            switch upper(element.value(i).ATTRIBUTE.name)
                                case 'ZXX'
                                    data(1,1,j) = value(1) + 1i*value(2);
                                case 'ZXY'
                                    data(1,2,j) = value(1) + 1i*value(2);
                                case 'ZYX'
                                    data(2,1,j) = value(1) + 1i*value(2);
                                case 'ZYY'
                                    data(2,2,j) = value(1) + 1i*value(2);
                                otherwise
                                    % do nothing
                            end
                        elseif strcmp(type,'T')
                            switch upper(element.value(i).ATTRIBUTE.name)
                                case 'TX'
                                    data(1,1,j) = value(1) + 1i*value(2);
                                case 'TY'
                                    data(2,1,j) = value(1) + 1i*value(2);
                                otherwise
                                    % do nothing
                            end
                        else
                            warning(['mttf.read_xml is not yet coded to read in the data type ' element.value.ATTRIBUTE.name]);
                        end
                    end
                end       
                if isfield(info.Data.Period(j),variance)
                    element = info.Data.Period(j).(variance);
                    if ~isfield(element,'value')
                        warning([info.Site.Id ': Empty variance for ' type ' for period #' num2str(j)]);
                    else
                        for i = 1:length(element.value)
                            if isnumeric(element.value(i).CONTENT)
                                value = element.value(i).CONTENT;
                            else
                                value = str2num(element.value(i).CONTENT); %#ok<*ST2NM>
                            end
                            if strcmp(type,'Z')
                                switch upper(element.value(i).ATTRIBUTE.name)
                                    case 'ZXX'
                                        std(1,1,j) = sqrt(value/2);
                                    case 'ZXY'
                                        std(1,2,j) = sqrt(value/2);
                                    case 'ZYX'
                                        std(2,1,j) = sqrt(value/2);
                                    case 'ZYY'
                                        std(2,2,j) = sqrt(value/2);
                                    otherwise
                                        % do nothing
                                end
                            elseif strcmp(type,'T')
                                switch upper(element.value(i).ATTRIBUTE.name)
                                    case 'TX'
                                        std(1,1,j) = sqrt(value/2);
                                    case 'TY'
                                        std(2,1,j) = sqrt(value/2);
                                    otherwise
                                        % do nothing
                                end
                            else
                                warning(['mttf.read_xml is not yet coded to read in the data type ' element.value.ATTRIBUTE.name]);
                            end
                        end
                    end
                else
                    std(1:2,1:2,j) = NaN;
                    metadata.PREDICTED = 1;
                end
                if nargout > 6
                    %read INVSIGCOV into SIG_S and RESIDCOV into SIG_E
                    SIG_S = info.Data.Period(j).(invsigcov);
                    SIG_E = info.Data.Period(j).(residcov);
                end
            end
            
            data = squeeze(data);
            std = squeeze(std);
            if nargout > 6
                SIG_S = squeeze(SIG_S);
                SIG_E = squeeze(SIG_E);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [periods,data,std,latlon,metadata] = read_edi(cfile,dataType)
            % [periods,data,std,latlon,metadata] = read_edi(cfile,dataType)
            %
            % Usage: [periods,Z,Zstd,latlon,metadata] = read_edi(cfile,1)
            %        [periods,T,Tstd,latlon,metadata] = read_edi(cfile,2)
            %
            % Based on a simple code snippet from Paul Bedrosian.
            % Reads the EDI format data file.
            % Only reads the impedance Z and the vertical field TFs T
            % from the file, or optionally one of these.
            %
            % dataType = 1 for Full_Impedance
            %          = 2 for Full_Vertical_Components
            % by default, read the impedance.
            %
            % metadata is a structure that may include:
            %   SITEID=AK1A
            %   ORIENT=[0 90 0 0 90]
            %   PROJECT=1D Synthetic response
            %   SURVEY=1D 3D comp
            %   YEAR=2017
            %   PROCESSEDBY=Paul Bedrosian
            %   PROCESSINGSOFTWARE=MT1D_pick
            %   PROCESSINGTAG=AK1A
            %   SITENAME=Adirondack Mountains, NY, USA
            %   RUNLIST=synth
            %   REMOTEREF=Robust Remote Reference
            %   REMOTESITE=NA
            %   UNITS='[mV/km]/[nT]'
            %   SIGNCONVENTION=exp(+ i\omega t)
            %   QUALITYRATING=5
            %   GOODFROMPERIOD=10
            %   GOODTOPERIOD=1e5
            %   QUALITYCOMMENT=''
            %   WARNINGFLAG=0
            %   WARNINGCOMMENT=''
            %   ORTHOGONAL=1
            %
            % theta0 is channel orientation relative to geographic North
            %
            % the output should always have +ve sign convention assuming
            % that there is enough information in the file.
            
            dummy=1e+32; % dummy data value to be excluded.
                                    
            if nargin < 2
                dataType = 1;
            end
           
            if dataType == 1
                type = 'Z';
                nchout = 2;
            elseif dataType == 2
                type = 'T';
                nchout = 1;
            end
            
            if strfind(cfile,'.edi')
                fid=fopen(cfile,'r');
            else
                fid=fopen(strcat(cfile,'.edi'),'r');
            end
            
            tline=fgetl(fid);
            while (isempty(strfind(tline,'DATAID')))
                tline=fgetl(fid);
            end
            tmp = sscanf(tline,'%s');
            DAT.siteid = tmp(9:end-1);

            while (isempty(strfind(tline,'>=DEFINEMEAS'))) % read in all metadata
                tline=fgetl(fid);
                i = strfind(tline,'=');
                field = strtrim(tline(1:i-1));
                info = strtrim(tline(i+1:end));
                if ~isempty(field) && ~strcmp(field,'>') && isempty(strfind(field,' '))
                    metadata.(field) = info; 
                end
            end
            
            signconvention = 1; % default is +ve
            if isfield(metadata,'SIGNCONVENTION')
                if strfind(metadata.SIGNCONVENTION,'+')
                    signconvention = 1;
                elseif strfind(metadata.SIGNCONVENTION,'-')
                    signconvention = -1;
                end
            end

            while (isempty(strfind(tline,'REFLAT')) && (isempty(strfind(tline,'LATITUDE'))))  % find out the site latitude
                tline=fgetl(fid);
            end
            if (~isempty(strfind(tline,'REFLAT'))) % Use for edi files that supply reflat/reflon as dms.
                tline = strtrim(tline);
                ind=strfind(tline,'=');
                tmp=sscanf(tline(ind+1:end),'%s');
                tline=fgetl(fid);
                tline = strtrim(tline);
                ind=strfind(tline,'=');
                tmp2=sscanf(tline(ind+1:end),'%s');
                
                dum=sscanf(tmp,'%d:%d:%f');
                dum2=sscanf(tmp2,'%d:%d:%f');
                DAT.stdec = [sign(dum(1))*(abs(dum(1))+dum(2)/60+dum(3)/3600) sign(dum2(1))*(abs(dum2(1))+dum2(2)/60+dum2(3)/3600)];
            end
            
            while (isempty(strfind(tline,'NFREQ'))) % find out how many freqs we have
                tline=fgetl(fid);
            end
            tmp=sscanf(tline,'%s');
            
            nf=str2num(tmp(7:end));

            tline=fgetl(fid);
            while (isempty(strfind(tline,'>FREQ'))) % read in the frequencies
                tline=fgetl(fid);
            end
            DAT.fr=fscanf(fid,'%e\n',nf)';
            
            tline=fgetl(fid);            
            while (isempty(strfind(tline,'>ZROT'))) % read in the first rotation value
                tline=fgetl(fid);
            end
            DAT.theta0=fscanf(fid,'%e',1);
    
            tline=fgetl(fid);
            while (isempty(strfind(tline,'>ZXXR'))) % read in the impedances
                tline=fgetl(fid);
            end
            for j=1:12
                DAT.z(j,:)=fscanf(fid,'%e\n',nf)';
                tline=fgetl(fid);
            end
            
            % not all sites have tipper data
            while (isempty(strfind(tline,'>TXR'))) % read in the tipper data
                tline=fgetl(fid);
            end
            for j=1:6
                DAT.t(j,:)=fscanf(fid,'%e\n',nf)';
                fgetl(fid);
            end
            
            fclose(fid);
            
            % do some clean up within the site loop.
            tmp=DAT.z; % replace any dummy Z values with NaN
            tmp(abs(tmp)==dummy)=NaN;
            DAT.z=tmp;
            tmp=DAT.t; % replace any dummy T values with NaN
            tmp(abs(tmp)==dummy)=NaN;
            DAT.t=tmp;
            if(DAT.fr(1) <= DAT.fr(2)) % EDI files should run from high to low frequency
                DAT.fr=fliplr(DAT.fr);
                DAT.z=fliplr(DAT.z);
                DAT.t=fliplr(DAT.t);
            end
            if (max(abs(DAT.t(1,:))) <= 1e-2) % we have effectively zero tipper, set to identically zero so flagging works
                DAT.t=zeros(6,size(DAT.t,2));
            end
            
            periods = 1./DAT.fr;
            latlon = DAT.stdec;
            data = (NaN + 1i*NaN) * zeros(2,nchout,nf);
            std = NaN * zeros(2,nchout,nf);

            if strcmp(type,'Z')
                for j = 1:nf
                    data(1,1,j) = DAT.z(1,j) + 1i*DAT.z(2,j);
                    data(1,2,j) = DAT.z(4,j) + 1i*DAT.z(5,j);
                    data(2,1,j) = DAT.z(7,j) + 1i*DAT.z(8,j);
                    data(2,2,j) = DAT.z(10,j) + 1i*DAT.z(11,j);
                    std(1,1,j) = sqrt(DAT.z(3,j)/2);
                    std(1,2,j) = sqrt(DAT.z(6,j)/2);
                    std(2,1,j) = sqrt(DAT.z(9,j)/2);
                    std(2,2,j) = sqrt(DAT.z(12,j)/2);
                end
            elseif strcmp(type,'T')
                 for j = 1:nf
                    data(1,1,j) = DAT.z(1,j) + 1i*DAT.z(2,j);
                    data(1,2,j) = DAT.z(4,j) + 1i*DAT.z(5,j);
                    data(2,1,j) = DAT.z(7,j) + 1i*DAT.z(8,j);
                    data(2,2,j) = DAT.z(10,j) + 1i*DAT.z(11,j);
                    std(1,1,j) = sqrt(DAT.z(3,j)/2);
                    std(1,2,j) = sqrt(DAT.z(6,j)/2);
                    std(2,1,j) = sqrt(DAT.z(9,j)/2);
                    std(2,2,j) = sqrt(DAT.z(12,j)/2);
                 end
            end
            
            metadata.SITEID = DAT.siteid;
            metadata.ORIENT = [0 90 0 0 90];
            metadata.THETA0 = DAT.theta0;
            metadata.ORTHOGONAL = 1; % assume true for EDIs
            if signconvention < 0 % take complex conjugate of all  
                data = conj(data);
            end
            metadata.SIGNCONVENTION = 'exp(+ i\omega t)';
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [T,apres,phase] = read_birrp(cfile)
            % Reads output of BIRRP robust procedure by Alan Chave.
            %
            % reads specified file name and all related files to output the
            % complete apres and phase structures for the 4 impedance components;
            % can deal with missing values. sorts by increasing period.
            %
            % <file name> xxx.mtres.all.*r#.rp  (xxx=kak, kny, or mmb)
            % * and # are either of 1 or 2.  If *=1, it means Ex (the northward component
            % of the geoelectric field).  *=2 corresponds to Ey (the eastward component
            % of the geoelectric field).  # is for the geomagnetic field as the same numbering as
            % the electric field.  So, if a file name includes '1r2', it means the response for Ex/By.
            %
            % <format> Each line contains data as, from the left column to the right,
            % period (second),  frequency (Hz), apparent resistivity (ohm m), 95% confidence
            % limit of the apparent resistivity (ohm m), phase (degree), 95% confidence
            % limit of the phase (degree), the equivalent degrees of freedom, the fraction
            % of the data used after coherence thresholding (should be 1.0), the fraction
            % of the data used after final bounded influence weighting, and
            % an integer convergence flag (should be 0 which means proper convergence)            
            
            fileroot = cfile(1:end-7);
            if ~exist([fileroot '.1r1.rp'],'file')
                error(['BIRRP data file ' fileroot ' for Zxx not found']);
            else
                data{1,1} = load([fileroot '.1r1.rp'],'-ascii');
            end
            if ~exist([fileroot '.1r2.rp'],'file')
                error(['BIRRP data file ' fileroot ' for Zxy not found']);
            else
                data{1,2} = load([fileroot '.1r2.rp'],'-ascii');
            end
            if ~exist([fileroot '.2r1.rp'],'file')
                error(['BIRRP data file ' fileroot ' for Zyx not found']);
            else
                data{2,1} = load([fileroot '.2r1.rp'],'-ascii');
            end
            if ~exist([fileroot '.2r2.rp'],'file')
                error(['BIRRP data file ' fileroot ' for Zyy not found']);
            else
                data{2,2} = load([fileroot '.2r2.rp'],'-ascii');
            end
            
            T = unique([data{1,1}(:,1); data{1,2}(:,1); data{2,1}(:,1); data{2,2}(:,1)]);
            rho = NaN*(zeros(2,2,length(T)));
            rho_se = NaN*zeros(2,2,length(T));
            phs = NaN*zeros(2,2,length(T));
            phs_se = NaN*zeros(2,2,length(T));
            for i = 1:2
                for j = 1:2
                    [~,iper,l]=intersect(T,data{i,j}(:,1));
                    for k = iper
                        rho(i,j,k) = data{i,j}(l,3);
                        rho_se(i,j,k) = data{i,j}(l,4);
                        phs(i,j,k) = data{i,j}(l,5);
                        phs_se(i,j,k) = data{i,j}(l,6);
                    end
                end
            end
            
            apres.xy = rho(1,2,:);
            apres.xy_se = rho_se(1,2,:);
            phase.xy = phs(1,2,:);
            phase.xy_se = phs_se(1,2,:);
            
            apres.yx = rho(2,1,:);
            apres.yx_se = rho_se(2,1,:);
            phase.yx = phs(2,1,:);
            phase.yx_se = phs_se(2,1,:);
            
            apres.xx = rho(1,1,:);
            apres.xx_se = rho_se(1,1,:);
            phase.xx = phs(1,1,:);
            phase.xx_se = phs_se(1,1,:);
            
            apres.yy = rho(2,2,:);
            apres.yy_se = rho_se(2,2,:);
            phase.yy = phs(2,2,:);
            phase.yy_se = phs_se(2,2,:);
            
            for fn = fieldnames(apres)'
                apres.(fn{1}) = squeeze(apres.(fn{1}));
            end

            for fn = fieldnames(phase)'
                phase.(fn{1}) = squeeze(phase.(fn{1}));
            end
            
            T = T';

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [impulse] = read_impulse(cfile)
            
            % [impulse] = read_impulse(cfile)
            %
            % The impulse is a structure with fields xx,xy,yx,yy
            % and some TF related and damping parameter metadata

            if exist(cfile,'file')
                % read from file, header first
                fid = fopen(cfile);
                tline = fgetl(fid);
                impulse.tfname = sscanf(tline,'%% TF name: %s');
                tline = fgetl(fid);
                location = sscanf(tline,'%% Geographic location: %f %f %f');
                impulse.lat = location(1);
                impulse.lon = location(2);
                if length(location) > 2
                    impulse.elev = location(3);
                else
                    impulse.elev = 0;
                end
                tline = fgetl(fid);
                impulse.orient = sscanf(tline,'%% Geographic orientation: %f\n');
                tline = fgetl(fid);
                impulse.delta_t = sscanf(tline,'%% Sampling rate (secs): %f\n');
                tline = fgetl(fid);
                impulse.nbounds = sscanf(tline,'%% Bounds: %d %d\n');
                tline = fgetl(fid);
                params = sscanf(tline,'%% Damping parameters: %d %e %d\n');
                impulse.option = params(1);
                impulse.lambda = params(2);
                impulse.errorbars = params(3);
                tline = fgetl(fid);
                format = sscanf(tline,'%% Format: %13c\n');
                if ~strcmp(format,'n xx xy yx yy')
                    error('Unknown format of impulse file.');
                end
                % read the impulse
                nlines = fscanf(fid,'%f\n',1);
                temp = textscan(fid,'%f %f %f %f %f',nlines,'delimiter','\n');
                fclose(fid);
                % save the impulse
                impulse.xx = temp{2};
                impulse.xy = temp{3};
                impulse.yx = temp{4};
                impulse.yy = temp{5};
            else
                error('Impulse file not found');
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = write_impulse(cfile,impulse)
            
            % [] = write_impulse(cfile,impulse)
            %
            % The impulse is a structure with fields xx,xy,yx,yy
            
            nmin = impulse.nbounds(1);
            nmax = impulse.nbounds(2);

            % Check that everything is self-consistent
            if nmax-nmin+1 ~= length(impulse.xx)
                error('The impulse length does not match the integer bounds');
            end
            
            % Fix orientation if needed
            if isempty(impulse.orient)
                orientation = 0;
            else
                orientation = impulse.orient;
            end
            
            % write to file, header first
            fid = fopen(cfile,'w');
            fprintf(fid,'%% TF name: %s\n',impulse.tfname);
            fprintf(fid,'%% Geographic location: %f %f %f\n',impulse.lat,impulse.lon,impulse.elev);
            fprintf(fid,'%% Geographic orientation: %f\n',orientation);
            fprintf(fid,'%% Sampling rate (secs): %f\n',impulse.delta_t);
            fprintf(fid,'%% Bounds: %d %d\n',impulse.nbounds);
            fprintf(fid,'%% Damping parameters: %d %e %d\n',impulse.option,impulse.lambda,impulse.errorbars);
            fprintf(fid,'%% Format: n xx xy yx yy\n');
            fprintf(fid,'%d\n',nmax-nmin+1);
            for i = 1:nmax-nmin+1
                fprintf(fid,'%d %15.9e %15.9e %15.9e %15.9e\n',nmin+i-1,impulse.xx(i),impulse.xy(i),impulse.yx(i),impulse.yy(i));
            end
            fclose(fid);
                    
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = plot_impulse(impulse,nbounds)
            
            % [] = plot_impulse(impulse,nbounds)
            %
            % impulse is a DTIR structure
            % see also read_impulse, write_impulse
            
            if nargin < 2
                nbounds = impulse.nbounds;
            end
            nmin = nbounds(1);
            nmax = nbounds(2);
            
            i1 = nmin - impulse.nbounds(1) + 1;
            i2 = nmax - impulse.nbounds(1) + 1;

            % colors consistent with apresplt
            symbol_colors = [ 0 0 .7; .7 0 0; 0 .7 0; .35 0 .35; 0 .35 .35 ; 0 0 0 ; .3 .3 0];
            
            figure;
            plot(nmin:nmax,squeeze(impulse.xx(i1:i2)),'color',symbol_colors(3,:),'linewidth',2); hold on;
            plot(nmin:nmax,squeeze(impulse.xy(i1:i2)),'color',symbol_colors(1,:),'linewidth',2); hold on;
            plot(nmin:nmax,squeeze(impulse.yx(i1:i2)),'color',symbol_colors(2,:),'linewidth',2); hold on;
            plot(nmin:nmax,squeeze(impulse.yy(i1:i2)),'color',symbol_colors(4,:),'linewidth',2); hold on;
%             plot(nmin:nmax,squeeze(impulse.xx(i1:i2)),'g','linewidth',2); hold on;
%             plot(nmin:nmax,squeeze(impulse.xy(i1:i2)),'r','linewidth',2); hold on;
%             plot(nmin:nmax,squeeze(impulse.yx(i1:i2)),'b','linewidth',2); hold on;
%             plot(nmin:nmax,squeeze(impulse.yy(i1:i2)),'k','linewidth',2); hold on;
            legend({'XX','XY','YX','YY'},'location','best','fontsize',18,'fontweight','bold');
            ylabel('mV / km / nT','fontsize',22,'fontweight','bold');
            xlabel('Time (samples)','fontsize',22,'fontweight','bold');
            set(gca,'fontsize',18,'fontweight','bold');

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = aprespltall(cdir,cext,period_range)
            
            % [] = aprespltall(cdir,ext,period_range)
            %
            % Same as EMTF/matlab/ZPLT/aprespltall.m except it does not cd
            % into the file directory. Plots in the file directory.
            %
            % Plots all files in cdir with extension = cext (e.g., xml).
            % Note: when reading XML files, often need to use absolute
            % path, or else the Matlab xmlread fails.
            %
            % By default plots Z-files.
            
            if nargin == 0
                cdir = uigetdir;
            end
            
            if nargin < 2
                cext = 'zss,zrr,zmm';
            end
            
            if nargin < 3
                period_range = [1 1e5];
            end
            
            comp=computer;
            if strcmp(comp(1:4),'PCWI')
                zid=strrep(zid,'/','\');
                cdir=strrep(cdir,'/','\');
            end
            
            cdir = [cdir '/'];
            filelist = dir(cdir);
            for i = 1:length(filelist)
                cfile = filelist(i).name;
                [~,fileroot,ext] = fileparts(cfile);
                if ~isempty(ext); ext = ext(2:end); end
                if ~isempty(strfind(cext,ext))
                    obj = mttf.read([cdir cfile]);
                    obj.apresplt(period_range,[1e-1 1e3]); %opb.impplt(2);
                    print(gcf,'-dpng',[cdir fileroot]);
                    close(gcf)
                end
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function result = isComplex(dataType)
            %  Usage : result = isComplex(dataType)
            %
            % gives a yes or no answer for complexity of a data type based
            % on a ModEM name
            
            switch dataType
                case {'Full_Impedance','Off_Diagonal_Impedance','Full_Vertical_Components'}
                    result = 1;
                case 'Full_Interstation_TF'
                    result = 1;
                case {'Off_Diagonal_Rho_Phase','Phase_Tensor'}
                    result = 0;
                case {'Ex_Field', 'Ey_Field', 'E_Field', 'Bx_Field', 'By_Field', 'Bz_Field', 'B_Field'}
                    result = 1;
                case 'Pole_Pole_DC_Rho'
                    result = 0;
                otherwise
                    error('This data type is not currently known in ModEM');
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [nComp,compChar,units,ochid] = components(dataType)
            %  Usage : [nComp,compChar,units,ochid] = components(dataType)
            %
            % specifies the number and names of components for each valid
            % ModEM data type (real and imaginary values are considered
            % separate components for counting)
            %
            % by default we always specify practical units, not SI units
            %
            % also include output channel IDs
            
            switch dataType
                case 'Full_Impedance'
                    nComp = 8;
                    compChar = ['ZXX';'ZXY';'ZYX';'ZYY'];
                    units = '[mV/km]/[nT]';
                    ochid = {'Ex' 'Ey'};
                case 'Off_Diagonal_Impedance'
                    nComp = 4;
                    compChar = ['ZXY';'ZYX'];
                    units = '[mV/km]/[nT]';
                    ochid = {'Ex' 'Ey'};
                case 'Full_Vertical_Components'
                    nComp = 4;
                    compChar = ['TX ';'TY '];
                    units = '[]';
                    ochid = {'Hz'};
                case 'Full_Interstation_TF'
                    nComp = 8;
                    compChar = ['MXX';'MXY';'MYX';'MYY'];
                    units = '[]';
                    ochid = {'Rx' 'Ry'};
                case 'Off_Diagonal_Rho_Phase'
                    nComp = 4;
                    compChar = ['RHOXY';'PHSXY';'RHOYX';'PHSYX'];
                    units = '[]';
                    ochid = {'Ex' 'Ey'};
                case 'Phase_Tensor'
                    nComp = 4;
                    compChar = ['PTXX';'PTXY';'PTYX';'PTYY'];
                    units = '[]';
                    ochid = {'Ex' 'Ey'};
                case 'Pole_Pole_DC_Rho'
                    nComp = 1;
                    compChar = 'RHODC';
                    units = '[]';
                    ochid = {''};
                case 'E_Field'
                    nComp = 4;
                    compChar = ['EX';'EY'];
                    units = '[mV/km]';
                    ochid = {'Ex' 'Ey'};
                case 'B_Field'
                    nComp = 6;
                    compChar = ['BX';'BY';'BZ'];
                    units = '[nT]';
                    ochid = {'Hx' 'Hy' 'Hz'};
                otherwise
                    error('This data type is not currently known in ModEM');
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [copyright] = newcopyright(title,authors,year,surveydoi)
            %  Usage : [copyright] = newcopyright(title,authors,year,surveydoi)
            %
            % a simply utility to help create new data copyright
            % information directly in Matlab
            
            copyright.Citation.Title = '';
            copyright.Citation.Authors = '';
            copyright.Citation.Year = '';
            copyright.Citation.SurveyDOI = 'doi:10.17611/DP/EMTF/MYPROJECT/MYSURVEY';
            
            if nargin >= 1; copyright.Citation.Title = title; end
            if nargin >= 2; copyright.Citation.Authors = authors; end
            if nargin >= 3; copyright.Citation.Year = year; end
            if nargin >= 4; copyright.Citation.SurveyDOI = ['doi:10.17611/DP/EMTF/' surveydoi]; end
            
            copyright.SelectedPublications = '';
            copyright.Acknowledgement = '';
            copyright.ConditionsOfUse = '';
            copyright.AdditionalInfo = '';
            
        end
        
        % **********************************************************************
        % Generates a general 2x2 rotation matrix that rotates any two horizontal channels from their original
        % layout (theta1, theta2) to a right hand orthogonal coordinate system with azimuth theta0.
        % Note that if and only if the original layout is also an orthogonal coordinate system, theta2 = theta1 + 90,
        % the determinant is 1 and this reduces to the more traditional form.
        % To reverse, use rot2inv.
        % 
        % All rotation codes are translated from EMTF-FCU
        % (c) A. Kelbert, 2017
        
        function R = rot2( theta1, theta2, theta0 )
            
            R(1,1) = cosd(theta1 - theta0);
            R(1,2) = cosd(theta2 - theta0);
            R(2,1) = sind(theta1 - theta0);
            R(2,2) = sind(theta2 - theta0);
            
        end
        
        % **********************************************************************
        % Generates the inverse of a general 2x2 rotation matrix that rotates any two horizontal channels
        % back to their original layout (theta1, theta2) from a right hand orthogonal coordinate system with azimuth theta0.
        % Same as matinv2( rot2(theta1, theta2, theta0) ).
        % 
        % All rotation codes are translated from EMTF-FCU
        % (c) A. Kelbert, 2017
        
        function Rinv = rot2inv( theta1, theta2, theta0 )
            
            det = sind(theta2 - theta1);
            
            Rinv(1,1) = sind(theta2 - theta0)/det;
            Rinv(1,2) = - cosd(theta2 - theta0)/det;
            Rinv(2,1) = - sind(theta1 - theta0)/det;
            Rinv(2,2) = cosd(theta1 - theta0)/det;
            
        end
        
        % **********************************************************************
        % Sets up the rotation matrices for input and output channels, respectively,
        % for a general EM data type. To be used to rotate from (and to) an original
        % site layout, which doesn't have to be orthogonal. Use this with EDI files
        % when ROT=NONE and with general Z-files.
        % Implicit assumptions: if a single output channel is found, do not rotate.
        % For an odd number of output channels, assume that Hz comes first.
        % Otherwise, rotate pairwise. Azimuth defaults to zero (geographic coordinates).
        % If azimuth is specified, rotates to a right hand orthogonal coordinate system
        % defined by that azimuth. If fwdORinv = 'INV', rotates back from the orthogonal
        % coordinate system to the original site layout.
        % Assumes 2 input channels, Hx and Hy, always (they don't need to be orthogonal,
        % but they are horizontal).
        % NOTE: The channels themselves are NOT rotated. We keep the original orientations
        %       for reference, we might want to go back to these later.
        %       If we want them rotated (e.g., for Z-files), there's a separate routine.
        % 
        % All rotation codes are translated from EMTF-FCU
        % (c) A. Kelbert, 2017

        function [U,V] = setup_rotation(orientation, old_azimuth, new_azimuth, fwdORinv)
            
            if (nargin <= 3)
                fwdORinv = 'FWD';
            end
            
            theta0 = new_azimuth;
            
            % create the rotation matrix for the two input channels first;
            % assuming Hx and Hy are the input channels, - pretty much always the case
            itheta(1) = orientation(1) + old_azimuth; ithetanew(1) = new_azimuth;
            itheta(2) = orientation(2) + old_azimuth; ithetanew(2) = new_azimuth + 90.0;
            if strcmpi(fwdORinv,'FWD')
                U = transpose(mttf.rot2inv(itheta(1),itheta(2),theta0));
            else
                U = transpose(mttf.rot2(itheta(1),itheta(2),theta0));
            end
            
            % now do output channels;
            % we rotate output channels two at a time, and we implicitly ignore the vertical magnetic field
            nchout = length(orientation) - 2;
            otheta = zeros(1,nchout);
            othetanew = zeros(1,nchout);
            V = zeros(nchout,nchout);
            if mod(nchout,2) == 1
                V(1,1) = 1.0;
                otheta(1) = 0.0;
                othetanew(1) = 0.0;
                i1 = 2;
            else
                i1 = 1;
            end
            for i=i1:2:nchout
                otheta(i) = orientation(i+2) + old_azimuth; othetanew(i) = new_azimuth;
                otheta(i+1) = orientation(i+3) + old_azimuth; othetanew(i+1) = new_azimuth + 90.0;
                if strcmpi(fwdORinv,'FWD')
                    V(i:i+1,i:i+1) = mttf.rot2(otheta(i),otheta(i+1),theta0);
                else
                    V(i:i+1,i:i+1) = mttf.rot2inv(otheta(i),otheta(i+1),theta0);
                end
            end
            
        end
        
    end
end
