classdef globalemdata
    %   read, write and redefine errors on global EM data
    %   this currently works with H (frequency domain magnetic fields)
    %   and C & D responses
    %   (c) Anna Kelbert, April 2014
    %
    % Gary's comment on HDZ array format:
    % these are three   components  of the complex magnetic field (in geomagnetic
    % coordinates HDZ(:,1,:) =    North ; HDZ(:,2,:) = East; HDZ(:,3,:) =  Down)
    % relative to a reference P10 signal computed  from fitting mid-latitude  observatories in
    % the time domain (i.e., this is one of the suggestions of Schmucker in his letter to you).
    % HDZsd gives corresponding standard errors, from the univariate TF estimation. and gives
    % at least some idea of quality for each component.   By rights, fitting of any parameters for
    % auroral correction should account for this (some sites are noisy or have minimal   data)
    % [nlat,ncomp,nper] = size(HDZ);

    properties
        header      % file header
        periods     % in days
        type        % H, C or D
        d           % data structure (array of periods)
        lat
        lon
        gmlat
        gmlon
        codes
        signConvention = 1  % 1 or -1
        uncorrected = 1  % true or false
        predicted = 0  % true or false
        rmsTotal
        rmsInfo
        rms
    end
    
    properties (SetAccess = protected)
        % This contains several different blocks. Each is a different way
        % of saving the same information; this is needed to work with
        % various historic scripts in a straightforward way without
        % additional code modifications...
        % BLOCK #1: Gary's data processing software
        % [nobs,ncomp,nper] (North,East,Down)
        HDZ
        HDZsd
        % BLOCK #2: Anna's source correction etc software
        % [nobs,nper] (East,South,Down)
        Hp
        Ht
        Hr
        Hpe
        Hte
        Hre
    end
    
    methods
        function [obj] = globalemdata(varargin)
            %   class constructor
            %   [obj] =
            %   globalemdata(cfile,dataType,predicted)
            
            if nargin == 0
                return
            end
            
            if nargin>0
                data = varargin{1};
            end
            if nargin>1
                obj.type = varargin{2};
            end
            if nargin>2
                obj.predicted = varargin{3};
            end

            if ischar(data)
                cfile = data;
                obj = read(cfile,obj.type,obj.predicted);
                return;
            elseif isstruct(data)
                obj.periods = zeros(1,length(data));
                for iper=1:length(data)
                    obj.periods(iper) = data(iper).period(1);
                end
                obj.codes = data(1).obs;
                obj.gmlat = data(1).lat;
                obj.gmlon = data(1).lon;
                for iper = 1:length(data)-1
                    [newcodes,ind]=setdiff(data(iper+1).obs,obj.codes);
                    obj.codes = [obj.codes; newcodes];
                    if ~isempty(ind)
                        obj.gmlat = [obj.gmlat; data(iper+1).lat(ind)];
                        obj.gmlon = [obj.gmlon; data(iper+1).lon(ind)];
                    end
                end
                [obj.lat,obj.lon]=gg2gm(90-obj.gmlat,obj.gmlon,-1);
                obj.lat=90-obj.lat; % Lat = 90 - Colat
            end
                        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = sortHDZ(obj,sortParameter)
            % sort by GM latitude or by observatory code
            % HDZ array is [nobs,ncomp,nper]
            % NOTE: data structure not sorted!
            
            if nargin < 2
                sortParameter = 'gmlat';
            end

            switch lower(sortParameter)
                case 'gmlat'
                    % sort all data by GM latitude
                    [~,ind] = sort(obj.gmlat);
                case 'codes'
                    % sort by alphabet
                    [~,ind] = sort(obj.codes);
                otherwise
                    error(['Don''t know how to sort by ' sortParameter])
            end
               
            obj.gmlat = obj.gmlat(ind);
            obj.gmlon = obj.gmlon(ind);
            obj.lat = obj.lat(ind);
            obj.lon = obj.lon(ind);
            obj.codes = obj.codes(ind);
            obj.HDZ = obj.HDZ(ind,:,:);
            obj.HDZsd = obj.HDZsd(ind,:,:);    
            
            % now sort the corresponding arrays (East, South and down)
            obj.Hp = squeeze(obj.HDZ(:,2,:));
            obj.Ht = - squeeze(obj.HDZ(:,1,:));
            obj.Hr = squeeze(obj.HDZ(:,3,:));
            obj.Hpe = squeeze(obj.HDZsd(:,2,:));
            obj.Hte = squeeze(obj.HDZsd(:,1,:));
            obj.Hre = squeeze(obj.HDZsd(:,3,:));
            
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function newobj = updateHDZ(obj,newHDZ,newHDZsd)
            % updates all data structures
            
            if nargin < 2
                error('Please supply new HDZ structure');
            end
            if isempty(newHDZ)
                error('New HDZ array cannot be empty');
            elseif min(size(newHDZ)== size(obj.HDZ))==0
                error('New HDZ array must have the same size as the old one');
            end
            
            newobj = obj;
            newobj.HDZ = newHDZ;
            if nargin > 2
                newobj.HDZsd = newHDZsd;
            end
            newobj = newobj.HDZ2data;
            
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function newobj = subset(obj,parameter,newrange)
            % this general subset function can at present be used to
            % subset by period range or GM latitude range, or codes
            % example usage:
            % newobj = obj.subset('period',[5 107])
            % newobj = obj.subset('codes',{'BEL','FUR','KAK'})
            % newobj = obj.subset('gmlat',[-60 60])
            
            newobj = obj;
            
            switch lower(parameter)
                case 'period'
                    ind = find(obj.periods>=newrange(1) & obj.periods<=newrange(2));
                    % subset the data
                    newobj.d = obj.d(ind);
                    newobj.periods = obj.periods(ind);
                    newobj = newobj.data2HDZ;
                case 'gmlat'
                    if isscalar(newrange)
                        newrange = [-newrange newrange];
                    end
                    ind = find(obj.gmlat>=newrange(1) & obj.gmlat<=newrange(2));
                    % subset the data
                    newobj.codes = obj.codes(ind);
                    newobj.gmlat = obj.gmlat(ind);
                    newobj.gmlon = obj.gmlon(ind);
                    newobj.lat = obj.lat(ind);
                    newobj.lon = obj.lon(ind);
                    for iper = 1:length(obj.d)
                        [~,ind]=intersect(obj.d(iper).obs,newobj.codes);
                        newobj.d(iper).period = obj.d(per).period(ind);
                        newobj.d(iper).obs = obj.d(per).obs(ind);
                        newobj.d(iper).lon = obj.d(per).lon(ind);
                        newobj.d(iper).lat = obj.d(per).lat(ind);
                        if isfield(newobj.d(iper),'H')
                            newobj.d(iper).H = obj.d(per).H(ind,:);
                            newobj.d(iper).Herr = obj.d(per).Herr(ind,:);
                        end
                        if isfield(newobj.d(iper),'C')
                            newobj.d(iper).C = obj.d(per).C(ind,:);
                            newobj.d(iper).Cerr = obj.d(per).Cerr(ind,:);
                        end
                        if isfield(newobj.d(iper),'D')
                            newobj.d(iper).D = obj.d(per).D(ind,:);
                            newobj.d(iper).Derr = obj.d(per).Derr(ind,:);
                        end
                    end
                    newobj = newobj.data2HDZ;
                case 'codes'
                    ind = intersect(obj.codes,newrange);
                    % subset the data
                    newobj.codes = obj.codes(ind);
                    newobj.gmlat = obj.gmlat(ind);
                    newobj.gmlon = obj.gmlon(ind);
                    newobj.lat = obj.lat(ind);
                    newobj.lon = obj.lon(ind);
                    for iper = 1:length(obj.d)
                        [~,ind]=intersect(obj.d(iper).obs,newobj.codes);
                        newobj.d(iper).period = obj.d(per).period(ind);
                        newobj.d(iper).obs = obj.d(per).obs(ind);
                        newobj.d(iper).lon = obj.d(per).lon(ind);
                        newobj.d(iper).lat = obj.d(per).lat(ind);
                        if isfield(newobj.d(iper),'H')
                            newobj.d(iper).H = obj.d(per).H(ind,:);
                            newobj.d(iper).Herr = obj.d(per).Herr(ind,:);
                        end
                        if isfield(newobj.d(iper),'C')
                            newobj.d(iper).C = obj.d(per).C(ind,:);
                            newobj.d(iper).Cerr = obj.d(per).Cerr(ind,:);
                        end
                        if isfield(newobj.d(iper),'D')
                            newobj.d(iper).D = obj.d(per).D(ind,:);
                            newobj.d(iper).Derr = obj.d(per).Derr(ind,:);
                        end
                    end
                    newobj = newobj.data2HDZ;
                otherwise
                    error(['Don''t know how to subset globalemdata by ' parameter])
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function newobj = sort(obj,sortParameter,sortDirection)
            % this general sort function can only at present be used to
            % sort by period
            
            if nargin < 2
                sortParameter = 'period';
                sortDirection = 'ascending';
            end
            
            newobj = obj;
            
            switch lower(sortParameter)
                case 'period'
                    if lower(sortDirection(1))=='a'
                        [newobj.periods,ind] = sort(obj.periods);
                    elseif lower(sortDirection(1))=='d'
                        [newobj.periods,ind] = sort(obj.periods,2,'descend');
                    else
                        error('Sort must be ascending or descending');
                    end
                    % sort the data structure
                    newobj.d = obj.d(ind);
                    % sort BLOCK #1
                    if ~isempty(obj.HDZ) && ~isempty(obj.HDZsd)
                        newobj.HDZ = obj.HDZ(:,:,ind);
                        newobj.HDZsd = obj.HDZsd(:,:,ind);
                    end
                    % sort BLOCK #2
                    if ~isempty(obj.Hp)
                        newobj.Hp = obj.Hp(:,ind);
                    end
                    if ~isempty(obj.Ht)
                        newobj.Ht = obj.Ht(:,ind);
                    end
                    if ~isempty(obj.Hr)
                        newobj.Hr = obj.Hr(:,ind);
                    end
                    if ~isempty(obj.Hpe)
                        newobj.Hpe = obj.Hpe(:,ind);
                    end
                    if ~isempty(obj.Hte)
                        newobj.Hte = obj.Hte(:,ind);
                    end
                    if ~isempty(obj.Hre)
                        newobj.Hre = obj.Hre(:,ind);
                    end
                otherwise
                    error(['Don''t know how to sort by ' sortParameter ' yet... sort not implemented'])
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = HDZ2data(obj)
            % recompute all other arrays after HDZ is updated
            
            if isempty(obj.HDZ)
                error('HDZ array does not exist in globalemdata object');
            end                
            
            % uncorrected data and errors (East, South and down)
            obj.Hp = squeeze(obj.HDZ(:,2,:));
            obj.Ht = - squeeze(obj.HDZ(:,1,:));
            obj.Hr = squeeze(obj.HDZ(:,3,:));
            obj.Hpe = squeeze(obj.HDZsd(:,2,:));
            obj.Hte = squeeze(obj.HDZsd(:,1,:));
            obj.Hre = squeeze(obj.HDZsd(:,3,:));
            
            % create data structure (starting from longest period to shortest)
            nper = length(obj.periods);
            nobs = length(obj.codes);
            if isempty(obj.d)
                obj.d(1:nper) = struct();
            end
            for iper = 1:nper
                for iobs = 1:nobs
                    obj.d(iper).period(iobs) = obj.periods(iper);
                    obj.d(iper).obs(iobs) = obj.codes(iobs);
                    obj.d(iper).lon(iobs) = obj.gmlon(iobs);
                    obj.d(iper).lat(iobs) = obj.gmlat(iobs);
                    obj.d(iper).H(iobs,1) = obj.Hp(iobs,iper);
                    obj.d(iper).H(iobs,2) = obj.Ht(iobs,iper);
                    obj.d(iper).H(iobs,3) = obj.Hr(iobs,iper);
                    obj.d(iper).Herr(iobs,1) = obj.Hpe(iobs,iper);
                    obj.d(iper).Herr(iobs,2) = obj.Hte(iobs,iper);
                    obj.d(iper).Herr(iobs,3) = obj.Hre(iobs,iper);
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = data2HDZ(obj)
            % recompute HDZ array from data structure
            % just to be safe, always sort by period first...
                        
            if isempty(obj.d)
                error('Data structure does not exist in globalemdata object');
            end
            
            obj = obj.sort('period','ascending');
            nper = length(obj.periods);
            nobs = length(obj.codes);
            obj.HDZ = nan(nobs,3,nper)+1i*nan(nobs,3,nper);
            obj.HDZsd = nan(nobs,3,nper);
            for iper = 1:nper
                [~,ista]=intersect(obj.codes,obj.d(iper).obs);
                if isfield(obj.d(iper),'H')
                    obj.HDZ(ista,1,iper) = - obj.d(iper).H(:,2); % Ht
                    obj.HDZ(ista,2,iper) =   obj.d(iper).H(:,1); % Hp
                    obj.HDZ(ista,3,iper) =   obj.d(iper).H(:,3); % Hr
                    obj.HDZsd(ista,1,iper) = obj.d(iper).Herr(:,2); % Ht
                    obj.HDZsd(ista,2,iper) = obj.d(iper).Herr(:,1); % Hp
                    obj.HDZsd(ista,3,iper) = obj.d(iper).Herr(:,3); % Hr
                end
            end
            
            % uncorrected data and errors (East, South and down)
            obj.Hp = squeeze(obj.HDZ(:,2,:));
            obj.Ht = - squeeze(obj.HDZ(:,1,:));
            obj.Hr = squeeze(obj.HDZ(:,3,:));
            obj.Hpe = squeeze(obj.HDZsd(:,2,:));
            obj.Hte = squeeze(obj.HDZsd(:,1,:));
            obj.Hre = squeeze(obj.HDZsd(:,3,:));

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = computeCD(obj)
            % use data structure to compute C & D responses

            if isempty(obj.d)
                obj = obj.HDZ2data;
            end
            
            R_e = 6371;
            % compute C and D responses (Ht is South, Hr is down)
            %Cratio = (R_e/2) * obj.Hr./obj.Ht;
            Cresp = (R_e/2) * diag(tan((90-obj.gmlat)*(pi/180))) * obj.Hr./obj.Ht;
            %Dratio = (R_e/2) * obj.Hp./obj.Ht;
            Dresp = (R_e/2) * diag(sin((90-obj.gmlat)*(pi/180))) * obj.Hp./obj.Ht;
            % compute C and D response error estimates
            c1 = 1./obj.Ht;
            c2 = obj.Hr./(obj.Ht.^2);
            Cerr = - (R_e/2) * diag(tan((90-obj.gmlat)*(pi/180))) * ...
                sqrt((abs(c1).*obj.Hre).^2 + (abs(c2).*obj.Hte).^2);
            d1 = 1./obj.Ht;
            d2 = obj.Hp./(obj.Ht.^2);
            Derr = - (R_e/2) * diag(sin((90-obj.gmlat)*(pi/180))) * ...
                sqrt((abs(d1).*obj.Hpe).^2 + (abs(d2).*obj.Hte).^2);
            
            for iper = 1:length(obj.d)
                [~,ista]=intersect(obj.codes,obj.d(iper).obs);
                obj.d(iper).C(ista) = Cresp(ista,iper);
                obj.d(iper).Cerr(ista) = abs(Cerr(ista,iper));
                obj.d(iper).D(ista) = Dresp(ista,iper);
                obj.d(iper).Derr(ista) = abs(Derr(ista,iper));
            end

        end
                        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function obj = addErrorFloor(obj,f)
%             % set the error floor to fraction f of the data;
%             % save the real errors in "dataErrors"
%             for i=1:obj.nPeriods
%                 if isempty(obj.d{i}.dataErrors)
%                     obj.d{i}.dataErrors = obj.d{i}.Zerr;
%                 end
%             end
%             obj.d = setErrorFloor3D(obj.d,f);
%                        
%         end
%  
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function obj = removeErrorFloor(obj)
%             % reverts data errors to original values
%             for i=1:obj.nPeriods
%                 if ~isempty(obj.d{i}.dataErrors)
%                     obj.d{i}.Zerr = obj.d{i}.dataErrors;
%                 end
%             end
%                        
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [h] = plot(obj,fig,lims,P)
            % Sample usage: h = obj.plot(gcf,lims,P)
            % Use this to plot the sites on any map or to create a new map
            % of the sites
            
            proj = 'miller';
            regional = 1;
            if nargin >= 4
                if isfield(P,'proj')
                    proj = P.proj;
                end
            end
            if strcmp(proj,'hammer')
                regional = 0;
            end
            
            if nargin < 3
                lims = modelplot.getLimits;
            end
            lonmin = lims.lonmin;
            lonmax = lims.lonmax;
            latmin = lims.latmin;
            latmax = lims.latmax;
                        
%             posn = [1,1,12,7];
%             set(fig,'Position',100*posn,...
%                 'PaperPosition',posn,...
%                 'PaperOrientation','Portrait',...
%                 'Color',[1 1 1]);            
%             m_proj(proj,...
%                 'long',[lonmin lonmax],...
%                 'lat', [latmin latmax]);
%             hold on;
%             cbndry  = [0.1 0.1 0.1];
%             m_coast('color',cbndry);              % Coastline...
%             m_coast('speckle','color',cbndry);    % with speckle added
%             hold on;
%             if regional
%                 m_grid('box','fancy','tickdir','in',... %'box','on',...
%                     'xtick',2,'linestyle','none',...
%                     'XaxisLocation','top',...
%                     'fontsize',28,'fontweight','demi');
%             else
%                 m_grid('box','on','tickdir','in',...
%                     'xtick',12,'ytick',[-45 -30 0 30 45],...
%                     'XaxisLocation','middle','xticklabels','',...
%                     'fontsize',28,'fontweight','demi');
%             end
            
            if nargin > 3
                obj.plot_sites(fig,obj,0,P);
            else
                obj.plot_sites(fig,obj);
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function write(obj,cfile,dataType)
            % 
            % computes the output data type as needed and writes to file
            
            if nargin < 3
                dataType = obj.type;
            end
            
            outobj = obj.sort('period','descending');
            switch upper(dataType)
                case 'H'
                    write_data(cfile,outobj.d,outobj.header,'H');
                case {'C','D'}
                    outobj = outobj.computeCD;
                    write_data(cfile,outobj.d,outobj.header,dataType);
                otherwise
                    error(['Can''t write data type ' dataType ' to file: method unknown']);
            end

        end
    end
    
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = read(cfile,dataType,predicted)
            %
            % Usage:  [obj] = read(cfile,dataType,predicted)
            %
            %  Reads an globalemdata object (backward compatible with old
            %  scripts)
            %  Note that in data files, periods are sorted from longest to
            %  shortest. This is done intentionally, so that when the 
            %  modeling is done in serial, the previous solution (with 
            %  hopefully a simpler and smoother structure) can be
            %  efficiently reused as a starting point for the following 
            %  period. On input, we always revert the order. We always
            %  want the periods to be sorted in ascending order for
            %  consistency with everything else.
            
            if nargin==0
                cfile = 'RTF210.mat';
                dataType = 'HDZ';
            elseif nargin == 1
                dataType = 'H';
            elseif ~exist(cfile,'file')
                error(['Data file ' cfile ' not found']);
            end
                        
            switch upper(dataType)
                case {'C','D'}
                    [data,info] = read_data(cfile,dataType);
                    obj = globalemdata(data,dataType);
                    obj.periods = fliplr(obj.periods);
                    obj.d = fliplr(data);
                    obj.header = info;
                case 'H'
                    [data,info] = read_hdata(cfile);
                    obj = globalemdata(data,dataType);
                    obj.periods = fliplr(obj.periods);
                    obj.d = fliplr(data);
                    obj.header = info;
                    obj = obj.data2HDZ;
                case 'HDZ'
                    load(cfile);
                    obj = globalemdata;
                    obj.HDZ = HDZ;
                    obj.HDZsd = HDZsd;
                    obj.lat = lat;
                    obj.lon = lon;
                    obj.gmlat = gmlat;
                    obj.gmlon = gmlon;
                    obj.codes = sta;
                    obj.type = 'H';
                    obj.header = '# Uncorrected magnetic fields from Egbert, 2012';
                    obj = obj.HDZ2data;                    
                otherwise
                    error(['Can''t read global data type ' dataType ': method unknown']);
            end
            
            if nargin > 2
                obj.predicted = predicted;
                if predicted
                    % in the future, will fix this! - shouldn't call an
                    % outside function and read the file again - this is a
                    % quick fix to get it to do what I want
                    [obj.rmsTotal, obj.rmsInfo] = read_rms(cfile,upper(dataType));
                    %% Average over observatory bins
                    if strcmp(upper(dataType),'H')
                        ncomp = 3;
                        comp = ['Hp';'Ht';'Hr'];
                    else
                        ncomp = 1;
                        comp = upper(dataType);
                    end
                    count = zeros(1,ncomp);
                    nper = length(obj.rmsInfo.per);
                    for k = 1:ncomp
                        for i = 1:length(obj.codes)
                            [~,ind] = intersect(obj.rmsInfo.obs,obj.codes(i));
                            count(k) = nper - sum(isnan(obj.rmsInfo.err(ind,:,k)));
                            obj.rms(i,k) = sqrt(nansum(obj.rmsInfo.err(ind,:,k))/count(k));
                        end
                    end
                end
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_sites(h,siteinfo,gm2gg,P)
            % On a given figure handle, plot some sites read from a *.coords
            % or a *.hout or *.cout file with data fits
            % This will only work for global and full impedance MT...
            % in general, for MT use mtdata and mtdataplot objects.
                        
            if ischar(siteinfo)
                % read the file first
                cfile = siteinfo;
                if findstr(cfile,'.hout') | findstr(cfile,'.cout')
                    if findstr(cfile,'.hout')
                        type = 'H';
                    elseif findstr(cfile,'.cout')
                        type = 'C';
                    end
                    [rms, info, obsrms] = read_rms(cfile,type,1); %#ok<ASGLU>
                    site.codes = obsrms.codes;
                    site.lat   = obsrms.lat;
                    site.lon   = obsrms.lon;
                    site.rms   = sqrt(sum(obsrms.rms.'.^2)/3);
                elseif findstr(cfile,'.coords')
                    % do not convert to geographic while reading -
                    % converting later in this function already
                    if nargin < 3
                        gm2gg = 0;
                    end
                    site = globalemdata.read_coords(cfile); 
                elseif findstr(cfile,'.dat') 
                    % assume MT - then compute site info externally through
                    % mtdatacompare.misfit - both real and predicted data
                    % needed in this case
                    %predobj = mtdata.read(cfile,'list','[mV/km]/[nT]','Full_Impedance');
                    %pdataobj = mtdataplot(predobj);
                    %site.lon = pdataobj.v.lon;
                    %site.lat = pdataobj.v.lat;
                    error('For MT data, use mtdatacompare.misfit to create sites to plot');
                else
                    error('Unknown file type for plotting site coordinates');
                end
            else
                site = siteinfo;
            end
            
            if nargin > 2
                if strcmpi(gm2gg,'gm2gg') || gm2gg==1
                    [site.lat,site.lon]=gg2gm(90 - site.lat,site.lon,'gm2gg');
                    site.lat=90-site.lat; % Lat = 90 - Colat
                end
            end
                    
            if ~isa(site,'globalemdata') && ~isfield(site,'codes')
                error('Please supply a structure of sites with codes, lat, lon & optional rms');
            end
            cmap = 'jet';
            cmap_dots = 'green';
            size_of_dots = 8;
            site_shape = 'o';
            if nargin >= 4
                if isfield(P,'cmap')
                    cmap = P.cmap;
                end
                if isfield(P,'sitecolor')
                    cmap_dots = P.sitecolor;
                end
                if isfield(P,'sitesize')
                    size_of_dots = P.sitesize;
                end
            end
            site_color = zeros(length(site.codes),3);
            site_round = zeros(length(site.codes),3);
            if (isa(site,'globalemdata') && ~isempty(site.rms)) || isfield(site,'rms')
                maxrms=3.0;
                if ~isempty(strfind(cmap_dots,'grey'))
                    map=cbrewer('seq','Greys',9);[mc,~]=size(map);
                elseif ~isempty(strfind(cmap_dots,'pink'))
                    map=cbrewer('seq','RdPu',9);[mc,~]=size(map);
                else
                    map=cbrewer('seq','YlGn',9);[mc,~]=size(map);
                end
                for k=1:length(site.rms)
                    ik=max(0,round(((site.rms(k)-1)/maxrms)*mc))+1;
                    ik=min(ik,mc);
                    site_color(k,:)=map(ik,:);
                    site_round(k,:) = [0. 0. 0.];
                end
            elseif ischar(cmap)
                for k=1:length(site.codes)
                    site_color(k,:) = [0.9 0.9 0.9];
                    site_round(k,:) = [0. 0. 0.];
                end
            else
                for k=1:length(site.codes)
                    site_color(k,:) = [0. 1. 0.]; %'g';
                    site_round(k,:) = [0. 0. 0.];
                end
            end
            if nargin >= 4
                if isfield(P,'site_color')
                    for k=1:length(site.codes)
                        site_color(k,:) = P.site_color;
                    end
                end
                if isfield(P,'site_round')
                    for k=1:length(site.codes)
                        site_round(k,:) = P.site_round;
                    end
                end
                if isfield(P,'site_shape')
                    site_shape = P.site_shape;
                end
                if isfield(P,'size_of_dots')
                    size_of_dots = P.size_of_dots;
                end
            end
            
            figure(h); hold on;
            for iobs=1:length(site.codes)
                m_plot(site.lon(iobs),site.lat(iobs),['w' site_shape],...
                    'markerfacecolor',site_color(iobs,:),...
                    'markeredgecolor',site_round(iobs,:),'markersize',size_of_dots);
                m_plot(site.lon(iobs)-360,site.lat(iobs),['w' site_shape],...
                    'markerfacecolor',site_color(iobs,:),...
                    'markeredgecolor',site_round(iobs,:),'markersize',size_of_dots);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function site = read_coords(cfile,gm2gg)
            % Reads an ascii file with colat,lon coordinates with two
            % header lines; used e.g., for global observatory *.coords
            % Optionally, converts from geomagnetic to geographic
            fid = fopen(cfile, 'r');
            textscan(fid,'%s',1,'delimiter','\n');
            textscan(fid,'%f',1,'delimiter','\n');
            temp=textscan(fid, '%4s %f %f');
            fclose(fid);
            site.codes=temp{1};
            site.lat=90-temp{2};
            site.lon=temp{3};
            if nargin > 1
                if strcmpi(gm2gg,'gm2gg') || gm2gg==1
                    gmlat = 90-site.lat;
                    gmlon = site.lon;
                    [site.lat,site.lon]=gg2gm(gmlat,gmlon,-1);
                    site.lat=90-site.lat; % Lat = 90 - Colat
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [codes,len] = codes_by_region(area)
            % Defines subsets of observatories based on a general geographic area
            % (without regard to political boundaries! - be careful in publication)
            % Regions: global
            %          EastChinaSea
            %          SeaOfJapan
            %          Japan
            %          China
            %          India
            %          MiddleEast
            %          Russia
            %          Africa
            %          SouthAmerica
            %          Alaska
            %          Canada
            %          NorthAmerica
            %          Australia
            %          Antarctica
            %          IndianOcean
            %          PacificOcean
            %          NorthernEurope
            %          Spain
            %          CentralEurope       
            
            if nargin < 1
                region = 'global';
            end
            
            region.EastChinaSea = {'YSS','VLA','MMB','ESA','MIZ','KAK','KNZ','HTY','SSO',};
            region.SeaOfJapan = {'BMT','SSH','LNP','KNY','CBI','GUA'};
            region.japan = {'YSS','MMB','ESA','MIZ','KAK','KNZ','HTY','SSO','KNY','CBI','GUA'};
            region.china = {'GLM','LZH','CDP','QIX','BMT','WHN','SSH','LNP','MUT','GZH','QGZ','PHU','THJ','SIL','SHL'};
            region.india = {'TRD','TIR','ETT','ANN','PND','VSK','HYB','NGP','ABG','UJJ','SAB','KRC'};
            region.middleeast = {'MLT','ELT','AMT','BGY','QSB','TFS','ASH','TKT','AAA','KSH'};
            region.russia = {'LNN','BOX','MOS','KZN','ARS','SVD','KGD','NVS','POD','IRT','MZL','YAK','MGD'};
            region.africa = {'HER','HBK','LMM','TAN','TSU','LUA','MFP','BNG','AAE','TAM','GUI','MBD','ASC'};
            region.southamerica = {'AIA','LIV','ARC','SGE','PST','TRW','PIL','VSS','LQA','HUA','SJG','PAB','KOU','TTB'};%'FUG'
            region.alaska = {'CWE','BRW','CMO','SHU','SIT'};
            region.canada = {'MBC','RES','ALE','THL','CBB','YKC','MEA','GLN','FCC','BLC','PBQ','IQA','GDH','NAQ'};
            region.northamerica = {'TEO','DLR','BSL','DAL','TUC','BOU','FRN','VIC','NEW','FRD','...','OTT','STJ'};
            region.australia = {'TUN','ING','TND','KDU','PMG','CTA','ASP','LRM','GNA','TOO','CNB','EYR'};
            region.antarctica = {'VNA','SNA','RBD','NVL','MOL','MAW','DVS','MIR','CSY','DRV','MCQ','SBA'};
            region.indianocean = {'CZT','PAF','AMS'};
            region.pacificocean = {'MID','HON','PPT','API'};
            region.northerneurope = {'LRV','VAL','HAD','ESK','LER','WIT','WNG','BFE','RSV','HLP','NUR','LOV','UPS','DOB',...
                'LYC','ABK','KIR','TRO','BJN','HRN','SOD','MMK','ARK'};
            region.spain = {'COI','SFS','SPT','TOL','ALM','EBR','LGR'};
            region.centraleurope = {'CLF','DOU','MAB','BFO','FUR','BDV','PRU','NGK','WIK','NCK','HRB','THY','AQU','PAG',...
                'PEG','IZN','ISK','...','GCK','ODE','KIV','MNK','LVV','BEL','SWI'};
            
            switch lower(area)
                case {'eastchinasea','seaofjapan','japan','china','india','middleeast','russia',...
                        'africa','southamerica','alaska','canada','northamerica',...
                        'australia','antarctica','indianocean','pacificocean',...
                        'northerneurope'}
                    codes = region.(lower(area));
                    len = length(codes);
                case 'global'
                    codes = region;
                    len = length(region.japan)+length(region.china)+length(region.india)+length(region.spain)+ ...
                        length(region.middleeast)+length(region.russia)+length(region.africa)+length(region.alaska)+ ...
                        length(region.southamerica)+length(region.northamerica)+length(region.canada)+ ...
                        length(region.australia)+length(region.antarctica)+length(region.pacificocean)+ ...
                        length(region.indianocean)+length(region.northerneurope)+length(region.centraleurope);
                otherwise
                    disp('Region not defined');
                    codes = [];
                    len = 0;            
            end            
        end
        
    end    
end