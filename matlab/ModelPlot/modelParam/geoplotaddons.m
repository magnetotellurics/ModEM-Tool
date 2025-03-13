classdef geoplotaddons
    %   functions to improve geographic plotting
    %   (a) Anna Kelbert, 2014-2016
    
    properties
        coords   = 'geographic'; % geographic or geomagnetic
    end

    methods
        
        function [h] = map(obj,lims,P)
            % Sample usage: h = obj.map(lims,P) or obj.map('USA')
            % Optional structure P contains plotting parameters.
            % Same as pcolor but used for plotting an empty map with data.
            % Full example:
            % [lims,posn] = latlontools.getLimits('SouthAmerica');
            % obj = llmodel;
            % obj.map(lims);
            % obj.plot_boundaries(gcf,'tectonic_boundaries_by_zone')
            % obj.plot_boundaries(gcf,'gshhg_global_country_borders.shp')
            % obj.plot_boundaries(gcf,'gshhg_usa_states_pol.shp')
            % obj.plot_boundaries(gcf,'physio.shp')
            % set(gcf,'PaperPosition',posn);        
            
            proj = 'miller';
            regional = 1;
            clong = 180;
            if nargin >= 3
                if isfield(P,'proj')
                    proj = P.proj;
                end
                if isfield(P,'clong')
                    clong = P.clong;
                end
            end
            
            if nargin == 1
                error('Please specify the lat/lon limits on input to map()');
            elseif ischar(lims)
                lims = latlontools.getLimits(lims);
            elseif isnumeric(lims)
                bounds = lims; clear lims;
                lims.lonmin = bounds(1);
                lims.lonmax = bounds(2);
                lims.latmin = bounds(3);
                lims.latmax = bounds(4);
            elseif isstruct(lims)
                % do nothing
            end
            
            lonmin = lims.lonmin;
            lonmax = lims.lonmax;
            latmin = lims.latmin;
            latmax = lims.latmax;
                
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
            m_proj(proj,...
                'long',[lonmin lonmax],...
                'lat', [latmin latmax]);
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
                    'xtick',4,'linestyle','none',...
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
            
        end
        
    end
    
    methods(Static)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_points(h,pts,shape,color,type,value)
            % On a given figure handle, plot some points read from a kml,
            % an xls file (e.g., volcanoes), or a shapefile
            % Optionally, subset by any structure field for plotting.
            % E.g., for volcanoes, "type" options are:  
            %                       Name
            %                       VolcanoType
            %                       RockType
            %                       TectonicSetting
            %                       Region
            %                       SubRegion           
            % Example use:
            % pts = read_volcanoes
            % modelplot.plot_points(h,pts,'^',[1 0 0],'RockType','Rhyolite');
            
            size_of_dots = 6;
            site_shape = '^';
            site_color = repmat([1 0 0],length(pts),1);
            site_round = repmat([0 0 0],length(pts),1);
            if nargin >= 3
                site_shape = shape;
            end
            if nargin >= 4
                if length(color) == length(pts)
                    site_color = color;
                elseif size(color,1) == 1
                    for ipt = 1:length(pts)
                        site_color(ipt,:) = color;
                    end
                else
                    error('The number of colors to plot does not match the number of points');
                end
            end
            ind = [];
            if nargin < 6
                type = 'Name';
                value = '';
            end
            for j=1:length(pts)
                if ~isfield(pts(j),type)
                    error(['No field ' type ' in the points structure']);
                elseif ~isempty(strfind(char(pts(j).(type)),value))
                    ind = [ind j];
                end
            end 
            
            figure(h);
            for ipt = 1:length(ind)
                ptslon = latlontools.lon360(pts(ind(ipt)).X);
                ptslat = pts(ind(ipt)).Y;
                m_plot(ptslon,ptslat,['w' site_shape],'linewidth',1,...
                    'markerfacecolor',site_color(ipt,:),...
                    'markeredgecolor',site_round(ipt,:),'markersize',size_of_dots);
                m_plot(ptslon-360,ptslat,['w' site_shape],'linewidth',1,...
                    'markerfacecolor',site_color(ipt,:),...
                    'markeredgecolor',site_round(ipt,:),'markersize',size_of_dots);
            end
                        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_boundaries(h,pts,width,color)
            % On a given figure handle, plot some boundaries read from a kml or a
            % shapefile
            %
            % Examples:
            % obj.plot_boundaries(gcf,'tectonic_boundaries_by_zone')
            % obj.plot_boundaries(gcf,'gshhg_global_country_borders.shp')
            % obj.plot_boundaries(gcf,'gshhg_usa_states_pol.shp')
            % obj.plot_boundaries(gcf,'physio.shp')
            
            if ischar(pts)
                % read the file first; assume Matlab format by default
                cfile = pts;
                if ~isempty(strfind(cfile,'.shp'))
                    pts = modelplot.read_shapefile(cfile);
                elseif ~isempty(strfind(cfile,'.kml'))
                    pts = modelplot.read_kml(cfile);
                else
                    load(cfile);
                    pts = placemark;
                end
            end
            if nargin < 3
                width = 2;
            end
            if nargin < 4
                cbndry = 0.5*ones(1,3);
            else
                cbndry = color;
            end
            
            plot_shapefile = 'n';
            plot_pts_by_zone = 'n';
            %purple = [0.5 0 0.5];
            magenta = [1 0 1];
            %red = [1 0 0];
            gray = [0.6 0.6 0.6];
            
            figure(h);
            if isfield(pts(1),'pts') % it's a structure of boundaries
                plot_pts_by_zone = 'y';
                if nargin < 4
                    % default colors for tectonic boundaries
                    cbndry = 0.5*ones(length(pts),3);
                    cbndry(1,:) = gray; % continental convergent boundary (purple) [1 0.75 0.8]
                    cbndry(2,:) = gray; % continental rift [0.25 0.88 0.82]
                    cbndry(3,:) = gray; % continental transform fault
                    cbndry(4,:) = gray; % oceanic convergent boundary (purple) [1 0.75 0.8]
                    cbndry(5,:) = gray; % oceanic spreading rift [0.25 0.88 0.82]
                    cbndry(6,:) = gray; % oceanic transform fault
                    cbndry(7,:) = magenta; % subduction zone (magenta)
                end
            elseif isfield(pts(1),'ncst')
                plot_shapefile = 'y';
            end
                        
            if plot_pts_by_zone=='y'
                for izone = 1:length(pts)
                    for i = 1:length(pts(izone).pts)
                        ptslon = latlontools.lon360(pts(izone).pts(i).X);
                        ptslat = pts(izone).pts(i).Y;
                        if ~ (max(abs(diff(ptslon))) > 180)
                            m_line(ptslon,ptslat,'linewidth',width,'color',cbndry(izone,:));
                            m_line(ptslon-360,ptslat,'linewidth',width,'color',cbndry(izone,:));
                        end
                    end
                end
            elseif plot_shapefile=='n'
                for i = 1:length(pts)
                    ptslon = latlontools.lon360(pts(i).X);
                    ptslat = pts(i).Y;
                    if ~ (max(abs(diff(ptslon))) > 180)
                        m_line(ptslon,ptslat,'linewidth',width,'color',cbndry);
                        m_line(ptslon-360,ptslat,'linewidth',width,'color',cbndry);
                    end
                end
            else
                for i = 1:length(pts.ncst)
                    ptslon = latlontools.lon360(pts.ncst{i}(:,1));
                    ptslat = pts.ncst{i}(:,2);
                    m_line(ptslon,ptslat,'linewi',width,'color',cbndry);
                    m_line(ptslon-360,ptslat,'linewi',width,'color',cbndry);
                end                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_gm_latitude(h,bounds)
            % On a given figure handle, shade values that exceed a certain
            % geomagnetic latitude
            
            if nargin < 1
                bounds = [-90 90];
            elseif isscalar(bounds)
                bounds = [-bounds bounds];
            end
            
            [lpatch_lat, lpatch_lon] = meshgrid(bounds(1), lon);
            [lpatch_lat, lpatch_lon] = gg2gm(90-lpatch_lat,lpatch_lon,-1);
            lpatch_lat=90-lpatch_lat; % Lat = 90 - Colat
            in=find(lpatch_lon>180+clong);lpatch_lon(in)=lpatch_lon(in)-360;
            [lpatch_lon, i] = sort(lpatch_lon);
            lpatch_lat = lpatch_lat(i);
            dlat = 1;
            lat_ext1 = -90:dlat:lpatch_lat(1); lat_ext2 = lpatch_lat(end):-dlat:-90;
            lpatch_lat = [lat_ext1 lpatch_lat' lat_ext2];
            lon_ext1 = ones(1,length(lat_ext1))*lpatch_lon(1); lon_ext2 = ones(1,length(lat_ext2))*lpatch_lon(end);
            lpatch_lon = [lon_ext1 lpatch_lon' lon_ext2];
            
            [upatch_lat, upatch_lon] = meshgrid(bounds(2), lon);
            [upatch_lat, upatch_lon] = gg2gm(90-upatch_lat,upatch_lon,-1);
            upatch_lat=90-upatch_lat; % Lat = 90 - Colat
            in=find(upatch_lon>180+clong);upatch_lon(in)=upatch_lon(in)-360;
            [upatch_lon, i] = sort(upatch_lon);
            upatch_lat = upatch_lat(i);
            dlat = 1;
            lat_ext1 = 90:-dlat:upatch_lat(1); lat_ext2 = upatch_lat(end):dlat:90;
            upatch_lat = [lat_ext1 upatch_lat' lat_ext2];
            lon_ext1 = ones(1,length(lat_ext1))*upatch_lon(1); lon_ext2 = ones(1,length(lat_ext2))*upatch_lon(end);
            upatch_lon = [lon_ext1 upatch_lon' lon_ext2];
            
            figure(h);
            m_patch(upatch_lon,upatch_lat,[0.7 0.7 0.7]);
            m_patch(lpatch_lon,lpatch_lat,[0.7 0.7 0.7]);

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function placemark = read_kml(kmlfile)
            % Reads a kml file with boundaries for plotting
            % e.g. pts = read_kml('EarthsTectonicPlates.kml');
            if ~isempty(strfind(kmlfile,'.kml'))
                fileroot = kmlfile(1:end-4);
            else
                fileroot = kmlfile;
            end
            if exist([fileroot '.mat'],'file')
                load([fileroot '.mat']);
                return
            end
            t = xmltree([fileroot '.kml']);
            placemarks = find(t,'/kml/Document/Folder/Placemark');
            placemark(1:length(placemarks)) = {};
            for i=1:length(placemarks)
                geom = find(t,placemarks(i),'name','MultiGeometry'); %#ok<*GTARG>
                name = get(t,children(t,find(t,placemarks(i),'name','name')),'value');
                if ~isempty(geom)
                    lines = find(t,geom,'name','LineString');
                    placemark(i).name = name;
                    for j=1:length(lines)
                        attr = attributes(t,'get',lines(j),1); id = attr.val;
                        coords = get(t,children(t,find(t,lines(j),'name','coordinates')),'value');
                        temp = strsplit(coords,' ');
                        placemark(i).pts(j).Geometry = 'LineString';
                        placemark(i).pts(j).Id = attr.val;
                        for k=1:length(temp)
                            loc = strsplit(char(temp(k)),',');
                            placemark(i).pts(j).X(k) = str2double(loc(1));
                            placemark(i).pts(j).Y(k) = str2double(loc(2));
                        end
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function placemark = read_shapefile(shpfile)
            % Reads a shp file with boundaries for plotting
            % e.g. pts = read_shapefile('physio.shp');
            if findstr(shpfile,'.shp')
                fileroot = shpfile(1:end-4);
            else
                fileroot = shpfile;
            end
            if exist([fileroot '.mat'],'file')
                load([fileroot '.mat']);
                return
            end
            placemark = m_shaperead(fileroot);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [placemark,ind] = read_volcanoes(type,value)
            % Reads the GVP_volcano_list.xlsx file for plotting
            % e.g. [pts,ind] = read_volcanoes('RockType','Rhyolite');
            %      "type" options:  Name
            %                       VolcanoType
            %                       RockType
            %                       TectonicSetting
            %                       Region
            %                       SubRegion
            [~,~,alldata] = xlsread('GVP_volcano_list.xlsx');
            %nonsense = alldata(1,:);
            %headers = alldata(2,:);
            gvp = alldata(3:end,:);
            placemark(1:length(gvp)) = struct();
            for j=1:length(gvp)  
                placemark(j).Geometry = 'Point';
                placemark(j).Id = cell2mat(gvp(j,1));
                placemark(j).Name = cellstr(gvp(j,2));                
                placemark(j).X = cell2mat(gvp(j,10)); % longitude
                placemark(j).Y = cell2mat(gvp(j,9)); % latitude
                placemark(j).Z = cell2mat(gvp(j,11)); % elevation
                placemark(j).VolcanoType = cellstr(gvp(j,4));                
                placemark(j).RockType = cellstr(gvp(j,12));                
                placemark(j).TectonicSetting = cellstr(gvp(j,13));                
                placemark(j).Region = cellstr(gvp(j,7));                
                placemark(j).Subregion = cellstr(gvp(j,8));
            end
            ind = [];
            if nargin < 2
                type = 'Name';
                value = '';
            end
            for j=1:length(placemark)
                if ~isempty(strfind(char(placemark(j).(type)),value))
                    ind = [ind j];
                end
            end                
        end
        end
    
end
