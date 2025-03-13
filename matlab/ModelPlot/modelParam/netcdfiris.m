classdef netcdfiris < handle
    %   This is intended as an extension of netcdf class to read and write
    %   IRIS Earth Model Collaboration files; Matlab doesn't seem to allow 
    %   extension of netcdf so make it an independent class.
    %   Read and write NetCDF files in the format longitude, latitude, depth
    %   and multiple variables in the format X(longitude, latitude, depth).
    %   In some cases, dimensions in NetCDF file are different - then, we 
    %   fix them upon reading so that the output conforms to the same standard.
    %   For writing to NetCDF files, always use this dimension order.
    %   If needed, run permute(X,[2 1 3]) before writing to file.
    %   --
    %   (c) Anna Kelbert, 10-25 March 2014. 
    %   College of Earth, Ocean and Atmospheric Sciences,
    %   Oregon State University, Corvallis, OR 97331-5503
    %   (c) Anna Kelbert, modified 28 Dec 2022 for more general input.
    %   (c) Anna Kelbert, modified 30 Jan 2023 to add model variables.
    %   This allows direct for I/O of the model file with ModEM code.
    
    properties
        lims
        latitude
        longitude
        depth
        vars
        myind
        header
        modelvar
        geospatial
        cfile
 	end

    methods
        function [obj] = netcdfiris(cfile)
            %   class constructor
            obj.cfile = cfile;
            ncid = netcdfiris.open(obj.cfile);
            obj.lims = netcdfiris.getLimits(ncid);
            [obj.header,obj.geospatial] = netcdfiris.getHeader(ncid);
            [obj.modelvar] = netcdfxyz.getModelVariables(ncid);
            [obj.longitude,obj.latitude,obj.depth] = netcdfiris.getPoints(ncid);
            [obj.vars,obj.myind] = netcdfiris.getValue(ncid);
            netcdfiris.close(ncid);
        end
        
        function write(obj,cfile)
            %ncid = netcdfiris.open(cfile);
            ncid = netcdf.create(cfile,'NC_NOCLOBBER');
            netcdf.endDef(ncid);
            netcdfiris.putHeader(ncid,obj.header,obj.geospatial);
            netcdfxyz.putModelVariables(ncid,obj.modelvar);
            netcdfiris.putPoints(ncid,obj.longitude,obj.latitude,obj.depth);
            netcdfiris.putValue(ncid,obj.vars);
            netcdfiris.close(ncid);
        end
        
    end
        
    methods(Static)
        
        function [ncid] = open(cfile)

            netcdf.setDefaultFormat('NC_FORMAT_NETCDF4_CLASSIC');

            if exist(cfile,'file')==0
                % open an empty file for writing
                %ncid = netcdf.open(cfile,'NC_WRITE');
                ncid = netcdf.create(cfile,'NC_NOCLOBBER');
                netcdf.endDef(ncid);
                
            else
                % open for reading only
                ncid = netcdf.open(cfile,'NC_NOWRITE');
            end
        end
        
        function close(ncid)
            
            netcdf.close(ncid);

        end
        
        function [header,geospatial] = getHeader(ncid)
            
            ncglobal = netcdf.getConstant('NC_GLOBAL');
            
            header.title = netcdf.getAtt(ncid,ncglobal,'title');
            header.id = netcdf.getAtt(ncid,ncglobal,'id');
            try
                header.data_revision = netcdf.getAtt(ncid,ncglobal,'data_revision');
            catch exception
                if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                    header.data_revision = '';
                end
            end
            header.summary = netcdf.getAtt(ncid,ncglobal,'summary');
            header.keywords = netcdf.getAtt(ncid,ncglobal,'keywords');
            header.Conventions = netcdf.getAtt(ncid,ncglobal,'Conventions');
            try
                header.Metadata_Conventions = netcdf.getAtt(ncid,ncglobal,'Metadata_Conventions');
            catch exception
                if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                    header.Metadata_Conventions = '';
                end
            end
            try
                header.author_name = netcdf.getAtt(ncid,ncglobal,'author_name');
            catch
                try
                    header.author_name = netcdf.getAtt(ncid,ncglobal,'creator_name');
                catch exception
                    if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                        header.author_name = '';
                    end
                end
            end            
            try
                header.author_url = netcdf.getAtt(ncid,ncglobal,'author_url');
            catch
                try
                    header.author_url = netcdf.getAtt(ncid,ncglobal,'creator_url');
                catch exception
                    if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                        header.author_url = '';
                    end
                end
            end
            try
                header.author_email = netcdf.getAtt(ncid,ncglobal,'author_email');
            catch
                try
                    header.author_email = netcdf.getAtt(ncid,ncglobal,'creator_email');
                catch exception
                    if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                        header.author_email = '';
                    end
                end
            end
            header.institution = netcdf.getAtt(ncid,ncglobal,'institution');
            try
                header.repository_name = netcdf.getAtt(ncid,ncglobal,'repository_name');
            catch
                header.repository_name = '';
            end
            try
                header.repository_institution = netcdf.getAtt(ncid,ncglobal,'repository_institution');
            catch
                header.repository_institution = '';
            end
            try
                header.repository_pid = netcdf.getAtt(ncid,ncglobal,'repository_pid');
            catch
                header.repository_pid = '';
            end
            try
                header.acknowledgement = netcdf.getAtt(ncid,ncglobal,'acknowledgment');
            catch
                header.acknowledgement = '';
            end
            try
                header.references = netcdf.getAtt(ncid,ncglobal,'references');
            catch
                try
                    header.references = netcdf.getAtt(ncid,ncglobal,'reference');
                catch exception
                    if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                        header.references = '';
                    end
                end
            end
            try
                header.history = netcdf.getAtt(ncid,ncglobal,'history');
            catch
                header.history = '';
            end
            try
                header.comment = netcdf.getAtt(ncid,ncglobal,'comment');
            catch
                header.comment = '';
            end
            
            % Read geospatial information if present, otherwise skip
            geospatial.exists = 1;
            try
                netcdf.getAtt(ncid,ncglobal,'geospatial_lat_min');
            catch
                geospatial.exists = 0;
            end
            if geospatial.exists
                geospatial.lat_min = netcdf.getAtt(ncid,ncglobal,'geospatial_lat_min');
                geospatial.lat_max = netcdf.getAtt(ncid,ncglobal,'geospatial_lat_max');
                geospatial.lat_units = netcdf.getAtt(ncid,ncglobal,'geospatial_lat_units');
                geospatial.lat_resolution = netcdf.getAtt(ncid,ncglobal,'geospatial_lat_resolution');
                geospatial.lon_min = netcdf.getAtt(ncid,ncglobal,'geospatial_lon_min');
                geospatial.lon_max = netcdf.getAtt(ncid,ncglobal,'geospatial_lon_max');
                geospatial.lon_units = netcdf.getAtt(ncid,ncglobal,'geospatial_lon_units');
                geospatial.lon_resolution = netcdf.getAtt(ncid,ncglobal,'geospatial_lon_resolution');
                geospatial.vertical_min = netcdf.getAtt(ncid,ncglobal,'geospatial_vertical_min');
                geospatial.vertical_max = netcdf.getAtt(ncid,ncglobal,'geospatial_vertical_max');
                geospatial.vertical_units = netcdf.getAtt(ncid,ncglobal,'geospatial_vertical_units');
                geospatial.vertical_positive = netcdf.getAtt(ncid,ncglobal,'geospatial_vertical_positive');
            end
                        
        end

        function [modelVariables] = getModelVariables(ncid)
            
            ncglobal = netcdf.getConstant('NC_GLOBAL');
            modelVariables.primary_coords = netcdf.getAtt(ncid,ncglobal,'model_primary_coords');
            modelVariables.corner_description = netcdf.getAtt(ncid,ncglobal,'model_corner_description');
            modelVariables.corner_location = netcdf.getAtt(ncid,ncglobal,'model_corner_location');
            modelVariables.corner_latitude = netcdf.getAtt(ncid,ncglobal,'model_corner_latitude');
            modelVariables.corner_longitude = netcdf.getAtt(ncid,ncglobal,'model_corner_longitude');
            modelVariables.corner_depth = netcdf.getAtt(ncid,ncglobal,'model_corner_depth');
            modelVariables.rotation_units = netcdf.getAtt(ncid,ncglobal,'model_rotation_units');
            modelVariables.rotation_angle = netcdf.getAtt(ncid,ncglobal,'model_rotation_angle');
                        
        end
        
        function [lims] = getLimits(ncid)
            
            ncglobal = netcdf.getConstant('NC_GLOBAL');
            lims.lonmin = str2double(netcdf.getAtt(ncid,ncglobal,'geospatial_lon_min'));
            lims.lonmax = str2double(netcdf.getAtt(ncid,ncglobal,'geospatial_lon_max'));
            lims.latmin = str2double(netcdf.getAtt(ncid,ncglobal,'geospatial_lat_min'));
            lims.latmax = str2double(netcdf.getAtt(ncid,ncglobal,'geospatial_lat_max'));
            lims.depthmin = str2double(netcdf.getAtt(ncid,ncglobal,'geospatial_vertical_min'));
            lims.depthmax = str2double(netcdf.getAtt(ncid,ncglobal,'geospatial_vertical_max'));
            
        end
        
        function [longitude, latitude, depth] = getPoints(ncid)
            
            longitude = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'longitude'));
            longitude = double(longitude);
            
            latitude = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latitude'));
            latitude = double(latitude);
            
            depth = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'depth'));
            depth = double(depth);
                        
        end
        
        function [vars,myind] = getValue(ncid,nout)
            
            [ndims,nvars] = netcdf.inq(ncid);

            %allow for the order of variables to be mixed up with dims
            nvars = nvars - ndims;
            if nargin < 2
                nout = nvars;
            end
            if nvars >= 1
                vars(1:nvars) = struct;
                netcdfind(1:nvars) = 0;
                i = 1;
                for j=1:nvars+ndims
                    % read variables and dimensions; exit if done - this
                    % is only really needed for some invalid netcdf files
                    try
                       readme = netcdf.inqVar(ncid,j);
                    catch exception
                        if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                            break
                        end
                    end
                    % now skip if this is a dimension, otherwise read on
                    if strcmp(readme,'longitude') || strcmp(readme,'latitude') || strcmp(readme,'depth')
                        continue
                    end
                    vars(i).short_name = netcdf.inqVar(ncid,j);
                    vars(i).long_name = netcdf.getAtt(ncid,j,'long_name');
                    if strcmp(netcdf.inqAttName(ncid,j,1),'display_name')
                        vars(i).display_name = netcdf.getAtt(ncid,j,'display_name');
                    else
                        ind = strfind(vars(i).long_name,'in')-1;
                        if isempty(ind)
                            ind = length(vars(i).long_name);
                        end
                        vars(i).display_name = vars(i).long_name(1:ind);
                    end
                    try
                        vars(i).units = netcdf.getAtt(ncid,j,'units');
                    catch
                        ind = strfind(vars(i).display_name,'in')+3;
                        if isempty(ind)
                            vars(i).units = '';
                        else
                            vars(i).units = vars(i).long_name(ind:end);
                        end
                    end
                    try
                        vars(i).missing_value = netcdf.getAtt(ncid,j,'missing_value');
                    catch exception
                        if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                            vars(i).missing_value = NaN;
                        end
                    end
                    vars(i).value = netcdf.getVar(ncid,j);
                    nanind = vars(i).value==vars(i).missing_value;
                    vars(i).value(nanind) = NaN;
                    vars(i).value = double(vars(i).value);
                    disp('Output converted to double for plotting in Matlab');
                    
                    netcdfind(i) = j;
                    i = i+1;
                end
            end
            
            % choose the index for output
            myind = 1;
            if nvars > 1
                disp(sprintf('Reading the following NetCDF variables:')) %#ok<*DSPS>
                for i=1:nvars
                    disp(sprintf('%d: %s\t%s\n',i,vars(i).short_name,vars(i).long_name))
                end
                if nargin > 1
                    if nout == 1
                        reply = input('Please make an integer choice [1]: ', 's');
                        if ~isempty(reply); myind = str2double(reply); end
                        disp(['Reading the variable ' vars(myind).short_name ' from NetCDF file']);
                    end
                end
            else
                disp(['Reading the variable ' vars(myind).short_name ' from NetCDF file']);
            end
            
            % check dimensions based on the first variable
            [~,~,dimids] = netcdf.inqVar(ncid,netcdfind(1));
            if length(dimids)<3
                % reading a 1D or 2D variable, likely because the order of
                % variables is mixed up. This should not happen, throw a warning
                warning('The variables are less than 3D - error reading NetCDF');
                return
            end
            disp(['1st dimension: ',netcdf.inqDim(ncid,dimids(1))]);
            disp(['2nd dimension: ',netcdf.inqDim(ncid,dimids(2))]);
            disp(['3rd dimension: ',netcdf.inqDim(ncid,dimids(3))]);
            [~,~,londimid] = netcdf.inqVar(ncid,netcdf.inqVarID(ncid,'longitude'));
            [~,~,latdimid] = netcdf.inqVar(ncid,netcdf.inqVarID(ncid,'latitude'));
            [~,~,depdimid] = netcdf.inqVar(ncid,netcdf.inqVarID(ncid,'depth'));
            newdimids = [londimid latdimid depdimid];
            if sum(dimids == newdimids) == 3 
                % order is X(longitude,latitude,depth) - do nothing
            else
                % compute the permutation matrix and rearrange the output
                disp('Output permuted to conform to X(longitude,latitude,depth)');
                [~,ii] = sort(dimids);
                [~,jj] = sort(newdimids);
                P = zeros(3,3);
                P(jj+(ii-1)*3) = 1;
                newdims = (P*[1 2 3]')';
                for i=1:nvars
                    vars(i).value = permute(vars(i).value,newdims);
                end
            end
                        
        end
       
        % accepts variables in the form X(longitude,latitude,depth)
        % if needed, call permute(X,[2 1 3]) before writing
        function putValue(ncid,vars,myind)
                        
            lonid = netcdf.inqVarID(ncid,'longitude');
            [~,~,londimid] = netcdf.inqVar(ncid,lonid);
            latid = netcdf.inqVarID(ncid,'latitude');
            [~,~,latdimid] = netcdf.inqVar(ncid,latid);
            depid = netcdf.inqVarID(ncid,'depth');
            [~,~,depdimid] = netcdf.inqVar(ncid,depid);
            
            % Verify variable size - not a full proof method (especially
            % for a regular cube) but checks for obvious errors nevertheless
            [~,nlon] = netcdf.inqDim(ncid,londimid);
            [~,nlat] = netcdf.inqDim(ncid,latdimid);
            [~,ndep] = netcdf.inqDim(ncid,depdimid);
            if sum(size(vars(1).value) == [nlon nlat ndep]) ~= 3
                disp('Error writing the variable(s) to NetCDF: dimensions are inconsistent');
                disp('Correct format is: X(longitude,latitude,depth)');
                disp('If needed, call permute(X,[2 1 3]) before writing to NetCDF');
            end

            % Put file in define mode.
            netcdf.reDef(ncid);
            varid = zeros(1,length(vars));
            for i = 1:length(vars)
                if nargin > 2
                    if i ~= myind; continue; end
                end
                %%%
                varname = matlab.lang.makeValidName(vars(i).short_name,'ReplacementStyle','delete');
                varid(i) = netcdf.defVar(ncid,varname,'float',[londimid latdimid depdimid]);
            end            
            % Take file out of define mode.
            netcdf.endDef(ncid);

            for i = 1:length(vars)
                if nargin > 2
                    if i ~= myind; continue; end
                end
                disp('Input converted to single to write as IRIS NetCDF file');
                myvalue = single(vars(i).value);               
                missing_value = 99999.0;
                if isfield(vars(i),'missing_value')
                    missing_value = vars(i).missing_value;
                end
                nanind = isnan(myvalue); 
                myvalue(nanind) = missing_value;
                netcdf.putVar(ncid,varid(i),myvalue);
                netcdf.reDef(ncid);
                netcdf.putAtt(ncid,varid(i),'long_name',vars(i).long_name);
                netcdf.putAtt(ncid,varid(i),'display_name',vars(i).display_name);
                netcdf.putAtt(ncid,varid(i),'units',vars(i).units);
                netcdf.putAtt(ncid,varid(i),'missing_value',missing_value);
                % NetCDF4 doesn't seem to allow the _FillValue attribute
                %netcdf.putAtt(ncid,varid(i),'_FillValue',missing_value);
                netcdf.endDef(ncid);
            end            

            % Synchronize the file to disk.
            netcdf.sync(ncid)

        end
  
        function putPoints(ncid, longitude, latitude, depth)
            
            % Put file in define mode.
            netcdf.reDef(ncid);
            
            londimid = netcdf.defDim(ncid, 'longitude', length(longitude));
            latdimid = netcdf.defDim(ncid, 'latitude', length(latitude));
            depdimid = netcdf.defDim(ncid, 'depth', length(depth));

            lonid = netcdf.defVar(ncid,'longitude','float',londimid);
            latid = netcdf.defVar(ncid,'latitude','float',latdimid);            
            depid = netcdf.defVar(ncid,'depth','float',depdimid);
            
            % Take file out of define mode.
            netcdf.endDef(ncid);

            netcdf.putVar(ncid,lonid,single(longitude));
            netcdf.reDef(ncid);
            netcdf.putAtt(ncid,lonid,'long_name','Longitude; positive east');
            netcdf.putAtt(ncid,lonid,'units','degrees_east');
            netcdf.putAtt(ncid,lonid,'standard_name','longitude');
            netcdf.endDef(ncid);

            netcdf.putVar(ncid,latid,single(latitude));
            netcdf.reDef(ncid);
            netcdf.putAtt(ncid,latid,'long_name','Latitude; positive north');
            netcdf.putAtt(ncid,latid,'units','degrees_north');
            netcdf.putAtt(ncid,latid,'standard_name','latitude');
            netcdf.endDef(ncid);

            netcdf.putVar(ncid,depid,single(depth));
            netcdf.reDef(ncid);
            netcdf.putAtt(ncid,depid,'long_name','depth below earth surface');
            netcdf.putAtt(ncid,depid,'units','km');
            netcdf.putAtt(ncid,depid,'positive','down');
            netcdf.endDef(ncid);

            % Synchronize the file to disk.
            netcdf.sync(ncid)

        end
        
        function putModelVariables(ncid,modelVariables)
            
            ncglobal = netcdf.getConstant('NC_GLOBAL');
           
            % Put file in define mode.
            netcdf.reDef(ncid);
            netcdf.putAtt(ncid,ncglobal,'model_primary_coords',modelVariables.primary_coords);
            netcdf.putAtt(ncid,ncglobal,'model_corner_description',modelVariables.corner_description);
            netcdf.putAtt(ncid,ncglobal,'model_corner_location',modelVariables.corner_location);
            netcdf.putAtt(ncid,ncglobal,'model_corner_latitude',modelVariables.corner_latitude);
            netcdf.putAtt(ncid,ncglobal,'model_corner_longitude',modelVariables.corner_longitude);
            netcdf.putAtt(ncid,ncglobal,'model_corner_depth',modelVariables.corner_depth);
            netcdf.putAtt(ncid,ncglobal,'model_rotation_units',modelVariables.rotation_units);
            netcdf.putAtt(ncid,ncglobal,'model_rotation_angle',modelVariables.rotation_angle);
            % Take file out of define mode.
            netcdf.endDef(ncid);
            
            % Synchronize the file to disk.
            netcdf.sync(ncid)
                        
        end
        
        function putHeader(ncid,header,geospatial,limits)
            
            ncglobal = netcdf.getConstant('NC_GLOBAL');
            
            % Put file in define mode.
            netcdf.reDef(ncid);
            netcdf.putAtt(ncid,ncglobal,'title',header.title);
            netcdf.putAtt(ncid,ncglobal,'id',header.id);
            netcdf.putAtt(ncid,ncglobal,'data_revision',header.data_revision);
            netcdf.putAtt(ncid,ncglobal,'summary',header.summary);
            netcdf.putAtt(ncid,ncglobal,'keywords',header.keywords);
            netcdf.putAtt(ncid,ncglobal,'Conventions',header.Conventions);
            netcdf.putAtt(ncid,ncglobal,'Metadata_Conventions',header.Metadata_Conventions);
            netcdf.putAtt(ncid,ncglobal,'author_name',header.author_name);
            netcdf.putAtt(ncid,ncglobal,'author_url',header.author_url);
            netcdf.putAtt(ncid,ncglobal,'author_email',header.author_email);
            netcdf.putAtt(ncid,ncglobal,'institution',header.institution);
            netcdf.putAtt(ncid,ncglobal,'repository_name',header.repository_name);
            netcdf.putAtt(ncid,ncglobal,'repository_institution',header.repository_institution);
            netcdf.putAtt(ncid,ncglobal,'repository_pid',header.repository_pid);
            netcdf.putAtt(ncid,ncglobal,'acknowledgment',header.acknowledgement);
            netcdf.putAtt(ncid,ncglobal,'references',header.references);
            netcdf.putAtt(ncid,ncglobal,'history',header.history);
            netcdf.putAtt(ncid,ncglobal,'comment',header.comment);
            % Take file out of define mode.
            netcdf.endDef(ncid);
            
            % Put file in define mode.
            netcdf.reDef(ncid);
            if nargin < 3
                netcdf.putAtt(ncid,ncglobal,'geospatial_lat_min',sprintf('%6.1f',limits.latmin));
                netcdf.putAtt(ncid,ncglobal,'geospatial_lat_max',sprintf('%6.1f',limits.latmax));
            else
                netcdf.putAtt(ncid,ncglobal,'geospatial_lat_min',geospatial.lat_min);
                netcdf.putAtt(ncid,ncglobal,'geospatial_lat_max',geospatial.lat_max);
            end
            netcdf.putAtt(ncid,ncglobal,'geospatial_lat_units',geospatial.lat_units);
            netcdf.putAtt(ncid,ncglobal,'geospatial_lat_resolution',geospatial.lat_resolution);
            if nargin < 3
                netcdf.putAtt(ncid,ncglobal,'geospatial_lon_min',sprintf('%6.1f',limits.lonmin));
                netcdf.putAtt(ncid,ncglobal,'geospatial_lon_max',sprintf('%6.1f',limits.lonmax));
            else
                netcdf.putAtt(ncid,ncglobal,'geospatial_lon_min',geospatial.lon_min);
                netcdf.putAtt(ncid,ncglobal,'geospatial_lon_max',geospatial.lon_max);
            end
            netcdf.putAtt(ncid,ncglobal,'geospatial_lon_units',geospatial.lon_units);
            netcdf.putAtt(ncid,ncglobal,'geospatial_lon_resolution',geospatial.lon_resolution);
            if nargin < 3
                netcdf.putAtt(ncid,ncglobal,'geospatial_vertical_min',sprintf('%4.0f',limits.depthmin));
                netcdf.putAtt(ncid,ncglobal,'geospatial_vertical_max',sprintf('%4.0f',limits.depthmax));
            else
                netcdf.putAtt(ncid,ncglobal,'geospatial_vertical_min',geospatial.vertical_min);
                netcdf.putAtt(ncid,ncglobal,'geospatial_vertical_max',geospatial.vertical_max);
            end
            netcdf.putAtt(ncid,ncglobal,'geospatial_vertical_units',geospatial.vertical_units);
            netcdf.putAtt(ncid,ncglobal,'geospatial_vertical_positive',geospatial.vertical_positive);
            % Take file out of define mode.
            netcdf.endDef(ncid);
            
            % Synchronize the file to disk.
            netcdf.sync(ncid)
                        
        end
        
        function fileHeader = initHeader(xmlfile)
            
            fileHeader = struct(...
                   'title','',...
                      'id','',...
           'data_revision','',...
                'summary',[''...
                           ' '...
                           ' '...
                           ''],...
                'keywords','',...
             'Conventions','CF-1.0',...
    'Metadata_Conventions','Unidata Dataset Discovery v1.0',...
             'author_name','',...
              'author_url','',...
            'author_email','',...
             'institution','',...
         'repository_name','',...
  'repository_institution','',...
          'repository_pid','',...
         'acknowledgement','',...
              'references','',...
                 'history','',...
                 'comment','');
             
             if nargin > 0
                 fileHeader = xml_read(xmlfile);
             end
             
        end

        
        function geospatialInfo = initGeospatial(limits,res)
            
    geospatialInfo = struct(...
              'lat_min', '-90.00',...
              'lat_max', '90.00',...
            'lat_units', 'degrees_north',...
       'lat_resolution', '2',...
              'lon_min', '-180.00',...
              'lon_max', '180.00',...
            'lon_units', 'degrees_east',...
       'lon_resolution', '2',...
         'vertical_min', '0',...
         'vertical_max', '3500',...
       'vertical_units', 'km',...
    'vertical_positive', 'down');

             if nargin > 0
                geospatialInfo.lon_min = sprintf('%6.1f',netcdfiris.lon180(limits.lonmin));
                geospatialInfo.lon_max = sprintf('%6.1f',netcdfiris.lon180(limits.lonmax));
                geospatialInfo.lat_min = sprintf('%6.1f',limits.latmin);
                geospatialInfo.lat_max = sprintf('%6.1f',limits.latmax);
                geospatialInfo.vertical_min = sprintf('%4.0f',limits.depthmin);
                geospatialInfo.vertical_max = sprintf('%4.0f',limits.depthmax);
             end
             
             if nargin > 1
                 if isscalar(res)
                     geospatialInfo.lon_resolution = sprintf('%4.2f',res);
                     geospatialInfo.lat_resolution = sprintf('%4.2f',res);
                 else
                     geospatialInfo.lon_resolution = sprintf('%4.2f',res(1));
                     geospatialInfo.lat_resolution = sprintf('%4.2f',res(2));
                 end
             end

        end
        
        function lon = lon180(lon)
            %   given list of longitudes, convert them to (-180,180]
            lon(lon>180) = lon(lon>180)-360;
        end
        
    end
    
end