function [status] = writeZ_3D(cfile, allData, header, units, isign, convert_to_apres)
%  Usage:  [status] = writeZ_3D(cfile, allData, header, units, isign, convert_to_apres);
%   write contents of cell array allData to file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   Site locations units have to match the model, i.e. use meters.
%   Last mod. in Dec 2016 to allow a data object instead of data structure.
%  (c) Anna Kelbert, 2011

% if no txType, assume old data file format that does not specify it
txType = '';
if isfield(allData{1},'txType') || isprop(allData{1},'txType')
    txType = allData{1}.txType;
end

if nargin <= 5
    convert_to_apres = 0;
end

[dataTypes, components, dt_comp_map] = determine_datatypes_n_comps(allData);

comps = components;
ncomp = length(comps);

%  compute the total number of frequencies and sites
nTx = length(allData);
allsites = allData{1}.siteChar;
allsitesloc = allData{1}.siteLoc;
if isfield(allData{1},'lat') || isprop(allData{1},'lat')
    GG = 1;
    lat = allData{1}.lat;
    lon = allData{1}.lon;
else
    GG = 0;
end
for j = 2:nTx
    for i = 1:length(allData{j}.siteChar)
        if isempty(find(contains(allsites,allData{j}.siteChar{i}), 1))
            allsites = [allsites; allData{j}.siteChar{i}];
            allsitesloc = [allsitesloc; allData{j}.siteLoc(i,:)];
            if GG
                lat = [lat; allData{j}.lat(i)];
                lon = [lon; allData{j}.lon(i)];
            end
        end
    end
end
nSites = length(allsites);
if ~ GG %  if GG locations are not supplied, set them to zero
    lat(1:nSites) = 0;
    lon(1:nSites) = 0;
end

%  regroup into observatory bins
info.data = nan(nSites,nTx,ncomp);
info.err = nan(nSites,nTx,ncomp);

for j = 1:nTx
    for i = 1:length(allsites)
        k=find(contains(allData{j}.siteChar,allsites{i}));
        if ~isempty(k)
            info.code{i} = allsites{i};
            info.loc(i,:) = allsitesloc(i,:);
            info.lon(i) = lon(i);
            info.lat(i) = lat(i);
            if isprop(allData{j},'TF')
                info.data(i,j,:) = allData{j}.TF(k,:);
                info.err(i,j,:) = allData{j}.TFerr(k,:);
            else
                info.data(i,j,:) = allData{j}.Z(k,:);
                info.err(i,j,:) = allData{j}.Zerr(k,:);
            end
            if sum(isnan(info.err(i,j,:)))>0 % fill in NaN errors 
                %info.err(i,j,:) = allData{j}.nanvalue;
                info.err(i,j,:) = NaN;
            end
            info.per(j) = allData{j}.T;
        end
    end
end
info.comp = comps;

fid = fopen(cfile,'w');

%  description <= 80 char in length
if nargin < 3
    header = 'Synthetic 3D MT data written in Matlab';
end

%  assume SI units by default as in ModEM; alternative is [mV/km]/[nT] 
if nargin < 4
    if isfield(allData{1},'units') || isprop(allData{1},'units')
        units = allData{1}.units;
    else
        units = '[V/m]/[T]';
    end
end

%  sign convention is -1 by default
if nargin < 5
    if isfield(allData{1},'signConvention') || isprop(allData{1},'signConvention')
        isign = allData{1}.signConvention;
    else
        isign = -1;
    end
end
if isign == -1
    signstr = 'exp(-i\omega t)';
else
    signstr = 'exp(+i\omega t)';
end

%  get the origin from the first period
origin = [0 0 0];
if isfield(allData{1},'origin') || isprop(allData{1},'origin')
    origin = allData{1}.origin;
end

%  get orientation from the first period
orientation = 0.0;
if isfield(allData{1},'orient') || isprop(allData{1},'orient')
    orientation = allData{1}.orient;
end

% compute apparent resistivity and phase for writing, if required
if convert_to_apres
    [apres, phase, dataType] = convert_impedance_to_apres(dataTypes, info);
end

% IMPORTANT: NEW LINE IN THE DATA FILE DEFINES THE TRANSMITTER TYPE
% (this is omitted when txType = '' for backward compatibility)
if ~isempty(txType)
    fprintf(fid,'+ %s\n',txType);
end
% END IMPORTANT: NEW LINE IN THE DATA FILE


for nTypes = 1 : length(dataTypes)
    if nargin < 4
        if isfield(allData{1},'units') || isprop(allData{1},'units')
            units = allData{1}.units;
        else
            units = '[V/m]/[T]';
        end
    end

    if isempty(units)
        units = '[V/m]/[T]';
    end

    % Comment line ignore by code, but describes each data in the item block
    if strcmp(dataTypes(nTypes), 'Off_Diagonal_Rho_Phase') || strcmp(dataTypes(nTypes), 'Phase_Tensor')
        comment = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Error';
    else
        comment = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error';
    end

    writeHeader(fid, header, comment, dataTypes(nTypes), signstr, units, orientation, origin, nTx, nSites)

    %  now write all impedances (or apparent resistivities!) line by line
    if convert_to_apres 
        writeApres(fid, nSites, nTx, info, phase, apres)
    elseif strcmp(dataTypes(nTypes), 'Off_Diagonal_Rho_Phase')
        [apres, phase] = split_into_apres_n_phase(info);
        writeApres(fid, nSites, nTx, info, phase, apres)
    elseif strcmp(dataTypes(nTypes), 'Phase_Tensor')
        writePhaseTensor(fid, nSites, nTx, info)
    else
        [comp_start, comp_end] = determine_comp_loc(info, dataTypes(nTypes), dt_comp_map);
        for k = 1:nSites
            for j = 1:nTx
                for i = comp_start : comp_end
                    if ~isnan(info.data(k,j,i))
                        fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                        fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                        fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                        fprintf(fid,'%s %15.6E %15.6E %15.6E\n',info.comp(i,:),real(info.data(k,j,i)),imag(info.data(k,j,i)),info.err(k,j,i)); % data
                    end
                end
            end
        end
    end

end

status = fclose(fid);
end

function [] = writeHeader(fid, header, comment, datatype, signstr, units, orientation, origin, nTx, nSites)
    fprintf(fid,'# %s\n', header);
    fprintf(fid,'# %s\n', comment);
    fprintf(fid,'> %s\n', datatype);
    fprintf(fid,'> %s\n', signstr);
    fprintf(fid,'> %s\n', units);
    fprintf(fid,'> %.2f\n', orientation);
    fprintf(fid,'> %.3f %.3f\n', origin(1:2));
    fprintf(fid,'> %d %d\n', nTx, nSites);
end

function [] = writeZ(fid, nSites, nTx, info)
    for k = 1 : nSites
        for j = 1 : nTx
            for i = nComps
                if ~isnan(info.data(nSite, tx, comp))
                    fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                    fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                    fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                    fprintf(fid,'%s %15.6E %15.6E %15.6E\n',info.comp(i,:),real(info.data(k,j,i)),imag(info.data(k,j,i)),info.err(k,j,i)); % data
                end
            end
        end

    end
end

% Convert Impedeance datatypes to apparent resistivity and phase
function [apres, phase, dataType] = convert_impedance_to_apres(dataType, info)
    rad_deg = 57.2958;

    if strcmp(dataType,'Off_Diagonal_Impedance')
        ixy = 1;
        iyx = 2;
    else
        ixy = 2;
        iyx = 3;
    end

    apres.xy = abs(info.data(:,:,ixy)).^2;
    apres.xy_se = info.err(:,:,ixy).^2; % variance of Z as in EDI file
    apres.yx = abs(info.data(:,:,iyx)).^2;
    apres.yx_se = info.err(:,:,iyx).^2; % variance of Z as in EDI file

    %phase.xy = rad_deg*atan((imag(info.data(:,:,ixy)))./real(info.data(:,:,ixy)));
    % for now, using atan2 which takes values (-pi,pi] as in ModEM
    phase.xy = rad_deg*atan2(imag(info.data(:,:,ixy)),real(info.data(:,:,ixy)));
    phase.xy_se = rad_deg*sqrt(apres.xy_se./apres.xy);

    %phase.yx = rad_deg*atan((imag(info.data(:,:,iyx)))./real(info.data(:,:,iyx)));
    % for now, using atan2 which takes values (-pi,pi] as in ModEM
    phase.yx = rad_deg*atan2(imag(info.data(:,:,iyx)),real(info.data(:,:,iyx))) + 180.;
    phase.yx_se = rad_deg*sqrt(apres.yx_se./apres.yx);

    % rescale apparent resistivity by period
    for l = 1:length(info.per)
        apres.yx(:,l) = apres.yx(:,l)*info.per(l)/5. ;
        apres.xy(:,l) = apres.xy(:,l)*info.per(l)/5. ;
        apres.yx_se(:,l) = sqrt(apres.yx_se(:,l).*apres.yx(:,l)*info.per(l)*4/5.);
        apres.xy_se(:,l) = sqrt(apres.xy_se(:,l).*apres.xy(:,l)*info.per(l)*4/5.);
    end

    dataType = 'Off_Diagonal_Rho_Phase';
end

% To write out Phase Tensor and Off_Diagonal_Rho_Phase, convert to apres/phase data structures
% that this file uses
function [apres, phase] = split_into_apres_n_phase(info)
    apres.xy = info.data(:,:,1);
    apres.xy_se = info.err(:,:,1);
    apres.yx = info.data(:,:,3);
    apres.yx_se = info.err(:,:,3);

    phase.xy = info.data(:,:,2);
    phase.xy_se = info.err(:,:,2);
    phase.yx = info.data(:,:,4);
    phase.yx_se = info.err(:,:,4);
end

function [] = writeApres(fid, nSites, nTx, info, phase, apres)
    for k = 1:nSites
        for j = 1:nTx
            if ~isnan(apres.xy(k,j))
                fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                fprintf(fid,'%s %15.6E %15.6E\n','RHOXY',apres.xy(k,j),apres.xy_se(k,j)); % data
            end
            if ~isnan(phase.xy(k,j))
                fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                fprintf(fid,'%s %15.6E %15.6E\n','PHSXY',phase.xy(k,j),phase.xy_se(k,j)); % data
            end
            if ~isnan(apres.yx(k,j))
                fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                fprintf(fid,'%s %15.6E %15.6E\n','RHOYX',apres.yx(k,j),apres.yx_se(k,j)); % data
            end
            if ~isnan(phase.yx(k,j))
                fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                fprintf(fid,'%s %15.6E %15.6E\n','PHSYX',phase.yx(k,j),phase.yx_se(k,j)); % data
            end
        end
    end
end

function [] = writePhaseTensor(fid, nSites, nTx, info)
    phase.xx = info.data(:,:,1);
    phase.xx_se = info.err(:,:,1);

    phase.xy = info.data(:,:,2);
    phase.xy_se = info.err(:,:,2);

    phase.yx = info.data(:,:,3);
    phase.yx_se = info.err(:,:,3);

    phase.yy = info.data(:,:,4);
    phase.yy_se = info.err(:,:,4);

    for k = 1:nSites
        for j = 1:nTx
            if ~isnan(phase.xx(k,j))
                fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                fprintf(fid,'%s %15.6E %15.6E\n','PTXX',phase.xx(k,j),phase.xx_se(k,j)); % data
            end
            if ~isnan(phase.xy(k,j))
                fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                fprintf(fid,'%s %15.6E %15.6E\n','PTXY',phase.xy(k,j),phase.xy_se(k,j)); % data
            end
            if ~isnan(phase.yx(k,j))
                fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                fprintf(fid,'%s %15.6E %15.6E\n','PTYX',phase.yx(k,j),phase.yx_se(k,j)); % data
            end
            if ~isnan(phase.yy(k,j))
                fprintf(fid,'%12.6E ',info.per(j)); % transmitter
                fprintf(fid,'%s %8.3f %8.3f ',info.code{k},info.lat(k),info.lon(k)); % receiver
                fprintf(fid,'%12.3f %12.3f %12.3f ',info.loc(k,:)); % receiver x,y,z
                fprintf(fid,'%s %15.6E %15.6E\n','PTYY',phase.yy(k,j),phase.yy_se(k,j)); % data
            end
        end
    end
end

function [dataTypes, components, dt_comp_map] = determine_datatypes_n_comps(allData)
    if isfield(allData{1}, 'type') || isprop(allData{1}, 'type')
        [dataTypes, components, dt_comp_map] = comps_from_type(allData);
    else
        [dataTypes, components, dt_comp_map] = types_from_comps(allData);
    end
end

function [dataTypes, components, dt_comp_map] = comps_from_type(allData)
    components = strings(0);
    dataTypes = strings(0);
    dt_comp_map = dictionary();

    for nperiod=1:length(allData)
        period_types = allData{nperiod}.type;
        if isa(period_types, "char")
            period_types = convertCharsToStrings(period_types);
        end

        for i = 1 : length(period_types)
            switch period_types(i)
                case 'Full_Impedance'
                    comps = ['ZXX'; 'ZXY'; 'ZYX'; 'ZYY'];
                case 'Off_Diagonal_Impedance'
                    comps = ['ZXY'; 'ZYX'];
                case 'Full_Vertical_Components'
                    comps = ['TX'; 'TY'];
                case 'Full_Interstation_TF'
                    comps = ['MXX'; 'MXY'; 'MYX'; 'MYY'];
                case 'Off_Diagonal_Rho_Phase'
                    comps = ['RHOXY'; 'PHSXY'; 'RHOYX'; 'PHSYX'];
                case 'Phase_Tensor'
                    comps = ['PTXX'; 'PTXY'; 'PTYX'; 'PTYY'];
                otherwise
                    error("Unkown datatype '%s'", period_types)
                    % Potentially determine by components
                    % Potentially fix by requesting dataType in allData
            end
            
            dt_comp_map(period_types(i)) = {comps};
    
            if isempty(intersect(dataTypes, period_types))
                dataTypes = [dataTypes, period_types];
            end
    
            if isempty(intersect(comps, components))
                % Add the list of components to comps
                components = vertcat(components, comps);
            end
        end
    end
end

function [dataTypes, components, dt_comp_map] = types_from_comps(allData)
    dataTypes = strings(0);
    components = strings(0);
    dt_comp_map = dictionary();

    for nperiod = 1 : length(allData)
        for ncomp = 1 : length(allData{nperiod}.compChar)
            comp = allData{nperiod}.compChar(ncomp,:);

            if ~isempty(intersect(comp, ['ZXX'; 'ZXY'; 'ZYX'; 'ZYY'], 'rows'))
                type = 'Full_Impedance';
                comps = ['ZXX'; 'ZXY'; 'ZYX'; 'ZYY'];

            elseif ~isempty(intersect(comp, ['ZXY'; 'ZYX'], 'rows'))
                type = 'Off_Diagonal_Impedance';
                comps = ['ZXY'; 'ZYX'];
            
            elseif ~isempty(intersect(comp, ['TX '; 'TY '], 'rows'))
                type = 'Full_Vertical_Components';
                comps = ['TX ' ; 'TY '];

            elseif ~isempty(intersect(comp, ['MXX'; 'MXY'; 'MYX'; 'MYY'], 'rows'))
                type = 'Full_Interstation_TF';
                comps = ['MXX' ; 'MXY' ; 'MYX' ; 'MYY'];

            elseif ~isempty(intersect(comp, ['RHOXY'; 'PHSXY'; 'RHOYX'; 'PHSYX'], 'rows'))
                type = 'Off_Diagonal_Rho_Phase';
                comps = ['RHOXY'; 'PHSXY'; 'RHOYX' ; 'PHSYX'];

            elseif ~isempty(intersect(comp, ['PTXX'; 'PTXY'; 'PTYX'; 'PTYY'], 'rows'))
                type = 'Phase_Tensor';
                comps = ['PTXX'; 'PTXY' ; 'PTYX' ; 'PTYY'];
            else
                error("Did not recognize this component")
            end

            dt_comp_map(type) = {comps};

            if isempty(intersect(dataTypes, type))
                dataTypes = [dataTypes, type];
            end

            if ~all(ismember(comps, components))
                components = vertcat(components, comps);
            end
        end
    end
end

% Utility function to determine order of datatypes inside info.data. Info.data
% contains and info.comp contains *all* datatypes (if there are multiple). This
% function will return the starting and ending location in info.data and
% info.comp of the requested datatype.
function [start_loc, end_loc] = determine_comp_loc(info, datatype, dt_comp_dict)
    start_loc = 1000;
    end_loc = 0;
    comps = dt_comp_dict(datatype);
    for i = 1 : size(comps{1}, 1)
        for j = 1 : length(info.comp)
            if strcmp(comps{1}(i,:), info.comp(j))
                if j < start_loc
                    start_loc = j;
                end

                if j > end_loc
                    end_loc = j;
                end
            end
        end
    end
end
