function [allData, header, units, isign, origin, info] = readZ_3D(cfile, newunits, onetype, nanvalue, newformat)
% readZ_3D - Read all data blocks of a ModEM 'list'format
%  Usage: [allData, header, units, isign, origin, info] = readZ_3D(cfile, newunits, onetype, nanvalue, newformat)
%
%  Read contents of a ModEM 'list' datafile into the cell array allData
% from cfile. Each cell represents one period and contains all information
% necessary to define data (loctions, values, error standard dev, components,
% dataType names) for a transmitter.
%
% Site locations units have to match the model, i.e. use meters.
%
% Input Arguments
%  cfile - Required - Filename to read from
%   string value
%  newunits - Optional - Units to convert data too either: '[mV/km]/[nT]', '[V/m]/[T]', '[V/m]/[A/m]', 'Ohm', or ''.
%       Unit conversion is performed by ImpUnits
%   string value
%  onetype - Optional - If present, only read that data type and skip others
%   string value
%  nanvalue - Optional - If present, use this to value in place of missing data or error values
%   float value or NaN
%  newformat - Optional - If present, read data with azimuth angles
%   boolean value
%
% Output Arguments
%   allData - Cell array, one cell per period, containing all data for all data types that were
%           in that had data for that period
%    cell array
%   header - The comment for the file (very first line)
%    string value
%   units - Units that the data has been converted too, or if 'newunits' was not present, the
%       units associated with the file
%    string value
%   isign - The sign of the data either -1 or 1
%    integer value
%   origin - [X,Y] values of data origin
%    double array
%   info - Cell array, one cell per data type containing all data of data type including
%     data, err, units, and site location (lat, lon, xyz location)
%    cell array
%
% See also readApres_3D, writeZ_3D, ImpUnits, mtdata
%  (c) Anna Kelbert, 2011-2013, 2020, 2023

if ~isfile(cfile)
    error("File '%s' could not be found.", cfile);
end

[fid, err_msg] = fopen(cfile,'r');
if fid == -1
    error("Error reading file '%s': %s", cfile, err_msg)
end

if nargin < 2
    newunits = '';
end

if nargin < 3
    onetype = '';
end

if nargin < 4
    nanvalue = NaN;
else
    if ~isnumeric(nanvalue)
        nanvalue = str2double(nanvalue);
    end
end

if nargin < 5
    newformat = 0;
end

origin = [NaN NaN];
nDataTypes = 0;
seenDataTypes = strings(0);

%  read the data: one block per data type
while 1
    % read the header info - with different compilers, various ends of lines
    header = '';
    l=fgetl(fid);                % read first line
    if l < 0; break; end         % end of file

    while ~contains(string(l),'>')       % loop until find the magic record
        header = [header l];
        l=fgetl(fid);
    end
    dataType = sscanf(l,'> %s');
    seenDataTypes = [seenDataTypes, dataType];

    blockinfo = textscan(fid,'> %s',5,'delimiter','\n');
    blockinfo = char(blockinfo{1});
    signstr  = blockinfo(1,:);
    if findstr(signstr,'-')
        isign = -1;
    else
        isign = +1;
    end
    typeUnits = strtrim(blockinfo(2,:));
    orientation = sscanf(blockinfo(3,:),'%f',1);
    neworigin = sscanf(blockinfo(4,:),'%f',2);
    tmp  = sscanf(blockinfo(5,:),'%d',2);
    nTx = tmp(1);
    nSites = tmp(2);
    % now read the data line by line
    switch strtrim(dataType)
        case 'Phase_Tensor'
            warning("Skipping read of datatype '%s'. Use 'readApres_3D(cfile, newunits, onetype)' for 'Phase_Tensor' datatypes", dataType);
            continue
        case 'Full_Impedance'
            nDataTypes = nDataTypes + 1;
            origin = neworigin; % prevents problem when origins differ
            if ~isempty(newunits)
                SI_factor = ImpUnits(typeUnits,newunits);
                units = newunits;
            else
                SI_factor = 1.0;
                units = typeUnits; % output the units of impedances
            end
            ncomp = 4;
            comp = ['ZXX';'ZXY';'ZYX';'ZYY'];
            info{nDataTypes}.data = nan(nSites,nTx,ncomp)+1i*nan(nSites,nTx,ncomp);
            info{nDataTypes}.err = nan(nSites,nTx,ncomp);
            info{nDataTypes}.type = 'Full_Impedance';
            info{nDataTypes}.units = units;
            if ~newformat
                % BEFORE ADDING AZIMUTH ANGLES
                data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            else
                % WITH THE AZIMUTH ANGLES
                data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f %f %f %f %f');
            end
            periods = unique(data{1});
            codes = sortrows(strtrim(char(unique(data{2}))));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist,'exact');
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            iskip = isequal(data{9}(ind),nanvalue) || ...
                                isequal(data{10}(ind),nanvalue) || ...
                                isequal(data{11}(ind),nanvalue);
                            if ~iskip
                                try
                                    info{nDataTypes}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                                    info{nDataTypes}.err(k,j,i) = data{11}(ind);
                                    info{nDataTypes}.lat(k) = data{3}(ind);
                                    info{nDataTypes}.lon(k) = data{4}(ind);
                                    info{nDataTypes}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                                catch
                                    tmp = char(data{2}(ind));
                                    warning(['Problem from line # ' num2str(ind(1)) ' site ' tmp(1,:) ' - possibly a duplicate site, using occurence 1']);
                                    info{nDataTypes}.data(k,j,i) = data{9}(ind(1)) + 1i* data{10}(ind(1));
                                    info{nDataTypes}.err(k,j,i) = data{11}(ind(1));
                                    info{nDataTypes}.lat(k) = data{3}(ind(1));
                                    info{nDataTypes}.lon(k) = data{4}(ind(1));
                                    info{nDataTypes}.loc(k,1:3) = [data{5}(ind(1)) data{6}(ind(1)) data{7}(ind(1))];
                                end
                            else
                                disp(['Skipping component ',strtrim(comp(i,:)),' for period ',num2str(j),' of ',codes(k,:)]);
                            end
                        end
                    end
                end
            end
            info{nDataTypes}.code = codes;
            info{nDataTypes}.per = periods;
            info{nDataTypes}.ncomp = 8;
            info{nDataTypes}.comp = comp;
        case 'Off_Diagonal_Impedance'
            nDataTypes = nDataTypes + 1;

            if isnan(origin)
                origin = neworigin; % prevents problem when origins differ
            end
            if ~isempty(newunits)
                SI_factor = ImpUnits(typeUnits,newunits);
                units = newunits;
            else
                SI_factor = 1.0;
                units = typeUnits; % output the units of impedances
            end
            ncomp = 2;
            comp = ['ZXY';'ZYX'];
            info{nDataTypes}.data = nan(nSites,nTx,ncomp);
            info{nDataTypes}.err = nan(nSites,nTx,ncomp);
            info{nDataTypes}.type = 'Off_Diagonal_Impedance';
            info{nDataTypes}.units = units;
            if ~newformat
                % BEFORE ADDING AZIMUTH ANGLES
                data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            else
                % WITH THE AZIMUTH ANGLES
                data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f %f %f %f %f');
            end            
            periods = unique(data{1});
            codes = sortrows(strtrim(char(unique(data{2}))));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            info{nDataTypes}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                            info{nDataTypes}.err(k,j,i) = data{11}(ind);
                            info{nDataTypes}.lat(k) = data{3}(ind);
                            info{nDataTypes}.lon(k) = data{4}(ind);
                            info{nDataTypes}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                        end
                    end
                end
            end
            info{nDataTypes}.code = codes;
            info{nDataTypes}.per = periods;
            info{nDataTypes}.ncomp = 4;
            info{nDataTypes}.comp = comp;
        case 'Off_Diagonal_Rho_Phase'
            warning("Will convert 'Off_Diagonal_Rho_Phase' into Full_Impedance. Use 'readApres_3D' to read 'Off_Diagonal_Rho_Phase' without conversion")
            nDataTypes = nDataTypes + 1;
            if isnan(origin)
                origin = neworigin; % prevents problem when origins differ
            end
            if ~isempty(newunits)
                SI_factor = ImpUnits(typeUnits,newunits);
                units = newunits;
            else
                SI_factor = 1.0;
                units = typeUnits; % output the units of impedances
            end
            ncomp = 4;
            comp = ['RHOXY';'PHSXY';'RHOYX';'PHSYX'];
            info{nDataTypes}.data = nan(nSites,nTx,ncomp);
            info{nDataTypes}.err = nan(nSites,nTx,ncomp);
            info{nDataTypes}.units = units;
            if ~newformat
                % BEFORE ADDING AZIMUTH ANGLES
                data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f %f %f');
            else
                % WITH THE AZIMUTH ANGLES
                data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f %f %f %f %f %f %f');
            end
            
            periods = unique(data{1});
            codes = sortrows(strtrim(char(unique(data{2}))));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            info{nDataTypes}.data(k,j,i) = data{9}(ind);
                            info{nDataTypes}.err(k,j,i) = data{10}(ind);
                            info{nDataTypes}.lat(k) = data{3}(ind);
                            info{nDataTypes}.lon(k) = data{4}(ind);
                            info{nDataTypes}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                        end
                    end
                end
            end
            % convert to Off_Diagonal_Impedance!!! (need a separate
            % function to read Apres+Phase and a separate class to store)
            for k = 1:length(codes)
                apres.xy = info{nDataTypes}.data(k,:,1);
                apres.xy_se = info{nDataTypes}.err(k,:,1);
                phase.xy = info{nDataTypes}.data(k,:,2);
                phase.xy_se = info{nDataTypes}.err(k,:,2);
                apres.yx = info{nDataTypes}.data(k,:,3);
                apres.yx_se = info{nDataTypes}.err(k,:,3);
                phase.yx = info{nDataTypes}.data(k,:,4);
                phase.yx_se = info{nDataTypes}.err(k,:,4);
                [Z,Zstd] = mttf.apres2imp(periods,apres,phase);
                info{nDataTypes}.data(k,:,1) = Z(1,1,:);
                info{nDataTypes}.data(k,:,2) = Z(1,2,:);
                info{nDataTypes}.data(k,:,3) = Z(2,1,:);
                info{nDataTypes}.data(k,:,4) = Z(2,2,:);
                info{nDataTypes}.err(k,:,1) = Zstd(1,1,:);
                info{nDataTypes}.err(k,:,2) = Zstd(1,2,:);
                info{nDataTypes}.err(k,:,3) = Zstd(2,1,:);
                info{nDataTypes}.err(k,:,4) = Zstd(2,2,:);
            end
            info{nDataTypes}.type = 'Full_Impedance';
            info{nDataTypes}.code = codes;
            info{nDataTypes}.per = periods;
            info{nDataTypes}.ncomp = 8;
            info{nDataTypes}.comp = ['ZXX';'ZXY';'ZYX';'ZYY'];
        case 'Full_Vertical_Components'
            nDataTypes = nDataTypes + 1;

            if isnan(origin)
                origin = neworigin; % prevents problem when origins differ
            end
            if ~isempty(newunits)
                SI_factor = ImpUnits(typeUnits, newunits);
                units = newunits;
            else
                SI_factor = 1.0;
                units = typeUnits;
            end

            ncomp = 2;
            comp = ['TX ';'TY '];
            info{nDataTypes}.type = 'Full_Vertical_Components';
            info{nDataTypes}.data = nan(nSites,nTx,ncomp);
            info{nDataTypes}.err = nan(nSites,nTx,ncomp);
            info{nDataTypes}.units = units;
            if ~newformat
                % BEFORE ADDING AZIMUTH ANGLES
                data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            else
                % WITH THE AZIMUTH ANGLES
                data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f %f %f %f %f');
            end
            periods = unique(data{1});
            codes = sortrows(char(unique(data{2})));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            try
                                info{nDataTypes}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                                info{nDataTypes}.err(k,j,i) = data{11}(ind);
                                info{nDataTypes}.lat(k) = data{3}(ind);
                                info{nDataTypes}.lon(k) = data{4}(ind);
                                info{nDataTypes}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                            catch
                                tmp = char(data{2}(ind));
                                warning(['Problem from line # ' num2str(ind(1)) ' site ' tmp(1,:) ' - possibly a duplicate site, using occurence 1']);
                                info{nDataTypes}.data(k,j,i) = data{9}(ind(1)) + 1i* data{10}(ind(1));
                                info{nDataTypes}.err(k,j,i) = data{11}(ind(1));
                                info{nDataTypes}.lat(k) = data{3}(ind(1));
                                info{nDataTypes}.lon(k) = data{4}(ind(1));
                                info{nDataTypes}.loc(k,1:3) = [data{5}(ind(1)) data{6}(ind(1)) data{7}(ind(1))];
                            end
                        end
                    end
                end
            end
            info{nDataTypes}.code = codes;
            info{nDataTypes}.per = periods;
            info{nDataTypes}.ncomp = 4;
            info{nDataTypes}.comp = comp;
        otherwise
            warning("This datatype (%s) is unkown to 'readZ_3D. Skipping", dataType)
            continue;
    end   
end
status = fclose(fid);
if status == -1
    warning("Error calling fclose of for file '%s'", cfile)
end

% If requested, only return one datatype
if ~isempty(onetype)
    tmp_info = {};
    for i = 1:length(info)
        if strcmp(onetype, strtrim(info{i}.type))
            tmp_info{1} = info{i};
            units = info{i}.units;
        end
    end
    info = tmp_info;

    if isempty(info)
        error_msg = sprintf("Requested onetype '%s' is not present in the file '%s'", onetype, cfile);
        resolve_message = sprintf("Please choose a datatype that is in this file: '[ %s ]'", join(seenDataTypes, " | "));
        error(error_msg + newline + newline + resolve_message);
    end
end

% If more than two data types are present, merge the periods
if isscalar(info)
    per = info{1}.per;
    ind1 = 1:length(per);
else
    per = info{1}.per;
    for i = 2:length(info)
        per = union(per, info{i}.per);
    end
    [tmp,ind1] = intersect(per,info{1}.per);
    [tmp,ind2] = intersect(per,info{2}.per);
end

% compute the maximum total number of components
ncomp = 0;
for i = 1:length(info)
    ncomp = ncomp + info{i}.ncomp/2;
end

if length(info) == 2
    % If we see full vertical components before full impedance, swap them so
    % that all data has them in the order described below
    if strcmp(info{1}.type, 'Full_Vertical_Components') && strcmp(info{2}.type, 'Full_Impedance')
        tmp_full_vert = info{1};
        info{1} = info{2};
        info{2} = tmp_full_vert;
    end
end


% For compatibility, convert to the (old) allData structure - most programs
% expect full impedances and vertical components (in this order) 
for j = 1:length(per)
    for k = 1:length(info)
        ind = find(info{k}.per==per(j));
        if ~isempty(ind); itx(k) = ind; else itx(k) = 0; end
    end    
    allsites = info{1}.code;
    allsitesloc = info{1}.loc;
    allsiteslat = info{1}.lat;
    allsiteslon = info{1}.lon;
    for k = 2:length(info)
        for i = 1:length(info{k}.code)
            if isempty(intersect(allsites,info{k}.code(i,:),'rows'))
                allsites = [allsites; info{k}.code(i,:)];
                allsitesloc = [allsitesloc; info{k}.loc(i,:)];
                allsiteslat = [allsiteslat; info{k}.lat(i)];
                allsiteslon = [allsiteslon; info{k}.lon(i)];
            end
        end
    end
    nsites = length(allsites);
    allData{j}.T = per(j);
    allData{j}.Cmplx = 1;
    allData{j}.signConvention = isign;
    allData{j}.nComp = ncomp;
    allData{j}.nanvalue = nanvalue;
    allData{j}.siteLoc = allsitesloc;
    for isite = 1:size(allsites,1) % quick fix as of 20190301. Clean up later
        allData{j}.siteChar{isite} = allsites(isite,:);
    end
    allData{j}.Z(1:nsites,1:ncomp) = NaN;
    allData{j}.Zerr(1:nsites,1:ncomp) = NaN;
    allData{j}.nComp = 2*ncomp;
    allData{j}.origin = [origin' 0];
    allData{j}.orient = orientation;
    allData{j}.lat = allsiteslat';
    allData{j}.lon = allsiteslon';
    allData{j}.type = strings(0);
    icomp1 = 1;
    for k = 1:length(info)
        if itx(k) > 0
            if allData{j}.Cmplx 
                icomp2 = icomp1 + info{k}.ncomp/2 - 1;
            else
                icomp2 = icomp1 + info{k}.ncomp - 1;
            end                
            [sites,irx] = intersect(allsites,info{k}.code,'rows');
            allData{j}.Z(irx,icomp1:icomp2) = SI_factor*squeeze(info{k}.data(:,itx(k),:));
            allData{j}.Zerr(irx,icomp1:icomp2) = SI_factor*squeeze(info{k}.err(:,itx(k),:));
            allData{j}.compChar(icomp1:icomp2,:) = info{k}.comp;
            allData{j}.type = [allData{j}.type, info{k}.type];
            allData{j}.units = units;
            icomp1 = icomp2+1;
        end
    end
end