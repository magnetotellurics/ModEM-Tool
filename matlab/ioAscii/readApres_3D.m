function [allData,header,units,isign,origin,info] = readApres_3D(cfile,newunits,onetype)
%  Usage:  [allData,header,units,isign,origin,info] = readApres_3D(cfile,newunits,onetype);
%   read contents of cell array allData from file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   Site locations units have to match the model, i.e. use meters.
%   If onetype is specified, only read that data type and skip others.
%  Same as readZ_3D but for the real apparent resistivity and phase.
%  (c) Anna Kelbert, 2011-2013


fid = fopen(cfile,'r');

if nargin < 3
    onetype = '';
end

origin = [NaN NaN];
nDataTypes = 1;

compCharMaxSize = 0;

%  read the data: one block per data type
while 1
    % read the header info
    tmp = textscan(fid,'# %s',1,'delimiter','\n');
    tmp = char(tmp{1});
    if isempty(tmp); break; end
    header = tmp;

    % read the block header
    tmp = textscan(fid,'# %s',1,'delimiter','\n');

    % IMPORTANT: NEW LINE IN THE DATA FILE DEFINES THE TRANSMITTER TYPE
    % (if missing, just skips over it - can still read old files)
    tmp = textscan(fid,'+ %s',1,'delimiter','\n');
    % END IMPORTANT: NEW LINE IN THE DATA FILE

    blockinfo = textscan(fid,'> %s',6,'delimiter','\n');
    blockinfo = char(blockinfo{1});
    dataType = blockinfo(1,:);
    signstr  = blockinfo(2,:);

    if findstr(signstr,'-')
        isign = -1;
    else
        isign = +1;
    end

    typeUnits = strtrim(blockinfo(3,:));
    orientation = sscanf(blockinfo(4,:),'%f',1);
    neworigin = sscanf(blockinfo(5,:),'%f',2);
    tmp  = sscanf(blockinfo(6,:),'%d',2);
    nTx = tmp(1);
    nSites = tmp(2);

    % now read the data line by line
    switch strtrim(dataType)
        case {'Full_Impedance','Off_Diagonal_Impedance','Full_Vertical_Components'}
            warning("Skipping read of datatype '%s'. Use 'readZ_3D(cfile,newunits,onetype)' for complex data types.", strtrim(dataType));
            continue
        case 'Off_Diagonal_Rho_Phase'
            if 5 > compCharMaxSize
                compCharMaxSize = 5;
            end

            if ~isempty(onetype) && ~strcmp(onetype,strtrim(dataType))
                continue
            end
            if isnan(origin)
                origin = neworigin; % prevents problem when origins differ
            end
            if nargin > 1
                SI_factor = ImpUnits(typeUnits,newunits);
                units = newunits;
            else
                SI_factor = 1.0;
                units = typeUnits; % output the units of impedances
            end
            ncomp = 4;
            comp = ['RHOXY';'PHSXY';'RHOYX';'PHSYX'];
            info{nDataTypes}.type = 'Off_Diagonal_Rho_Phase';
            info{nDataTypes}.data = nan(nSites,nTx,ncomp);
            info{nDataTypes}.err = nan(nSites,nTx,ncomp);
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f');

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
            info{nDataTypes}.code = codes;
            info{nDataTypes}.per = periods;
            info{nDataTypes}.ncomp = 4;
            info{nDataTypes}.comp = comp;
            nDataTypes = nDataTypes + 1;
        case 'Phase_Tensor'
            if 4 > compCharMaxSize
                compCharMaxSize = 4;
            end

            if ~isempty(onetype) && ~strcmp(onetype, strtrim(dataType))
                continue
            end

            if isnan(origin)
                origin = neworigin;
            end

            if nargin > 1
                SI_factor = ImpUnits(typeUnits, newunits);
                units = newunits;
            else
                SI_factor = 1.0;
                units = typeUnits;
            end

            ncomp = 4;
            comp = ['PTXX'; 'PTXY'; 'PTYX'; 'PTYY'];
            info{nDataTypes}.type = 'Phase_Tensor';
            info{nDataTypes}.data = nan(nSites, nTx, ncomp);
            info{nDataTypes}.err = nan(nSites, nTx, ncomp);
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f');
            periods = unique(data{1});
            codes = sortrows(strtrim(char(unique(data{2}))));
            for j = 1:length(periods)
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    irx = strmatch(strtrim(codes(k,:)), codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp
                        icomp = strmatch(strtrim(comp(i,:)), complist);
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
            info{nDataTypes}.code = codes;
            info{nDataTypes}.per = periods;
            info{nDataTypes}.ncomp = 4;
            info{nDataTypes}.comp = comp;
            nDataTypes = nDataTypes + 1;
        otherwise
            warning("This datatype (%s) is unkown to 'readApres_3D'. Skipping", dataType);
            continue;
    end   
end
status = fclose(fid);


if exist('info', 'var') == 0
    error("The requested onetype '%s' was not in the file '%s'", onetype, cfile);
end

% if two data types are present, merge the periods
if length(info) == 1
    per = info{1}.per;
    ind1 = 1:length(per);
else
    per = info{1}.per;
    for i = 2:length(info)
        per = union(per, info{1}.per);
    end
    [tmp,ind1] = intersect(per,info{1}.per);
    [tmp,ind2] = intersect(per,info{2}.per);
end

% compute the maximum total number of components
ncomp = 0;
for i = 1:length(info)
    ncomp = ncomp + info{i}.ncomp;
end

% for compatibility, convert to the (old) allData structure - most programs
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
    allData{j}.Cmplx = 0;
    allData{j}.units = units;
    allData{j}.signConvention = isign;
    allData{j}.nComp = ncomp;
    allData{j}.siteLoc = allsitesloc;

    for isite = 1:size(allsites,1)
        allData{j}.siteChar{isite} = allsites(isite,:);
    end

    allData{j}.Z(1:nsites,1:ncomp) = NaN;
    allData{j}.Zerr(1:nsites,1:ncomp) = NaN;
    allData{j}.origin = [origin' 0];
    allData{j}.orient = orientation;
    allData{j}.lat = allsiteslat';
    allData{j}.lon = allsiteslon';
    allData{j}.compChar = createArray([ncomp, compCharMaxSize], 'char');
    allData{j}.compChar(:,:,:) = ' ';
    allData{j}.type = strings(0);

    icomp1 = 1;
    for k = 1:length(info)
        if itx(k) > 0
            icomp2 = icomp1 + info{k}.ncomp - 1;                
            [sites,irx] = intersect(allsites,info{k}.code,'rows');
            allData{j}.Z(irx,icomp1:icomp2) = SI_factor*squeeze(info{k}.data(:,itx(k),:));
            allData{j}.Zerr(irx,icomp1:icomp2) = SI_factor*squeeze(info{k}.err(:,itx(k),:));
            allData{j}.compChar(icomp1:icomp2,1:size(info{k}.comp,2)) = info{k}.comp;
            icomp1 = icomp2+1;
            allData{j}.type = [allData{j}.type, info{k}.type];
        end
    end
end