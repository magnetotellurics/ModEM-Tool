function [allData,header,units,isign,origin,info] = readZ_3D(cfile,newunits,onetype,nanvalue,newformat)
%  Usage:  [allData,header,units,isign,origin,info] = readZ_3D(cfile,newunits,onetype,nanvalue,newformat);
%   read contents of cell array allData from file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   Site locations units have to match the model, i.e. use meters.
%   If onetype is specified, only read that data type and skip others.
%   If nanvalue is specified, use this value, if found in place of the data 
%      entry or standard deviation, to skip a line (not typical).
%      The comparison isn't exact but works.
%   If newformat is specified, then read with the azimuth angles.
%  (c) Anna Kelbert, 2011-2013, 2020, 2023


fid = fopen(cfile,'r');

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
    %tmp = textscan(fid,'# %s',1,'delimiter','\n');
    %tmp = char(tmp{1});
    %if isempty(tmp); break; end
    %header = tmp;
    %% read the block header
    %tmp = textscan(fid,'# %s',1,'delimiter','\n');
    %% IMPORTANT: NEW LINE IN THE DATA FILE DEFINES THE TRANSMITTER TYPE
    %% (if missing, just skips over it - can still read old files)
    %tmp = textscan(fid,'+ %s',1,'delimiter','\n');
    % END IMPORTANT: NEW LINE IN THE DATA FILE
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
        case 'Full_Impedance'
            if ~isempty(onetype) && ~strcmp(onetype,strtrim(dataType))
                % just keep reading... won't use
            end
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
            info{1}.data = nan(nSites,nTx,ncomp)+1i*nan(nSites,nTx,ncomp);
            info{1}.err = nan(nSites,nTx,ncomp);
            info{1}.type = 'Full_Impedance';
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
                                    info{1}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                                    info{1}.err(k,j,i) = data{11}(ind);
                                    info{1}.lat(k) = data{3}(ind);
                                    info{1}.lon(k) = data{4}(ind);
                                    info{1}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                                catch
                                    tmp = char(data{2}(ind));
                                    warning(['Problem from line # ' num2str(ind(1)) ' site ' tmp(1,:) ' - possibly a duplicate site, using occurence 1']);
                                    info{1}.data(k,j,i) = data{9}(ind(1)) + 1i* data{10}(ind(1));
                                    info{1}.err(k,j,i) = data{11}(ind(1));
                                    info{1}.lat(k) = data{3}(ind(1));
                                    info{1}.lon(k) = data{4}(ind(1));
                                    info{1}.loc(k,1:3) = [data{5}(ind(1)) data{6}(ind(1)) data{7}(ind(1))];
                                end
                            else
                                disp(['Skipping component ',strtrim(comp(i,:)),' for period ',num2str(j),' of ',codes(k,:)]);
                            end
                        end
                    end
                end
            end
            info{1}.code = codes;
            info{1}.per = periods;
            info{1}.ncomp = 8;
            info{1}.comp = comp;
        case 'Off_Diagonal_Impedance'
            if ~isempty(onetype) && ~strcmp(onetype,strtrim(dataType))
                % just keep reading... won't use
            end
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
            info{1}.data = nan(nSites,nTx,ncomp);
            info{1}.err = nan(nSites,nTx,ncomp);
            info{1}.type = 'Off_Diagonal_Impedance';
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
                            info{1}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                            info{1}.err(k,j,i) = data{11}(ind);
                            info{1}.lat(k) = data{3}(ind);
                            info{1}.lon(k) = data{4}(ind);
                            info{1}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                        end
                    end
                end
            end
            info{1}.code = codes;
            info{1}.per = periods;
            info{1}.ncomp = 4;
            info{1}.comp = comp;
            %info{1}.ncomp = 8;
            %info{1}.comp = ['ZXX';'ZXY';'ZYX';'ZYY'];
        case 'Off_Diagonal_Rho_Phase'
            warning("Will convert 'Off_Diagonal_Rho_Phase' into Full_Impedance. Use 'readApres_3D' to read 'Off_Diagonal_Rho_Phase' without conversion")
            if ~isempty(onetype) && ~strcmp(onetype,strtrim(dataType))
                % just keep reading... won't use
            end
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
            info{1}.data = nan(nSites,nTx,ncomp);
            info{1}.err = nan(nSites,nTx,ncomp);
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
                            info{1}.data(k,j,i) = data{9}(ind);
                            info{1}.err(k,j,i) = data{10}(ind);
                            info{1}.lat(k) = data{3}(ind);
                            info{1}.lon(k) = data{4}(ind);
                            info{1}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                        end
                    end
                end
            end
            % convert to Off_Diagonal_Impedance!!! (need a separate
            % function to read Apres+Phase and a separate class to store)
            for k = 1:length(codes)
                apres.xy = info{1}.data(k,:,1);
                apres.xy_se = info{1}.err(k,:,1);
                phase.xy = info{1}.data(k,:,2);
                phase.xy_se = info{1}.err(k,:,2);
                apres.yx = info{1}.data(k,:,3);
                apres.yx_se = info{1}.err(k,:,3);
                phase.yx = info{1}.data(k,:,4);
                phase.yx_se = info{1}.err(k,:,4);
                [Z,Zstd] = mttf.apres2imp(periods,apres,phase);
                info{1}.data(k,:,1) = Z(1,1,:);
                info{1}.data(k,:,2) = Z(1,2,:);
                info{1}.data(k,:,3) = Z(2,1,:);
                info{1}.data(k,:,4) = Z(2,2,:);
                info{1}.err(k,:,1) = Zstd(1,1,:);
                info{1}.err(k,:,2) = Zstd(1,2,:);
                info{1}.err(k,:,3) = Zstd(2,1,:);
                info{1}.err(k,:,4) = Zstd(2,2,:);
            end
            info{1}.type  = 'Full_Impedance';
            info{1}.code = codes;
            info{1}.per = periods;
            info{1}.ncomp = 8;
            info{1}.comp = ['ZXX';'ZXY';'ZYX';'ZYY'];
        case 'Full_Vertical_Components'
            if ~isempty(onetype) && ~strcmp(onetype,strtrim(dataType))
                % just keep reading... won't use
            end
            if isnan(origin)
                origin = neworigin; % prevents problem when origins differ
            end
            ncomp = 2;
            comp = ['TX ';'TY '];
            info{2}.type  = 'Full_Vertical_Components';
            info{2}.data = nan(nSites,nTx,ncomp);
            info{2}.err = nan(nSites,nTx,ncomp);
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
                                info{2}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                                info{2}.err(k,j,i) = data{11}(ind);
                                info{2}.lat(k) = data{3}(ind);
                                info{2}.lon(k) = data{4}(ind);
                                info{2}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                            catch
                                tmp = char(data{2}(ind));
                                warning(['Problem from line # ' num2str(ind(1)) ' site ' tmp(1,:) ' - possibly a duplicate site, using occurence 1']);
                                info{2}.data(k,j,i) = data{9}(ind(1)) + 1i* data{10}(ind(1));
                                info{2}.err(k,j,i) = data{11}(ind(1));
                                info{2}.lat(k) = data{3}(ind(1));
                                info{2}.lon(k) = data{4}(ind(1));
                                info{2}.loc(k,1:3) = [data{5}(ind(1)) data{6}(ind(1)) data{7}(ind(1))];
                            end
                        end
                    end
                end
            end
            info{2}.code = codes;
            info{2}.per = periods;
            info{2}.ncomp = 4;
            info{2}.comp = comp;
            if ~isempty(onetype) && strcmp(onetype,strtrim(dataType))
                tmp = info{2}; clear info; info{1} = tmp; units = []; SI_factor = 1;
            elseif ~isempty(onetype)
                tmp = info{1}; clear info; info{1} = tmp; % a clutch to fix later    
            end
        otherwise
            disp('Unknown data type');
            break;
    end   
end
status = fclose(fid);

% if two data types are present, merge the periods
if length(info) == 1
    per = info{1}.per;
    ind1 = 1:length(per);
else
    per = sort(union(info{1}.per,info{2}.per));
    [tmp,ind1] = intersect(per,info{1}.per);
    [tmp,ind2] = intersect(per,info{2}.per);
end

% compute the maximum total number of components
ncomp = 0;
for i = 1:length(info)
    ncomp = ncomp + info{i}.ncomp/2;
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
    allData{j}.Cmplx = 1;
    allData{j}.units = units;
    allData{j}.signConvention = isign;
    allData{j}.nComp = ncomp;
    %allData{j}.compChar(1:ncomp) = '';
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
            icomp1 = icomp2+1;
        end
    end
end