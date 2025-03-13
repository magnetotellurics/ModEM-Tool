function [header,units,isign,origin,info] = readZ_3Dinfo(cfile,newunits)
%  Usage:  [allData,header,units,isign,origin,info] = readZ_3D(cfile,newunits);
%   read contents of cell array allData from file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   Site locations units have to match the model, i.e. use meters.
%  (c) Anna Kelbert, 2011


fid = fopen(cfile,'r');

%  read the data: one block per data type
while 1
    % read the header info
    tmp = textscan(fid,'# %s',1,'delimiter','\n');
    tmp = char(tmp{1});
    if isempty(tmp); break; end
    header = tmp;
    % read the block header
    tmp = textscan(fid,'# %s',1,'delimiter','\n');
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
    origin = sscanf(blockinfo(5,:),'%f',2);
    tmp  = sscanf(blockinfo(6,:),'%d',2);
    nTx = tmp(1);
    nSites = tmp(2);
    % now read the data line by line
    switch strtrim(dataType)
        case 'Full_Impedance'
            if nargin > 1
                SI_factor = ImpUnits(typeUnits,newunits);
                units = newunits;
            else
                SI_factor = 1.0;
                units = typeUnits; % output the units of impedances
            end
            ncomp = 4;
            comp = ['ZXX';'ZXY';'ZYX';'ZYY'];
            info{1}.data = nan(nSites,nTx,ncomp);
            info{1}.err = nan(nSites,nTx,ncomp);
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            periods = unique(data{1});
            codes = sortrows(char(unique(data{2})));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(codes(k,:),char(data{2}(itx)));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(comp(i,:),char(data{8}(itx(irx))));
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
            info{1}.comp = comp;
        case 'Full_Vertical_Components'
            ncomp = 2;
            comp = ['TX ';'TY '];
            info{2}.data = nan(nSites,nTx,ncomp);
            info{2}.err = nan(nSites,nTx,ncomp);
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            periods = unique(data{1});
            codes = sortrows(char(unique(data{2})));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(codes(k,:),char(data{2}(itx)));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),strtrim(char(data{8}(itx(irx)))));
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            info{2}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                            info{2}.err(k,j,i) = data{11}(ind);
                            info{2}.lat(k) = data{3}(ind);
                            info{2}.lon(k) = data{4}(ind);
                            info{2}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                        end
                    end
                end
            end
            info{2}.code = codes;
            info{2}.per = periods;
            info{2}.comp = comp;
        otherwise
            disp('Unknown data type');
            break;
    end   
end
status = fclose(fid);

% for compatibility, convert to the (old) allData structure - assume full
% impedances and vertical components (in this order) with same periods etc
%for j = 1:length(info{1}.per)
%    irx = find(~isnan(info{1}.data(:,j,1)));
%    allData{j}.T = info{1}.per(j);
%    allData{j}.Cmplx = 1;
%    allData{j}.units = units;
%    allData{j}.signConvention = isign;
%    allData{j}.nComp = 8;
%    allData{j}.compChar = info{1}.comp;
%    allData{j}.siteLoc = info{1}.loc(irx,:);
%    allData{j}.siteChar = info{1}.code(irx,:);
%    allData{j}.Z = SI_factor*squeeze(info{1}.data(irx,j,:));
%    allData{j}.Zerr = SI_factor*squeeze(info{1}.err(irx,j,:));
%    if (length(info) > 1)
%        allData{j}.nComp = 12;
%        allData{j}.compChar(5:6,:) = info{2}.comp;
%        allData{j}.Z(:,5:6) = squeeze(info{2}.data(irx,j,:));
%        allData{j}.Zerr(:,5:6) = squeeze(info{2}.err(irx,j,:));
%    end
%    allData{j}.origin = [origin' 0];
%    allData{j}.orient = orientation;
%    allData{j}.lat = info{1}.lat(irx)';
%    allData{j}.lon = info{1}.lon(irx)';
%end

