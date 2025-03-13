function [status] = writeZ_3D(cfile,allData,header,units,isign,APRES)
%  Usage:  [status] = writeZ_3D(cfile,allData,header,units,isign,APRES);
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
    APRES = 0;
end

%  assume 2, 4 or 6 complex components in this order (to match the other scripts)
ncomp = allData{1}.nComp/2;
if ncomp == 2
    %comp = ['ZXY';'ZYX'];
    %dataType = 'Off_Diagonal_Impedance';
    comp = ['TX';'TY'];
    dataType = 'Full_Vertical_Components';
else
    comp = ['ZXX';'ZXY';'ZYX';'ZYY';'TX ';'TY '];
    dataType = 'Full_Impedance';
end

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
                info.err(i,j,:) = allData{j}.nanvalue;
            end
            info.per(j) = allData{j}.T;
        end
    end
end
info.comp = comp;

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
        units = '[V/m]/[T]'
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
% Z(1)  = dlog(abs(tempZ(1))**2*MU_0/omega)
% Z(2)  = atan2(ISIGN*dimag(tempZ(1)),real(tempZ(1)))*R2D
% Z(3)  = dlog(abs(tempZ(2))**2*MU_0/omega)
% Z(4)  = atan2(ISIGN*dimag(tempZ(2)),real(tempZ(2)))*R2D+180.0d0
if APRES
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

%  file header: the info line
fprintf(fid,'# %s\n',header);

%  data type header: comment, then six lines
comment = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error';
fprintf(fid,'# %s\n',comment);

% IMPORTANT: NEW LINE IN THE DATA FILE DEFINES THE TRANSMITTER TYPE
% (this is omitted when txType = '' for backward compatibility)
if ~isempty(txType)
    fprintf(fid,'+ %s\n',txType);
end
% END IMPORTANT: NEW LINE IN THE DATA FILE

fprintf(fid,'> %s\n',dataType);
fprintf(fid,'> %s\n',signstr);
fprintf(fid,'> %s\n',units);
fprintf(fid,'> %.2f\n',orientation);
fprintf(fid,'> %.3f %.3f\n',origin(1:2));
fprintf(fid,'> %d %d\n',nTx,nSites);

%  now write all impedances (or apparent resistivities!) line by line
if APRES
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
else
    for k = 1:nSites
        for j = 1:nTx
            for i = 1:min(ncomp,4)
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

if ncomp > 4

    %  repeat file header
    fprintf(fid,'# %s\n',header);

    %  data type header: comment, then six lines
    comment = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error';
    fprintf(fid,'# %s\n',comment);
    fprintf(fid,'> %s\n','Full_Vertical_Components');
    fprintf(fid,'> %s\n',signstr);
    fprintf(fid,'> %s\n','[]');
    fprintf(fid,'> %.2f\n',orientation);
    fprintf(fid,'> %.3f %.3f\n',origin(1:2));
    fprintf(fid,'> %d %d\n',nTx,nSites);

    %  now write all impedances line by line
    for k = 1:nSites
        for j = 1:nTx
            for i = 5:6
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

status = fclose(fid);
