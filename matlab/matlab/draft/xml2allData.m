% load all transfer function files
% list = dir('*.xml');
% nsites = length(list);
% files = { list.name };
files = findfiles('xml','MOCHA');
nsites = length(files);

T = zeros(200,1);
lat = zeros(nsites,1);
lon = zeros(nsites,1);
elev = zeros(nsites,1);
periods = struct([]);    
siteID = struct([]);
TF = struct([]);
TFVar = struct([]);

np = 0;

for i=1:nsites

    [pathstr,name,ext] = fileparts(files{i});
    disp(['Reading XML file ' name ext]);
    Data = xml_parseany(fileread(files{i}));
    
    siteID{i} = Data.SiteID{1}.CONTENT;
    lat(i) = str2double(Data.Location{1}.Latitude{1}.CONTENT);
    lon(i) = str2double(Data.Location{1}.Longitude{1}.CONTENT);
    elev(i) = str2double(Data.Location{1}.Elevation{1}.CONTENT);

    nperiods = str2double(Data.Frequencies{1}.ATTRIBUTE.count);
    
    for j=1:nperiods
        % save periods in seconds
        type = Data.Frequencies{1}.Frequency{j}.ATTRIBUTE.type;
        value = str2double(Data.Frequencies{1}.Frequency{j}.ATTRIBUTE.value);
        if strcmp(type,'period')
            periods{i}(j) = value;
        else
            periods{i}(j) = 1/value;
        end
        % relate to the common (for all sites) list of periods
        ip = find(T==periods{i}(j), 1);
        if isempty(ip)
           np = np + 1;
           T(np) = periods{i}(j);
           ip = np;
        end
        % save the transfer functions and variances
        TF{i,ip} = Data.Frequencies{1}.Frequency{j}.Z{1}.component;
        TFVar{i,ip} = Data.Frequencies{1}.Frequency{j}.Z.VAR{1}.value;
    end
    
end

nfreq = np;
T = T(1:np);

% store data in a 3D array
Zchar = ['Zxx';'Zxy';'Zyx';'Zyy';'Tx ';'Ty '];
Z = NaN * zeros(nsites,nfreq,6) + 1i*NaN * zeros(nsites,nfreq,6);
Zstd = NaN * zeros(nsites,nfreq,6);

for i = 1:nsites
    for j = 1:nfreq
        % transfer functions (values)
        for k = 1:length(TF{i,j})
            chanin = TF{i,j}{k}.ATTRIBUTE.input;
            chanout = TF{i,j}{k}.ATTRIBUTE.output;
            zreal = str2double(TF{i,j}{k}.real{1}.CONTENT);
            zimag = str2double(TF{i,j}{k}.imag{1}.CONTENT);
            switch [chanout '/' chanin]
                case 'Ex/Hx' % Zxx
                    Z(i,j,1) = zreal + 1i*zimag;
                case 'Ex/Hy' % Zxy
                    Z(i,j,2) = zreal + 1i*zimag;
                case 'Ey/Hx' % Zyx
                    Z(i,j,3) = zreal + 1i*zimag;
                case 'Ey/Hy' % Zyy
                    Z(i,j,4) = zreal + 1i*zimag;
                case 'Hz/Hx' % Tx
                    Z(i,j,5) = zreal + 1i*zimag;
                case 'Hz/Hy' % Ty
                    Z(i,j,6) = zreal + 1i*zimag;
            end
            % transfer functions (standard errors)
            chanin = TFVar{i,j}{k}.ATTRIBUTE.input;
            chanout = TFVar{i,j}{k}.ATTRIBUTE.output;
            zvar = str2double(TFVar{i,j}{k}.CONTENT);
            switch [chanout '/' chanin]
                case 'Ex/Hx' % Zxx
                    Zstd(i,j,1) = sqrt(zvar);
                case 'Ex/Hy' % Zxy
                    Zstd(i,j,2) = sqrt(zvar);
                case 'Ey/Hx' % Zyx
                    Zstd(i,j,3) = sqrt(zvar);
                case 'Ey/Hy' % Zyy
                    Zstd(i,j,4) = sqrt(zvar);
                case 'Hz/Hx' % Tx
                    Zstd(i,j,5) = sqrt(zvar);
                case 'Hz/Hy' % Ty
                    Zstd(i,j,6) = sqrt(zvar);
            end
        end % loop over transfer functions    
    end % frequencies
end % sites

% sort by period (increasing)
[sorted,ind] = sort(T);
T = sorted;
Z = Z(:,ind,:);
Zstd = Zstd(:,ind,:);

% sort approximately by latitude, then longitude
[sorted,ind] = sortrows([round(lat) lon],[1 2]);
sites = char(siteID{ind});
lat = lat(ind);
lon = lon(ind);
elev = elev(ind);
Z = Z(ind,:,:);
Zstd = Zstd(ind,:,:);

% delete duplicates, if present
% dupl = [0; ((abs(diff(lat)) < 1e-3) & (abs(diff(lon)) < 1e-3)); 0];
% ind = find(dupl==0);
% sites = sites(ind,:);
% lat = lat(ind);
% lon = lon(ind);
% elev = elev(ind);
% Z = Z(ind,:,:);
% Zstd = Zstd(ind,:,:);

% interpolate (replace NaNs with valid data, where possible)
for i = 1:nsites
    for k = 1:6
        temp = Z(i,:,k);
        ind = find(isnan(temp)==0);
        Z(i,:,k) = interp1(T(ind),temp(ind),T);
        temp = Zstd(i,:,k);
        ind = find(isnan(temp)==0);
        Zstd(i,:,k) = interp1(T(ind),temp(ind),T);        
    end
end

% make sure NaNs are still complex
ind = find(isnan(Z)==1);
Z(ind) = NaN + 1i*NaN;

% Save into an allData structure, that, historically, is a frequency
% array, each element looking like:
%            T: 102.4000
%        Cmplx: 1
%        units: 'Ohm'
%        nComp: 8
%     compChar: [8x7 char]
%      siteLoc: [109x3 double]
%     siteChar: [109x3 char]
%            Z: [109x4 double]
%         Zerr: [109x4 double]
% This can be written to an input data file for ModEM with
% writeZ_3D('file.imp',allData,'Data from XML files','[mV/km]/[nT]',+1)

allData = struct([]);
compChar = ['Re(Zxx)'; 'Im(Zxx)'; ...
    'Re(Zxy)'; 'Im(Zxy)'; 'Re(Zyx)'; 'Im(Zyx)'; ...
    'Re(Zyy)'; 'Im(Zyy)'; ...
    'Re(Tx) '; 'Im(Tx) '; ...
    'Re(Ty) '; 'Im(Ty) '];
nComp = length(compChar);

% Stick to the origin we have chosen previously for compatibility
origin = [42.015999 -116.476997 0];

% Compute coordinates in metres
%origin = [min(lat) min(lon) 0];
[lat_m,lon_m] = latlon2xy(lat,lon,origin(1),origin(2));

% Use Mapping Toolbox to compute coordinates in km
% origin = [min(lat) min(lon) 0];
% lat_km = distance(min(lat),min(lon),lat,min(lon),[6371 0]);
% lon_km = distance(min(lat),min(lon),min(lat),lon,[6371 0]);
% lat_m = 1000*lat_km;
% lon_m = 1000*lon_km;

for j = 1:nfreq
    allData{j} = struct('T',T(j),...
        'Cmplx',1,'units','[mV/km]/[nT]',...
        'nComp',nComp,'compChar',compChar,...
        'siteLoc',[lat_m lon_m -elev],'siteChar',char(sites),...
        'Z',squeeze(Z(:,j,:)),'Zerr',squeeze(Zstd(:,j,:)),...
        'origin',origin,'lat',lat,'lon',lon);
end
%%
% % now, create a subset of 6-8 frequencies for inversion ...
% % subset = [4 8 14 18 24 30 39];
% % subset = 1:3:42;
% subset = 3:3:42;
% allDataSubset = struct([]);
% for j = 1:length(subset)
%     allDataSubset{j} = allData{subset(j)};
% end
%writeZ_3D('YellowstoneWithTopography_42freq.imp',allData,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, original errors','[mV/km]/[nT]',+1)
%writeZ_3D('YellowstoneWithTopography_14freq.imp',allDataSubset,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, original errors','[mV/km]/[nT]',+1)

% also, create a data set & subset without topography
nsites = size(allData{1}.siteLoc,1);
siteLoc = allData{1}.siteLoc;
siteLoc = [siteLoc(:,1) siteLoc(:,2) zeros(nsites,1)];
allDataNoTopography = allData;
for j = 1:nfreq
    allDataNoTopography{j}.siteLoc = siteLoc;
end
% allDataSubsetNoTopography = struct([]);
% for j = 1:length(subset)
%     allDataSubsetNoTopography{j} = allDataNoTopography{subset(j)};
% end
writeZ_3D('Yellowstone_42freq_allData.dat',allDataNoTopography,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, original errors, no topography','[mV/km]/[nT]',+1)
%writeZ_3D('Yellowstone_14freq_allData.dat',allDataSubsetNoTopography,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, original errors, no topography','[mV/km]/[nT]',+1)

%%
% allDataErrFloor3 = setErrorFloor3D(allData,0.03);
% allDataErrFloor5 = setErrorFloor3D(allData,0.05);
% %writeZ_3D('YellowstoneWithTopography_ErrFl3_42freq.imp',allDataErrFloor3,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, error floor 3%','[mV/km]/[nT]',+1)
% %writeZ_3D('YellowstoneWithTopography_ErrFl5_42freq.imp',allDataErrFloor5,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, error floor 5%','[mV/km]/[nT]',+1)
% allDataNTErrFloor3 = setErrorFloor3D(allDataNoTopography,0.03);
% allDataNTErrFloor5 = setErrorFloor3D(allDataNoTopography,0.05);
% writeZ_3D('Yellowstone_ErrFl3_42freq.imp',allDataNTErrFloor3,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, error floor 3%, no topography','[mV/km]/[nT]',+1)
% writeZ_3D('Yellowstone_ErrFl5_42freq.imp',allDataNTErrFloor5,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, error floor 5%, no topography','[mV/km]/[nT]',+1)
% allDataSubsetErrFloor3 = setErrorFloor3D(allDataSubset,0.03);
% allDataSubsetErrFloor5 = setErrorFloor3D(allDataSubset,0.05);
% %writeZ_3D('YellowstoneWithTopography_ErrFl3_14freq.imp',allDataSubsetErrFloor3,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, error floor 3%','[mV/km]/[nT]',+1)
% %writeZ_3D('YellowstoneWithTopography_ErrFl5_14freq.imp',allDataSubsetErrFloor5,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, error floor 5%','[mV/km]/[nT]',+1)
% allDataSubsetNTErrFloor3 = setErrorFloor3D(allDataSubsetNoTopography,0.03);
% allDataSubsetNTErrFloor5 = setErrorFloor3D(allDataSubsetNoTopography,0.05);
% writeZ_3D('Yellowstone_ErrFl3_14freq.imp',allDataSubsetNTErrFloor3,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, error floor 3%, no topography','[mV/km]/[nT]',+1)
% writeZ_3D('Yellowstone_ErrFl5_14freq.imp',allDataSubsetNTErrFloor5,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, error floor 5%, no topography','[mV/km]/[nT]',+1)
% 
% allDataZonly = allDataSubsetNoTopography;
% for j = 1:length(allDataZonly)
%     allDataZonly{j}.nComp = 8;
%     allDataZonly{j}.compChar = ['Re(Zxx)'; 'Im(Zxx)'; ...
%     'Re(Zxy)'; 'Im(Zxy)'; 'Re(Zyx)'; 'Im(Zyx)'; ...
%     'Re(Zyy)'; 'Im(Zyy)'];
%     allDataZonly{j}.Z = allDataZonly{j}.Z(:,1:4);
%     allDataZonly{j}.Zerr = allDataZonly{j}.Zerr(:,1:4);
% end
% allDataZonlyErrFloor3 = setErrorFloor3D(allDataZonly,0.03);
% allDataZonlyErrFloor5 = setErrorFloor3D(allDataZonly,0.05);
% writeZ_3D('Yellowstone7NoTopographyZonly.imp',allDataZonly,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, no error floor (no tippers)','[mV/km]/[nT]',+1)
% writeZ_3D('Yellowstone7NoTopographyZonlyErrFl5.imp',allDataZonlyErrFloor5,'Earthscope 2006-2009 + Snake River Plain around Yellowstone from XML files, error floor 5% (no tippers)','[mV/km]/[nT]',+1)


% save
% save YellowstoneAllData allData

% Map from allData to 3D arrays (assumes same site locations for all freq.)
% Z = NaN * zeros(nsites,nfreq,6);
% Zstd = NaN * zeros(nsites,nfreq,6);
% for j = 1:nfreq
%     Z(:,j,:) = allData{j}.Z;
%     Zstd(:,j,:) = allData{j}.Zerr;
%     lat = allData{1}.lat;
%     lon = allData{1}.lon;
% end
% Undistort, then map back to allData
% Zobs = Z;
% [Zund,dist,Cinv1,Cinv2] = undistort(lon,lat,Zobs,Zstd,dist);
% allDataUndistort = struct([]);
% for j = 1:nfreq
%     allDataUndistort{j} = struct('T',T(j),...
%         'Cmplx',1,'units','[mV/km]/[nT]',...
%         'nComp',nComp,'compChar',compChar,...
%         'siteLoc',[lat_km lon_km elev],'siteChar',char(sites),...
%         'Z',squeeze(Zund(:,j,:)),'Zerr',squeeze(Zstd(:,j,:)),...
%         'origin',origin,'lat',lat,'lon',lon);
% end
% 
% save allData allData allDataUndistort Cinv1 Cinv2 Z Zund Zstd dist lon lat;