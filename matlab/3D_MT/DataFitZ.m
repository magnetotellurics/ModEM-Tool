function [rms,info,site] = DataFitZ(resp,data,comp,maxrms,fplot,sitenames,iper)

% [rms,info,site] = DataFitZ(resp,data,comp,maxrms,fplot,sitenames,iper)
%
% Reads two data files and outputs a pseudo-plot for the responses
% and the data misfit for one component.

%% Read real and synthetic data
%resp = readZ_3D(fresp,'[mV/km]/[nT]');
%data = readZ_3D(fdata,'[mV/km]/[nT]');

PLOT = 1;
if nargin < 5
    [pathstr, name] = fileparts(fresp);
    fplot = name;
end

nper = length(data);
ncomp = 6;
[tmp,icomp]=intersect(resp{1}.compChar,comp,'rows')
if isempty(icomp)
    disp(['Unknown component ' comp]);
    return
end
if strcmp(comp,'ZXY')
    clims = [0 30];
elseif strcmp(comp,'ZYX')
    clims = [-30 0];
else
    clims = [-15 15];
end

%% Make a list of all sites in the real data file
allsites = data{1}.siteChar;
allsitesloc = data{1}.siteLoc;
for j = 2:nper
    for i = 1:length(data{j}.siteChar)
        if isempty(intersect(allsites,data{j}.siteChar(i,:),'rows'))
            allsites = [allsites; data{j}.siteChar(i,:)];
            allsitesloc = [allsitesloc; data{j}.siteLoc(i,:)];
        end
    end
end
nsites = length(allsites);
x = allsitesloc(:,1);
y = allsitesloc(:,2);
[lat,lon] = xy2latlon(x,y,data{1}.origin(1),data{1}.origin(2),'m');

%% Regroup into observatory bins
info.data = nan(length(allsites),nper,ncomp)+1i*nan(length(allsites),nper,ncomp);
info.resp = nan(length(allsites),nper,ncomp)+1i*nan(length(allsites),nper,ncomp);
info.err = nan(length(allsites),nper,ncomp);
for j = 1:nper
    for i = 1:length(allsites)
        if nargin > 4
            in = 0;
            [dummy,k]=intersect(sitenames,allsites(i,:),'rows');
            if ~isempty(k)
                in = 1;
            end
        else
            in = 1;
        end
        [dummy,k]=intersect(data{j}.siteChar,allsites(i,:),'rows');
        if ~isempty(k) && in
            info.code(i,:) = allsites(i,:);
            info.loc(i,:) = allsitesloc(i,:);
            info.lon(i) = lon(i);
            info.lat(i) = lat(i);
            info.data(i,j,:) = data{j}.Z(k,:);
            info.err(i,j,:) = data{j}.Zerr(k,:);
            info.per(j) = data{j}.T;
        end
        [dummy,k]=intersect(resp{j}.siteChar,allsites(i,:),'rows');
        if ~isempty(k) && in
            info.resp(i,j,:) = resp{j}.Z(k,:);
        end
    end
end
info.res = 0.5*(real(info.data - info.resp)./info.err).^2 + ...
    0.5*(imag(info.data - info.resp)./info.err).^2;

%% Compute total RMS
count = 0;
misfit = 0;
for j = 1:nper
    for i = 1:length(allsites)
        count = count + ncomp - sum(isnan(info.err(i,j,:)));
        misfit = misfit + nansum(info.res(i,j,:));
    end
end
rms = sqrt(misfit/count);

if ~PLOT
    return
end

%% Average over frequency and component bins
for j = 1:nper
    freq.f(j) = 1/(info.per(j));
    freq.p(j) = info.per(j);
    freq.d(j) = info.per(j)/(24*3600);
    for k = 1:ncomp
        count = nsites - sum(isnan(info.err(:,j,k)));
        freq.rms(j,k) = sqrt(nansum(info.res(:,j,k))/count);
    end
end

%% Average over site bins
j = 1;
for i = 1:length(allsites)
    count = nper - sum(isnan(info.res(i,:,icomp)));
    if count > 0
        site.codes(j,:) = allsites(i,:);
        site.lon(j) = info.lon(i);
        site.lat(j) = info.lat(i);
        site.resp(j,:) = info.resp(i,:,icomp);
        site.rms(j,:) = sqrt(info.res(i,:,icomp));
        j = j+1;
    end
end    

%% Plot
if nargin < 7
    iper = 1:length(freq.p);
end
for k=iper
cper = sprintf('%03f',freq.p(k));
f=figure('Position',[300,300,1400,600],...
        'PaperPosition',[1,1,14,6],...
        'PaperOrientation','Portrait');
%figure(1), clf
%subplot('position',[0.1 0.1 0.6 0.8]);
%cmax = nanmax(nanmax(site.resp(:,k)));
%cmin = nanmin(nanmin(site.resp(:,k)));
subplot(2,3,1);
htitle = [num2str(freq.p(k)) ' secs; total RMS ' num2str(rms)];
InterpPlot(site.lon,site.lat,sqrt(freq.p(k))*real(site.resp(:,k)),['Re ' comp ' x sqrt(T)'],'',clims,5);
colormap(jet);
subplot(2,3,4);
InterpPlot(site.lon,site.lat,sqrt(freq.p(k))*imag(site.resp(:,k)),['Im ' comp ' x sqrt(T)'],'',clims,5);
colormap(jet);
subplot(2,3,[2 3 5 6]);
InterpPlot(site.lon,site.lat,site.rms(:,k),'RMS',htitle,[1 5],5);
%colormap(invhot);
if nargin > 2
    caxis([1 maxrms]);
end
% xlabel('Log(10) T (secs)');
% ylabel('RMS');
% legend(h,'Zxx','Zxy','Zyx','Zyy','Tx','Ty','Location','BestOutside');
print('-djpeg95','-r300',[fplot '_' comp '_' cper '.jpg']);
close(f);
end
