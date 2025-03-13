function [rms,info] = DataFit(fresp,fdata,maxrms,fplot,sitenames,prange)

% [rms,info] = DataFit(fresp,fdata,maxrms,fplot,sitenames,prange)
%
% Reads two data files and outputs the RMS misfit. Plots a summary.

%% Read real and synthetic data
resp = readZ_3D(fresp,'[mV/km]/[nT]');
data = readZ_3D(fdata,'[mV/km]/[nT]');

PLOT = 1;
if nargin < 4
    [pathstr, name] = fileparts(fresp);
    fplot = name;
end

nper = length(data);
ncomp = data{1}.nComp/2;

if nargin > 5
    minper = prange(1);
    maxper = prange(2);
else
    minper = 1e-4;
    maxper = 1e6;
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
    if info.per(j) < maxper && info.per(j) > minper
        for i = 1:length(allsites)
            count = count + ncomp - sum(isnan(info.err(i,j,:)));
            misfit = misfit + nansum(info.res(i,j,:));
        end
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
    count = ncomp*nper - sum(sum(isnan(info.res(i,:,:))));
    if count > 0
       site.codes(j,:) = allsites(i,:);
       site.lon(j) = info.lon(i);
       site.lat(j) = info.lat(i);
       site.rms(j) = sqrt(nansum(nansum(info.res(i,:,:)))/count);
       j = j+1;
    end
end    

%% Plot
figure('Position',[300,300,600,600],...
        'PaperPosition',[1,1,6,6],...
        'PaperOrientation','Portrait');
%figure(1), clf
%subplot('position',[0.1 0.1 0.6 0.8]);
subplot(7,1,1:5);
htitle = ['RMS ' num2str(rms) ' (averaged over ' num2str(nper) ' periods)'];
InterpPlot(site.lon,site.lat,site.rms,'RMS',htitle);
colormap(invhot);
if nargin > 2
    caxis([1 maxrms]);
end
%proj_hammer_sp(clong,'','patch',P);
%color_ind=ceil(min(site.rms,2.5)/(2.5)*63)+1;
%[x,y]= m_ll2xy(lon,lat);
% for i=1:length(site.rms)
%     if ~isnan(site.rms(i))
%         p=m_plot(lon(i),lat(i),'wo',...
%             'color',P(color_ind(i),:),...
%             'markerfacecolor',P(color_ind(i),:),...
%             'markeredgecolor','k',...
%             'linewidth',0.5);
%         set(p,'markersize',5);
%     end
% end
% cb = colorbar('location','southoutside');
% set(cb,'xtick',[0 0.25 0.5 0.75 1],'xticklabel',...
%     [0 0.5 1 1.5 2.5])

%subplot('position',[0.7 0.1 0.3 0.8]);
subplot(7,1,6:7);
h = zeros(1,6);
for i = 1:ncomp
    color = [(i-1)/ncomp (i-1)/ncomp (i-1)/ncomp];
    h(i)=semilogx(freq.p,freq.rms(:,i),'x-','color',color); hold on;
    set(h(i),'LineWidth',2);
end
xlabel('Log(10) T (secs)');
ylabel('RMS');
legend(h,'Zxx','Zxy','Zyx','Zyy','Tx','Ty','Location','BestOutside');
print('-djpeg95','-r300',[fplot '_rms.jpg']);
