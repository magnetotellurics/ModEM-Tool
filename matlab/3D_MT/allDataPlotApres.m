fdir = '.';
filt = [fdir '/*.dat'];
[filename, pathname] = uigetfile(filt, 'Impedance File');
fdata = [pathname filename];
fresp = [pathname filename];
filename = filename(1:end-4);
[rms,info] = DataFit(fresp,fdata);

%% save sites list
fid = fopen('sites.txt','w');
for i = 1:length(info.code)
    fprintf(fid,'%s %f %f\n',info.code(i,:),info.lat(i),info.lon(i));
end
fclose(fid);

%% compute apparent resistivity and phase
rad_deg = 57.2958;
apres.xy = abs(info.data(:,:,2)).^2;
apres.xy_re = abs(info.resp(:,:,2)).^2;
apres.xy_se = info.err(:,:,2).^2; % variance of Z as in EDI file
apres.yx = abs(info.data(:,:,3)).^2;
apres.yx_re = abs(info.resp(:,:,3)).^2;
apres.yx_se = info.err(:,:,3).^2; % variance of Z as in EDI file
phase.xy = rad_deg*atan((imag(info.data(:,:,2)))./real(info.data(:,:,2)));
phase.xy_re = rad_deg*atan((imag(info.resp(:,:,2)))./real(info.resp(:,:,2)));
phase.xy_se = rad_deg*sqrt(apres.xy_se./apres.xy);
phase.yx = rad_deg*atan((imag(info.data(:,:,3)))./real(info.data(:,:,3)));
phase.yx_re = rad_deg*atan((imag(info.resp(:,:,3)))./real(info.resp(:,:,3)));
phase.yx_se = rad_deg*sqrt(apres.yx_se./apres.yx);

%% rescale apparent resistivity by period
for l = 1:length(info.per)
  apres.yx(:,l) = apres.yx(:,l)*info.per(l)/5. ;
  apres.xy(:,l) = apres.xy(:,l)*info.per(l)/5. ;
  apres.yx_re(:,l) = apres.yx_re(:,l)*info.per(l)/5. ;
  apres.xy_re(:,l) = apres.xy_re(:,l)*info.per(l)/5. ;
  apres.yx_se(:,l) = sqrt(apres.yx_se(:,l).*apres.yx(:,l)*info.per(l)*4/5.);
  apres.xy_se(:,l) = sqrt(apres.xy_se(:,l).*apres.xy(:,l)*info.per(l)*4/5.);
end

%% also store the vertical transfer functions
tx = info.data(:,:,5);
tx_se = info.err(:,:,5);
ty = info.data(:,:,6);
ty_se = info.err(:,:,6);

%%
% OBS = 'MTG17';
% [str,i] = intersect(info.code,OBS,'rows');
% if isempty(i)
%     error(['No observatory in data file matching ' char(OBS)]);
% end
for i=1:length(info.lon)
OBS = char(info.code(i,:));
lon = num2str(info.lon(i));
lat = num2str(info.lat(i));
per = info.per;

f=figure; clf; subplot(3,1,1)
p(1)=loglog(per,apres.xy(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
errorbar(per,apres.xy(i,:),2*apres.xy_se(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
p(2)=loglog(per,apres.yx(i,:),'rx','linewidth',1.5); hold on
errorbar(per,apres.yx(i,:),2*apres.yx_se(i,:),'rx','linewidth',1.5); hold off
set(gca,'xlim',[5 10^5],'ylim',[0.1 20000]);
legend(p,'XY','YX');
ylabel('\rho_a (\Omega m)','fontweight','demi');
RMS = sqrt(nansum(info.res(i,:))/sum(~isnan(info.err(i,:))));
title([char(OBS) ' [ LON = ' lon '; LAT = ' lat '; RMS = ' num2str(RMS)  ' ]'],'fontweight','demi');
subplot(3,1,2)
p(1)=semilogx(per,phase.xy(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
errorbar(per,phase.xy(i,:),2*phase.xy_se(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
p(2)=semilogx(per,phase.yx(i,:),'rx','linewidth',1.5); hold on
errorbar(per,phase.yx(i,:),2*phase.yx_se(i,:),'rx','linewidth',1.5); hold off
set(gca,'xlim',[5 10^5],'ylim',[0 90]);
ylabel('\phi (Degrees)','fontweight','demi');
xlabel('Period (secs)','fontweight','demi');
subplot(3,1,3)
p(1)=semilogx(per,tx(i,:),'go','color',[0. 0.5 0.8],'linewidth',1.5); hold on
errorbar(per,tx(i,:),2*tx_se(i,:),'go','color',[0. 0.5 0.8],'linewidth',1.5); hold on
p(2)=semilogx(per,ty(i,:),'mx','linewidth',1.5); hold on
errorbar(per,ty(i,:),2*ty_se(i,:),'mx','linewidth',1.5); hold off
set(gca,'xlim',[5 10^5],'ylim',[-1 1]);
legend(p,'TX','TY');
ylabel('Vertical TFs','fontweight','demi');
xlabel('Period (secs)','fontweight','demi');
print(f,'-dpng',[pathname '/data/' filename '_' OBS '_apres']);
close(f);
end
%% 
rho2 = [apres.xy(i,:); apres.yx(i,:)].';
rho2_se = 2*[apres.xy_se(i,:); apres.yx_se(i,:)].';
ph2 = [phase.xy(i,:); phase.yx(i,:)].';
ph2_se = 2*[phase.xy_se(i,:); phase.yx_se(i,:)].';
pltall = ones(length(per),1);
[rho2_axes,ph2_axes] = pltrhom(length(per),pltall,per,rho2,rho2_se,ph2,ph2_se,lims,c_title_Modes,hfig);

%%
posn = [1,1,10,10];
fg = figure('Position',100*posn,...
        'PaperPosition',posn,...
        'PaperOrientation','Portrait',...
        'Color',[1 1 1]);

for k=8%,length(AllData1)
    x = AllData1{k}.siteLoc(:,1);
    y = AllData1{k}.siteLoc(:,2);
    origin = AllData1{k}.origin;
    [lat,lon] = xy2latlon(x,y,origin(1),origin(2));
    T = AllData1{k}.T;
    Z1 = AllData1{k}.Z;
    Z2 = AllData2{k}.Z;
    % Zxx
    clf
    rmin = min(real(Z1(:,1))); rmax = max(real(Z1(:,1)));
    imin = min(imag(Z1(:,1))); imax = max(imag(Z1(:,1)));
    subplot(2,2,1);
    %axes('Position',[1 1 plotsize plotsize]);
    %InterpPlot(lon,lat,real(Z1(:,1)),'Re(Zxx)',[num2str(T) ' secs; measured'],[-20 20],6);
    InterpPlot(lon,lat,real(Z1(:,1)),'Re(Zxx)',[num2str(T) ' secs; measured'],[rmin rmax],6);
    subplot(2,2,3);
    InterpPlot(lon,lat,imag(Z1(:,1)),'Im(Zxx)',[num2str(T) ' secs; measured'],[imin imax],6);
    subplot(2,2,2);
    InterpPlot(lon,lat,real(Z2(:,1)),'Re(Zxx)',[num2str(T) ' secs; computed'],[rmin rmax],6);
    subplot(2,2,4);
    InterpPlot(lon,lat,imag(Z2(:,1)),'Im(Zxx)',[num2str(T) ' secs; computed'],[imin imax],6);
    strper = sprintf('%03d',k);
    print('-r300','-djpeg',[filename '_Zxx_' strper '.jpg']);
    % Zxy
    clf
    rmin = min(real(Z1(:,2))); rmax = max(real(Z1(:,2)));
    imin = min(imag(Z1(:,2))); imax = max(imag(Z1(:,2)));
    subplot(2,2,1);
    InterpPlot(lon,lat,real(Z1(:,2)),'Re(Zxy)',[num2str(T) ' secs; measured'],[rmin rmax],6);
    subplot(2,2,3);
    InterpPlot(lon,lat,imag(Z1(:,2)),'Im(Zxy)',[num2str(T) ' secs; measured'],[imin imax],6);
    subplot(2,2,2);
    InterpPlot(lon,lat,real(Z2(:,2)),'Re(Zxy)',[num2str(T) ' secs; computed'],[rmin rmax],6);
    subplot(2,2,4);
    InterpPlot(lon,lat,imag(Z2(:,2)),'Im(Zxy)',[num2str(T) ' secs; computed'],[imin imax],6);
    strper = sprintf('%03d',k);
    print('-r300','-djpeg',[filename '_Zxy_' strper '.jpg']);
    % Zyx
    clf
    rmin = min(real(Z1(:,3))); rmax = max(real(Z1(:,3)));
    imin = min(imag(Z1(:,3))); imax = max(imag(Z1(:,3)));
    subplot(2,2,1);
    InterpPlot(lon,lat,real(Z1(:,3)),'Re(Zyx)',[num2str(T) ' secs; measured'],[rmin rmax],6);
    subplot(2,2,3);
    InterpPlot(lon,lat,imag(Z1(:,3)),'Im(Zyx)',[num2str(T) ' secs; measured'],[imin imax],6);
    subplot(2,2,2);
    InterpPlot(lon,lat,real(Z2(:,3)),'Re(Zyx)',[num2str(T) ' secs; computed'],[rmin rmax],6);
    subplot(2,2,4);
    InterpPlot(lon,lat,imag(Z2(:,3)),'Im(Zyx)',[num2str(T) ' secs; computed'],[imin imax],6);
    strper = sprintf('%03d',k);
    print('-r300','-djpeg',[filename '_Zyx_' strper '.jpg']);
    % Zyy
    clf
    rmin = min(real(Z1(:,4))); rmax = max(real(Z1(:,4)));
    imin = min(imag(Z1(:,4))); imax = max(imag(Z1(:,4)));
    subplot(2,2,1);
    InterpPlot(lon,lat,real(Z1(:,4)),'Re(Zyy)',[num2str(T) ' secs; measured'],[rmin rmax],6);
    subplot(2,2,3);
    InterpPlot(lon,lat,imag(Z1(:,4)),'Im(Zyy)',[num2str(T) ' secs; measured'],[imin imax],6);
    subplot(2,2,2);
    InterpPlot(lon,lat,real(Z2(:,4)),'Re(Zyy)',[num2str(T) ' secs; computed'],[rmin rmax],6);
    subplot(2,2,4);
    InterpPlot(lon,lat,imag(Z2(:,4)),'Im(Zyy)',[num2str(T) ' secs; computed'],[imin imax],6);
    strper = sprintf('%03d',k);
    print('-r300','-djpeg',[filename '_Zyy_' strper '.jpg']);
    % Tx
    clf
    rmin = min(real(Z1(:,5))); rmax = max(real(Z1(:,5)));
    imin = min(imag(Z1(:,5))); imax = max(imag(Z1(:,5)));
    subplot(2,2,1);
    InterpPlot(lon,lat,real(Z1(:,5)),'Re(Tx)',[num2str(T) ' secs; measured'],[rmin rmax],6);
    subplot(2,2,3);
    InterpPlot(lon,lat,imag(Z1(:,5)),'Im(Tx)',[num2str(T) ' secs; measured'],[imin imax],6);
    subplot(2,2,2);
    InterpPlot(lon,lat,real(Z2(:,5)),'Re(Tx)',[num2str(T) ' secs; computed'],[rmin rmax],6);
    subplot(2,2,4);
    InterpPlot(lon,lat,imag(Z2(:,5)),'Im(Tx)',[num2str(T) ' secs; computed'],[imin imax],6);
    strper = sprintf('%03d',k);
    print('-r300','-djpeg',[filename '_Tx_' strper '.jpg']);
    % Ty
    clf
    rmin = min(real(Z1(:,6))); rmax = max(real(Z1(:,6)));
    imin = min(imag(Z1(:,6))); imax = max(imag(Z1(:,6)));
    subplot(2,2,1);
    InterpPlot(lon,lat,real(Z1(:,6)),'Re(Ty)',[num2str(T) ' secs; measured'],[rmin rmax],6);
    subplot(2,2,3);
    InterpPlot(lon,lat,imag(Z1(:,6)),'Im(Ty)',[num2str(T) ' secs; measured'],[imin imax],6);
    subplot(2,2,2);
    InterpPlot(lon,lat,real(Z2(:,6)),'Re(Ty)',[num2str(T) ' secs; computed'],[rmin rmax],6);
    subplot(2,2,4);
    InterpPlot(lon,lat,imag(Z2(:,6)),'Im(Ty)',[num2str(T) ' secs; computed'],[imin imax],6);
    strper = sprintf('%03d',k);
    print('-r300','-djpeg',[filename '_Ty_' strper '.jpg']);
end
