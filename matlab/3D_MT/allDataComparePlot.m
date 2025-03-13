fdir = '.';
filt = [fdir '/*.dat'];
[filename, pathname] = uigetfile(filt, 'Impedance File: Real Data');
fdata = [pathname filename];
%data = readZ_3D(fdata,'[mV/km]/[nT]');
[filename, pathname] = uigetfile(filt, 'Impedance File: Responses');
fresp = [pathname filename];
%resp = readZ_3D(fresp,'[mV/km]/[nT]');
filename = filename(1:end-4);
[pathstr, fplot] = fileparts(fresp);
% fid = fopen('../final_sites.txt','r');
% tmp = textscan(fid,'%s %f %f');
% fclose(fid);
% sitenames = char(tmp{1});
% [rms,info] = DataFit(fresp,fdata,5,fplot,sitenames);
%[rms,info] = DataFit(fresp,fdata,5,fplot);
[rms,info] = DataFit(fresp,fdata,4);

info.data(isnan(info.data)) = NaN+1i*NaN;
info.resp(isnan(info.resp)) = NaN+1i*NaN;

% save site list in tab separated *.txt file
fid = fopen([fplot '.txt'],'w');
fprintf(fid,'site\tlatitude\tlongitude\trms\tsurvey\n');
for i = 1:length(info.code)
    if findstr(info.code(i,:),'SR')
        survey = 1;
    else
        survey = 0;
    end
    ncomp = size(info.res,3);
    nper = size(info.res,2);
    count = ncomp*nper - sum(sum(isnan(info.res(i,:,:))));
    if count > 0
        rms = sqrt(nansum(nansum(info.res(i,:,:)))/count);
        fprintf(fid,'%s\t%f\t%f\t%f\t%f\n',info.code(i,:),info.lat(i),info.lon(i),rms,survey);
    end
end
fclose(fid);

%%
total_rms = sqrt( nansum(nansum(nansum(info.res))) / sum(sum(sum(~isnan(info.res)))) );

ES(1,length(info.code)) = 0;
for i = 1:length(info.code)
    if isempty(findstr(info.code(i,:),'SR'))
        ES(i) = 1;
    end
end
ES_rms = sqrt( nansum(nansum(nansum(info.res(ES==1,:,:)))) / sum(sum(sum(~isnan(info.res(ES==1,:,:))))) );
SR_rms = sqrt( nansum(nansum(nansum(info.res(ES==0,:,:)))) / sum(sum(sum(~isnan(info.res(ES==0,:,:))))) );

ZXX_rms = sqrt( nansum(nansum(info.res(:,:,1))) / sum(sum(~isnan(info.res(:,:,1)))) );
ZXY_rms = sqrt( nansum(nansum(info.res(:,:,2))) / sum(sum(~isnan(info.res(:,:,2)))) );
ZYX_rms = sqrt( nansum(nansum(info.res(:,:,3))) / sum(sum(~isnan(info.res(:,:,3)))) );
ZYY_rms = sqrt( nansum(nansum(info.res(:,:,4))) / sum(sum(~isnan(info.res(:,:,4)))) );
Z_rms = sqrt( nansum(nansum(nansum(info.res(:,:,1:4)))) / sum(sum(sum(~isnan(info.res(:,:,1:4))))) );
ESZ_rms = sqrt( nansum(nansum(nansum(info.res(ES==1,:,1:4)))) / sum(sum(sum(~isnan(info.res(ES==1,:,1:4))))) );
TX_rms = sqrt( nansum(nansum(info.res(:,:,5))) / sum(sum(~isnan(info.res(:,:,5)))) );
TY_rms = sqrt( nansum(nansum(info.res(:,:,6))) / sum(sum(~isnan(info.res(:,:,6)))) );
T_rms = sqrt( nansum(nansum(nansum(info.res(:,:,5:6)))) / sum(sum(sum(~isnan(info.res(:,:,5:6))))) );

disp( ' Subset  |  RMS'); 
disp( '------------------'); 
disp([' ALL     |  ' num2str(total_rms)]); 
disp(['  ES     |  ' num2str(ES_rms)]); 
disp(['  SR     |  ' num2str(SR_rms)]); 
disp([' ZXX     |  ' num2str(ZXX_rms)]); 
disp([' ZXY     |  ' num2str(ZXY_rms)]); 
disp([' ZYX     |  ' num2str(ZYX_rms)]); 
disp([' ZYY     |  ' num2str(ZYY_rms)]); 
disp([' ALL Z   |  ' num2str(Z_rms)]); 
disp(['  ES Z   |  ' num2str(ESZ_rms)]); 
disp(['  TX     |  ' num2str(TX_rms)]); 
disp(['  TY     |  ' num2str(TY_rms)]); 
disp([' ALL T   |  ' num2str(T_rms)]); 

%%
% resp = readZ_3D(fresp,'[mV/km]/[nT]');
% data = readZ_3D(fdata,'[mV/km]/[nT]');
% [rms,info,site] = DataFitZ(resp,data,'ZXY',5,fplot,sitenames);

%%
OBS = 'WYYS3';%'WYH18','WYYS3','WYJ18','WYI18','MTG17';
[str,i] = intersect(info.code,OBS,'rows');
if isempty(i)
    error(['No observatory in data file matching ' char(OBS)]);
end
%%
for i=1:length(info.lon)
OBS = char(info.code(i,:));
lon = num2str(info.lon(i));
lat = num2str(info.lat(i));
per = info.per;
range = 100; %20;
xlims = [0 1000]; %[5 10^5];
ncomp = size(info.data,3)/2;

f=figure('Position',[200,200,1400,900],...
        'PaperPosition',[1,1,18,10],...
        'PaperOrientation','Portrait'); clf; 
% Diagonal TFs
subplot(2,3,1)
data = info.data(i,:,1).*sqrt(per);
resp = info.resp(i,:,1).*sqrt(per);
std2 = 2*info.err(i,:,1).*sqrt(per);
TFrms = sqrt(nansum(info.res(i,:,1))/sum(~isnan(info.err(i,:,1))));
TFinfo = [char(OBS) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
p(1)=semilogx(per,real(resp),'b-','linewidth',1.5); hold on
errorbar(per,real(data),std2,'bx','linewidth',1.5); hold on
p(2)=semilogx(per,imag(resp),'--','color',[0.8 0.5 0.],'linewidth',1.5); hold on
errorbar(per,imag(data),std2,'o','color',[0.8 0.5 0.],'linewidth',1.5); hold off
ymean = real(nanmedian(data));
set(gca,'xlim',xlims,'ylim',[ymean-range/2 ymean+range/2]);
set(gca,'xminortick','on','yminortick','on');
legend(p,'Re(ZXX)','Im(ZXX)');
ylabel('ZXX * sqrt(T)','fontweight','demi');
xlabel('Period (secs)','fontweight','demi');
title(TFinfo,'fontweight','demi');
subplot(2,3,5)
data = info.data(i,:,4).*sqrt(per);
resp = info.resp(i,:,4).*sqrt(per);
std2 = 2*info.err(i,:,4).*sqrt(per);
TFrms = sqrt(nansum(info.res(i,:,4))/sum(~isnan(info.err(i,:,4))));
TFinfo = [char(OBS) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
p(1)=semilogx(per,real(resp),'b-','linewidth',1.5); hold on
errorbar(per,real(data),std2,'bx','linewidth',1.5); hold on
p(2)=semilogx(per,imag(resp),'--','color',[0.8 0.5 0.],'linewidth',1.5); hold on
errorbar(per,imag(data),std2,'o','color',[0.8 0.5 0.],'linewidth',1.5); hold off
ymean = real(nanmedian(data));
set(gca,'xlim',xlims,'ylim',[ymean-range/2 ymean+range/2]);
set(gca,'xminortick','on','yminortick','on');
legend(p,'Re(ZYY)','Im(ZYY)');
ylabel('ZYY * sqrt(T)','fontweight','demi');
xlabel('Period (secs)','fontweight','demi');
title(TFinfo,'fontweight','demi');
% Off-diagonal TFs
subplot(2,3,2)
data = info.data(i,:,2).*sqrt(per);
resp = info.resp(i,:,2).*sqrt(per);
std2 = 2*info.err(i,:,2).*sqrt(per);
TFrms = sqrt(nansum(info.res(i,:,2))/sum(~isnan(info.err(i,:,2))));
TFinfo = [char(OBS) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
p(1)=semilogx(per,real(resp),'r-','linewidth',1.5); hold on
errorbar(per,real(data),std2,'rx','linewidth',1.5); hold on
p(2)=semilogx(per,imag(resp),'--','color',[0. 0.5 0.8],'linewidth',1.5); hold on
errorbar(per,imag(data),std2,'o','color',[0. 0.5 0.8],'linewidth',1.5); hold off
ymean = real(nanmedian(data));
set(gca,'xlim',xlims,'ylim',[ymean-range/2 ymean+range/2]);
set(gca,'xminortick','on','yminortick','on');
legend(p,'Re(ZXY)','Im(ZXY)');
ylabel('ZXY * sqrt(T)','fontweight','demi');
xlabel('Period (secs)','fontweight','demi');
title(TFinfo,'fontweight','demi');
subplot(2,3,4)
data = info.data(i,:,3).*sqrt(per);
resp = info.resp(i,:,3).*sqrt(per);
std2 = 2*info.err(i,:,3).*sqrt(per);
TFrms = sqrt(nansum(info.res(i,:,3))/sum(~isnan(info.err(i,:,3))));
TFinfo = [char(OBS) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
p(1)=semilogx(per,real(resp),'r-','linewidth',1.5); hold on
errorbar(per,real(data),std2,'rx','linewidth',1.5); hold on
p(2)=semilogx(per,imag(resp),'--','color',[0. 0.5 0.8],'linewidth',1.5); hold on
errorbar(per,imag(data),std2,'o','color',[0. 0.5 0.8],'linewidth',1.5); hold off
ymean = real(nanmedian(data));
set(gca,'xlim',xlims,'ylim',[ymean-range/2 ymean+range/2]);
set(gca,'xminortick','on','yminortick','on');
legend(p,'Re(ZYX)','Im(ZYX)');
ylabel('ZYX * sqrt(T)','fontweight','demi');
xlabel('Period (secs)','fontweight','demi');
title(TFinfo,'fontweight','demi');
% Vertical TFs
% subplot(2,3,3)
% data = info.data(i,:,5);
% resp = info.resp(i,:,5);
% std2 = 2*info.err(i,:,5);
% TFrms = sqrt(nansum(info.res(i,:,5))/sum(~isnan(info.err(i,:,5))));
% TFinfo = [char(OBS) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
% p(1)=semilogx(per,real(resp),'m-','linewidth',1.5); hold on
% errorbar(per,real(data),std2,'mx','linewidth',1.5); hold on
% p(2)=semilogx(per,imag(resp),'g--','linewidth',1.5); hold on
% errorbar(per,imag(data),std2,'go','linewidth',1.5); hold off
% ymean = 0; %real(nanmedian(data));
% set(gca,'xlim',xlims,'ylim',[ymean-0.5 ymean+0.5]);
% set(gca,'xminortick','on','yminortick','on');
% legend(p,'Re(TX)','Im(TX)');
% ylabel('TX','fontweight','demi');
% xlabel('Period (secs)','fontweight','demi');
% title(TFinfo,'fontweight','demi');
% subplot(2,3,6)
% data = info.data(i,:,6);
% resp = info.resp(i,:,6);
% std2 = 2*info.err(i,:,6);
% TFrms = sqrt(nansum(info.res(i,:,6))/sum(~isnan(info.err(i,:,6))));
% TFinfo = [char(OBS) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
% p(1)=semilogx(per,real(resp),'m-','linewidth',1.5); hold on
% errorbar(per,real(data),std2,'mx','linewidth',1.5); hold on
% p(2)=semilogx(per,imag(resp),'g--','linewidth',1.5); hold on
% errorbar(per,imag(data),std2,'go','linewidth',1.5); hold off
% ymean = 0; %real(nanmedian(data));
% set(gca,'xlim',xlims,'ylim',[ymean-0.5 ymean+0.5]);
% set(gca,'xminortick','on','yminortick','on');
% legend(p,'Re(TY)','Im(TY)');
% ylabel('TY','fontweight','demi');
% xlabel('Period (secs)','fontweight','demi');
% title(TFinfo,'fontweight','demi');
print(f,'-dpng',[pathname '/' filename '_fit_' OBS '_TFs']);
close(f);
end


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

% rescale apparent resistivity by period
for l = 1:length(info.per)
  apres.yx(:,l) = apres.yx(:,l)*info.per(l)/5. ;
  apres.xy(:,l) = apres.xy(:,l)*info.per(l)/5. ;
  apres.yx_re(:,l) = apres.yx_re(:,l)*info.per(l)/5. ;
  apres.xy_re(:,l) = apres.xy_re(:,l)*info.per(l)/5. ;
  apres.yx_se(:,l) = sqrt(apres.yx_se(:,l).*apres.yx(:,l)*info.per(l)*4/5.);
  apres.xy_se(:,l) = sqrt(apres.xy_se(:,l).*apres.xy(:,l)*info.per(l)*4/5.);
end

% also plot the vertical transfer functions
tx = info.data(:,:,5);
tx_se = info.err(:,:,5);
tx_re = info.resp(:,:,5);
ty = info.data(:,:,6);
ty_se = info.err(:,:,6);
ty_re = info.resp(:,:,6);

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

f=figure('Position',[200,200,450,900],...
        'PaperPosition',[1,1,5,10],...
        'PaperOrientation','Portrait'); clf; 
subplot(3,1,1)
p(1)=loglog(per,apres.xy_re(i,:),'-','color',[0. 0.5 0.8],'linewidth',1.5); hold on
errorbar(per,apres.xy(i,:),2*apres.xy_se(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
p(2)=loglog(per,apres.yx_re(i,:),'r-','linewidth',1.5); hold on
errorbar(per,apres.yx(i,:),2*apres.yx_se(i,:),'rx','linewidth',1.5); hold off
ymean = (nanmedian(apres.xy(i,:))+nanmedian(apres.yx(i,:)))/2;
%set(gca,'xlim',[5 10^5],'ylim',[0.1 20000]);
set(gca,'xlim',[5 10^5],'ylim',[ymean/20 ymean*20]);
set(gca,'xminortick','on','ytick',[1 10 100 1000 10000]);
legend(p,'XY','YX');
ylabel('\rho_a (\Omega m)','fontweight','demi');
RMS = sqrt(nansum(info.res(i,:))/sum(~isnan(info.err(i,:))));
title([char(OBS) ' [ LON = ' lon '; LAT = ' lat '; RMS = ' num2str(RMS)  ' ]'],'fontweight','demi');
subplot(3,1,2)
p(1)=semilogx(per,phase.xy_re(i,:),'-','color',[0. 0.5 0.8],'linewidth',1.5); hold on
errorbar(per,phase.xy(i,:),2*phase.xy_se(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
p(2)=semilogx(per,phase.yx_re(i,:),'r-','linewidth',1.5); hold on
errorbar(per,phase.yx(i,:),2*phase.yx_se(i,:),'rx','linewidth',1.5); hold off
ymean = (nanmean(phase.xy(i,:))+nanmean(phase.yx(i,:)))/2;
%set(gca,'xlim',[5 10^5],'ylim',[-90 90]);
set(gca,'xlim',[5 10^5],'ylim',[ymean-45 ymean+45]);
set(gca,'xminortick','on','yminortick','on');
ylabel('\phi (Degrees)','fontweight','demi');
xlabel('Period (secs)','fontweight','demi');
subplot(3,1,3)
p(1)=semilogx(per,tx_re(i,:),'g-','linewidth',1.5); hold on
errorbar(per,tx(i,:),2*tx_se(i,:),'go','linewidth',1.5); hold on
p(2)=semilogx(per,ty_re(i,:),'m-','linewidth',1.5); hold on
errorbar(per,ty(i,:),2*ty_se(i,:),'mx','linewidth',1.5); hold off
ymean = (nanmean(ty(i,:))+nanmean(tx(i,:)))/2;
%set(gca,'xlim',[5 10^5],'ylim',[-1 1]);
set(gca,'xlim',[5 10^5],'ylim',[ymean-0.5 ymean+0.5]);
set(gca,'xminortick','on','yminortick','on');
legend(p,'TX','TY');
ylabel('Vertical TFs','fontweight','demi');
xlabel('Period (secs)','fontweight','demi');
print(f,'-dpng',[pathname '/apres/' filename '_fit_' OBS '_apres']);
close(f);
end
% %% 
% rho2 = [apres.xy(i,:); apres.yx(i,:)].';
% rho2_se = 2*[apres.xy_se(i,:); apres.yx_se(i,:)].';
% ph2 = [phase.xy(i,:); phase.yx(i,:)].';
% ph2_se = 2*[phase.xy_se(i,:); phase.yx_se(i,:)].';
% pltall = ones(length(per),1);
% [rho2_axes,ph2_axes] = pltrhom(length(per),pltall,per,rho2,rho2_se,ph2,ph2_se,lims,c_title_Modes,hfig);
