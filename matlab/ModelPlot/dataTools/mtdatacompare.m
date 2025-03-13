classdef mtdatacompare < mtdata
    %   (c) Anna Kelbert, April 2014
    
    properties (SetAccess = protected) % already defined in mtdata
        % simple structure of all data regrouped into observatory bins;
        % most useful for plotting
        %v
        % for impedances, sometimes it is convenient to also compute 
        % and store the apparent resistivities and phases
        %apres
        %phase
    end
    
    methods
        
        function [obj] = mtdatacompare(obj1,obj2)
            %   class constructor
            if obj1.predicted 
                warning('The first object in constructor must contain real data')
            end
            if ~obj2.predicted
                warning('The second object in constructor must contain predicted data')
            end
            
            if ~ strcmp(obj2.units,obj1.units)
                obj2 = obj2.changeUnits(obj1.units);
            end
            
            if obj2.signConvention ~= obj1.signConvention
                obj2 = obj2.changeSignConvention;
            end         
            
            if ~ strcmp(obj2.type,obj1.type)
                warning('MT data types are incompatible');
            end
           
            type = obj1.type;
            obj = obj@mtdata(type);
            obj = obj.regroup(obj1,obj2);    
            obj.predicted = obj1.predicted && obj2.predicted;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,info,apres,phase] = regroup(obj,obj1,obj2,sitenames)
                                                            
            % Make a list of all sites in the real data file
            nper = obj1.nPeriods;
            if obj1.Cmplx
                ncomp = obj1.nComp/2;
            else
                ncomp = obj1.nComp;
            end
            allsites = obj1.d{1}.siteChar;
            allsitesloc = obj1.d{1}.siteLoc;
            lat = obj1.d{1}.lat;
            lon = obj1.d{1}.lon;
            for j = 2:nper
                for i = 1:length(obj1.d{j}.siteChar)
                    if ~contains(allsites,obj1.d{j}.siteChar{i})
                        allsites = [allsites; obj1.d{j}.siteChar{i}];
                        allsitesloc = [allsitesloc; obj1.d{j}.siteLoc(i,:)];
                        lat = [lat; obj1.d{j}.lat(i)];
                        lon = [lon; obj1.d{j}.lon(i)];
                    end
                end
            end
            nsites = length(allsites);
            x = allsitesloc(:,1);
            y = allsitesloc(:,2);
            %[lat,lon] = xy2latlon(x,y,obj1.d{1}.origin(1),obj1.d{1}.origin(2),'m');
            
            % Regroup into observatory bins
            info.data = nan(length(allsites),nper,ncomp)+1i*nan(length(allsites),nper,ncomp);
            info.resp = nan(length(allsites),nper,ncomp)+1i*nan(length(allsites),nper,ncomp);
            info.err = nan(length(allsites),nper,ncomp);
            for j = 1:nper
                j2 = findPeriod(obj2,obj1.periods(j));
                for i = 1:length(allsites)
                    if nargin > 3
                        in = 0;
                        k=find(contains(sitenames,allsites{i}));
                        if ~isempty(k)
                            in = 1;
                        end
                    else
                        in = 1;
                    end
                    k=find(contains(obj1.d{j}.siteChar,strtrim(allsites{i})));
                    if ~isempty(k) && in
                        info.code(i,:) = allsites{i};
                        info.loc(i,:) = allsitesloc(i,:);
                        info.lon(i) = lon(i);
                        info.lat(i) = lat(i);
                        info.data(i,j,:) = obj1.d{j}.TF(k,:);
                        info.err(i,j,:) = obj1.d{j}.TFerr(k,:);
                        info.per(j) = obj1.d{j}.T;
                    end
                    if j2>0
                        k=find(contains(obj2.d{j2}.siteChar,strtrim(allsites{i})));
                        if ~isempty(k) && in
                            info.resp(i,j,:) = obj2.d{j2}.TF(k,:);
                        end
                    end
                end
            end
            info.res = 0.5*(real(info.data - info.resp)./info.err).^2 + ...
                0.5*(imag(info.data - info.resp)./info.err).^2;
            
            info.comp = obj.compChar;
            
            % store the result in obj
            obj.v = info;
            
            % also compute apparent resistivities and phases if reasonable
            % info = regroup(obj);
            rad_deg = 57.2958;
            if strcmp(obj.type,'Off_Diagonal_Rho_Phase')
                apres.xy = squeeze(info.data(:,:,1));
                apres.xy_re = squeeze(info.resp(:,:,1));
                apres.xy_se = squeeze(info.err(:,:,1));
                phase.xy = squeeze(info.data(:,:,2));
                phase.xy_re = squeeze(info.resp(:,:,2));
                phase.xy_se = squeeze(info.err(:,:,2));
                apres.yx = squeeze(info.data(:,:,3));
                apres.yx_re = squeeze(info.resp(:,:,3));
                apres.yx_se = squeeze(info.err(:,:,3));
                phase.yx = squeeze(info.data(:,:,4));
                phase.yx_re = squeeze(info.resp(:,:,4));
                phase.yx_se = squeeze(info.err(:,:,4));
                return
            elseif strcmp(obj.type,'Full_Impedance')
                ixy = 2;
                iyx = 3;
            elseif strcmp(obj.type,'Off_Diagonal_Impedance')
                ixy = 1;
                iyx = 2;
            else
                warning('Unable to compute apparent resistivities and phases for this data type');
            end
            apres.xy = abs(info.data(:,:,ixy)).^2;
            apres.xy_re = abs(info.resp(:,:,ixy)).^2;
            apres.xy_se = info.err(:,:,ixy).^2; % variance of Z as in EDI file
            apres.yx = abs(info.data(:,:,iyx)).^2;
            apres.yx_re = abs(info.resp(:,:,iyx)).^2;
            apres.yx_se = info.err(:,:,iyx).^2; % variance of Z as in EDI file
            % for now, using atan2 which takes values (-pi,pi] as in ModEM
            phase.xy = rad_deg*atan2(imag(info.data(:,:,ixy)),real(info.data(:,:,ixy)));
            phase.xy_re = rad_deg*atan2(imag(info.resp(:,:,ixy)),real(info.resp(:,:,ixy)));
            phase.xy_se = rad_deg*sqrt(apres.xy_se./apres.xy);
            % for now, using atan2 which takes values (-pi,pi] as in ModEM
            phase.yx = rad_deg*atan2(imag(info.data(:,:,iyx)),real(info.data(:,:,iyx))) + 180.;
            phase.yx_re = rad_deg*atan2(imag(info.resp(:,:,iyx)),real(info.resp(:,:,iyx))) + 180.;
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
            
            % store the result in obj
            obj.apres = apres;
            obj.phase = phase;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = aprespltall(obj,fileName)
            % create and apparent resistivity and phase object from an
            % impedance object; can write the new object to file, or plot
            
            if isempty(obj.v)
                error('Call obj = mtdatacompare(data,pred) first to initialize');
            end
            
            for i=1:length(obj.v.lon)
                siteName = char(obj.v.code(i,:));
                if nargin > 1
                    f = apresplt(obj,siteName,fileName);
                    close(f);
                else
                    f = apresplt(obj,siteName);
                    input('Press enter to continue ...');
                    close(f);
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = apresplt(obj,siteName,fileName)
            % plot apparent resistivities and phases
            
            if isempty(obj.v)
                error('Call obj = mtdatacompare(data,pred) first to initialize');
            end
            
            [~,i] = intersect(obj.v.code,siteName,'rows');
            if isempty(i)
                error(['No site in MT data object matching ' char(siteName)]);
            end
            
            lon = num2str(obj.v.lon(i));
            lat = num2str(obj.v.lat(i));
            per = obj.v.per;

            f=figure('Position',[200 200 500 800],'PaperPosition',[2 2 5 8]); clf;
            subplot(3,1,[1 2]);
            p(1)=loglog(per,obj.apres.xy_re(i,:),'-','color',[0. 0.5 0.8],'linewidth',1.5); hold on
            errorbar(per,obj.apres.xy(i,:),2*obj.apres.xy_se(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
            p(2)=loglog(per,obj.apres.yx_re(i,:),'r-','linewidth',1.5); hold on
            errorbar(per,obj.apres.yx(i,:),2*obj.apres.yx_se(i,:),'rx','linewidth',1.5); hold off
            ymean = (nanmedian(obj.apres.xy(i,:))+nanmedian(obj.apres.yx(i,:)))/2;
            %set(gca,'xlim',[5 10^5],'ylim',[0.1 20000]);
            set(gca,'xlim',[5 10^5],'ylim',[ymean/20 ymean*20]);
            set(gca,'xminortick','on','ytick',[1 10 100 1000 10000]);
            legend(p,'XY','YX');
            ylabel('\rho_a (\Omega m)','fontsize',14,'fontweight','demi');
            RMS = sqrt(nansum(obj.v.res(i,:))/sum(~isnan(obj.v.err(i,:))));
            title([char(siteName) ' [ LON = ' lon '; LAT = ' lat '; RMS = ' num2str(RMS)  ' ]'],'fontweight','demi');
            subplot(3,1,3);
            p(1)=semilogx(per,obj.phase.xy_re(i,:),'-','color',[0. 0.5 0.8],'linewidth',1.5); hold on
            errorbar(per,obj.phase.xy(i,:),2*obj.phase.xy_se(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
            p(2)=semilogx(per,obj.phase.yx_re(i,:),'r-','linewidth',1.5); hold on
            errorbar(per,obj.phase.yx(i,:),2*obj.phase.yx_se(i,:),'rx','linewidth',1.5); hold off
            ymean = (nanmean(obj.phase.xy(i,:))+nanmean(obj.phase.yx(i,:)))/2;
            %set(gca,'xlim',[5 10^5],'ylim',[-90 90]);
            set(gca,'xlim',[5 10^5],'ylim',[ymean-45 ymean+45]);
            set(gca,'xminortick','on','yminortick','on');
            ylabel('\phi (Degrees)','fontsize',14,'fontweight','demi');
            xlabel('Period (secs)','fontsize',14,'fontweight','demi');

            if nargin > 2
                print(f,'-dpng',[fileName '_' siteName '_apres']);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = impplt(obj,siteName,fileName)
            % plot the transfer functions directly
            
            if isempty(obj.v)
                error('Call obj = mtdatacompare(data,pred) first to initialize');
            end
            
            i = find(contains(obj.v.code,siteName));
            if isempty(i)
                error(['No site in MT data object matching ' char(siteName)]);
            end
            
            lon = num2str(obj.v.lon(i));
            lat = num2str(obj.v.lat(i));
            per = obj.v.per;
            range = 20; % 100; 
            xlims = [5 10^5]; %[0 1000];
            %ncomp = size(obj.v.data,3)/2;
            
            if strcmp(obj.type,'Full_Impedance')
                f=figure('Position',[200,200,1400,900],...
                    'PaperPosition',[1,1,18,10],...
                    'PaperOrientation','Portrait'); clf;
                % Diagonal TFs
                subplot(2,2,1)
                data = obj.v.data(i,:,1).*sqrt(per);
                resp = obj.v.resp(i,:,1).*sqrt(per);
                std2 = 2*obj.v.err(i,:,1).*sqrt(per);
                TFrms = sqrt(nansum(obj.v.res(i,:,1))/sum(~isnan(obj.v.err(i,:,1))));
                TFinfo = [char(siteName) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
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
                subplot(2,2,4)
                data = obj.v.data(i,:,4).*sqrt(per);
                resp = obj.v.resp(i,:,4).*sqrt(per);
                std2 = 2*obj.v.err(i,:,4).*sqrt(per);
                TFrms = sqrt(nansum(obj.v.res(i,:,4))/sum(~isnan(obj.v.err(i,:,4))));
                TFinfo = [char(siteName) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
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
                subplot(2,2,2)
                data = obj.v.data(i,:,2).*sqrt(per);
                resp = obj.v.resp(i,:,2).*sqrt(per);
                std2 = 2*obj.v.err(i,:,2).*sqrt(per);
                TFrms = sqrt(nansum(obj.v.res(i,:,2))/sum(~isnan(obj.v.err(i,:,2))));
                TFinfo = [char(siteName) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
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
                subplot(2,2,3)
                data = obj.v.data(i,:,3).*sqrt(per);
                resp = obj.v.resp(i,:,3).*sqrt(per);
                std2 = 2*obj.v.err(i,:,3).*sqrt(per);
                TFrms = sqrt(nansum(obj.v.res(i,:,3))/sum(~isnan(obj.v.err(i,:,3))));
                TFinfo = [char(siteName) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
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
            elseif strcmp(obj.type,'Full_Vertical_Components')
                % Vertical TFs
                f=figure('Position',[200,200,400,900],...
                    'PaperPosition',[1,1,5,10],...
                    'PaperOrientation','Portrait'); clf;
                subplot(2,1,1)
                data = obj.v.data(i,:,1);
                resp = obj.v.resp(i,:,1);
                std2 = 2*obj.v.err(i,:,1);
                TFrms = sqrt(nansum(obj.v.res(i,:,1))/sum(~isnan(obj.v.err(i,:,1))));
                TFinfo = [char(siteName) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
                p(1)=semilogx(per,real(resp),'m-','linewidth',1.5); hold on
                errorbar(per,real(data),std2,'mx','linewidth',1.5); hold on
                p(2)=semilogx(per,imag(resp),'g--','linewidth',1.5); hold on
                errorbar(per,imag(data),std2,'go','linewidth',1.5); hold off
                ymean = 0; %real(nanmedian(data));
                set(gca,'xlim',xlims,'ylim',[ymean-0.5 ymean+0.5]);
                set(gca,'xminortick','on','yminortick','on');
                legend(p,'Re(TX)','Im(TX)');
                ylabel('TX','fontweight','demi');
                xlabel('Period (secs)','fontweight','demi');
                title(TFinfo,'fontweight','demi');
                subplot(2,1,2)
                data = obj.v.data(i,:,2);
                resp = obj.v.resp(i,:,2);
                std2 = 2*obj.v.err(i,:,2);
                TFrms = sqrt(nansum(obj.v.res(i,:,2))/sum(~isnan(obj.v.err(i,:,2))));
                TFinfo = [char(siteName) ' (' lon '; ' lat '); RMS = ' num2str(TFrms)];
                p(1)=semilogx(per,real(resp),'m-','linewidth',1.5); hold on
                errorbar(per,real(data),std2,'mx','linewidth',1.5); hold on
                p(2)=semilogx(per,imag(resp),'g--','linewidth',1.5); hold on
                errorbar(per,imag(data),std2,'go','linewidth',1.5); hold off
                ymean = 0; %real(nanmedian(data));
                set(gca,'xlim',xlims,'ylim',[ymean-0.5 ymean+0.5]);
                set(gca,'xminortick','on','yminortick','on');
                legend(p,'Re(TY)','Im(TY)');
                ylabel('TY','fontweight','demi');
                xlabel('Period (secs)','fontweight','demi');
                title(TFinfo,'fontweight','demi');
            end
            if nargin > 2
                print(f,'-dpng',[fileName '_fit_' siteName '_TFs']);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [rms,info,freq,site] = misfit(obj,maxrms,fplot,prange)
            
            % [rms,info,freq,site] = misfit(obj,maxrms,fplot,prange)
            %
            % Outputs the RMS misfit for two data objects. Plots a summary.
                        
            PLOT = 1;
            if nargin < 2
                PLOT = 0;
            end
            
            %obj = mtdatacompare(obj1,obj2);
            info = obj.v;
            nper = length(info.per);
            nsites = length(info.code);
            ncomp = length(info.comp);
            
            if nargin > 3
                minper = prange(1);
                maxper = prange(2);
            else
                minper = 1e-4;
                maxper = 1e6;
            end

            % Compute total RMS
            count = 0;
            misfit = 0;
            for j = 1:nper
                if info.per(j) < maxper && info.per(j) > minper
                    for i = 1:nsites
                        newcount = ncomp - sum(isnan(info.err(i,j,:)));
                        count = count + newcount;
                        if newcount > 0
                            misfit = misfit + nansum(info.res(i,j,:));
                        end
                    end
                end
            end
            rms = sqrt(misfit/count);
            
            if ~PLOT
                return
            end
            
            % Average over frequency and component bins
            for j = 1:nper
                freq.f(j) = 1/(info.per(j));
                freq.p(j) = info.per(j);
                freq.d(j) = info.per(j)/(24*3600);
                for k = 1:ncomp
                    count = nsites - sum(isnan(info.err(:,j,k)));
                    freq.rms(j,k) = sqrt(nansum(info.res(:,j,k))/count);
                end
            end
            
            % Average over site bins
            j = 1;
            for i = 1:nsites
                count = ncomp*nper - sum(sum(isnan(info.res(i,:,:))));
                if count > 0
                    site.codes(j,:) = info.code(i,:);
                    site.lon(j) = info.lon(i);
                    site.lat(j) = info.lat(i);
                    site.rms(j) = sqrt(nansum(nansum(info.res(i,:,:)))/count);
                    j = j+1;
                end
            end
            
            % Plot
            figure('Position',[300,300,600,600],...
                'PaperPosition',[1,1,6,6],...
                'PaperOrientation','Portrait');
            subplot(7,1,1:5);
            htitle = ['RMS ' num2str(rms) ' (averaged over ' num2str(nper) ' periods)'];
            InterpPlot(site.lon,site.lat,site.rms,'RMS',htitle);
            colormap(invhot);
            if nargin > 2
                caxis([1 maxrms]);
            end            
            subplot(7,1,6:7);
            h = zeros(1,6);
            for i = 1:ncomp
                color = [(i-1)/ncomp (i-1)/ncomp (i-1)/ncomp];
                h(i)=semilogx(freq.p,freq.rms(:,i),'x-','color',color); hold on;
                set(h(i),'LineWidth',2);
            end
            xlabel('Log(10) T (secs)');
            ylabel('RMS');
            legend(h,info.comp,'Location','BestOutside');
            if nargin > 2
                print('-djpeg95','-r300',[fplot '_rms.jpg']);
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = plot(obj,comp,iper,fplot)
            
            % [f] = plot(obj,comp,iper,fplot)
            %
            % Plots data1, data2, and difference as interpolated map.
            % Requires comp and iper.
            % If fplot is specified, saves a jpeg.
                                   
            %obj = mtdatacompare(obj1,obj2);
            %info = obj.v;
            [~, info] = misfit(obj);
            info.diff = info.data - info.resp;
            
            if nargin < 2
                disp('Usage: [f] = obj.plot(comp,iper,fplot)');
                return
            end
            
            nper = length(info.per);

            if nargin >= 2
                [~,icomp]=intersect(info.comp,comp,'rows');
                if isempty(icomp)
                    disp(['Unknown component ' comp]);
                    return
                end
            else
                icomp = 1;
            end
            
            if nargin >= 3
                if iper<=0 || iper>nper
                    disp(['Unknown frequency # ' iper]);
                    return
                end
            else
                iper = 1;
            end
            
            k = iper;
            datareal = sqrt(info.per(k))*real(info.data(:,k,icomp));
            dataimag = sqrt(info.per(k))*imag(info.data(:,k,icomp));
            predreal = sqrt(info.per(k))*real(info.resp(:,k,icomp));
            predimag = sqrt(info.per(k))*imag(info.resp(:,k,icomp));
            diffreal = sqrt(info.per(k))*real(info.diff(:,k,icomp));
            diffimag = sqrt(info.per(k))*imag(info.diff(:,k,icomp));
            
            climsre(1) = (nanmin(datareal)+nanmin(predreal))/2;
            climsre(2) = (nanmax(datareal)+nanmax(predreal))/2;
            rangere = climsre(2) - climsre(1);
            if isnan(rangere)
                error('One of the data sets is entirely NaN for this period.');
            end

            climsim(1) = (nanmin(dataimag)+nanmin(predimag))/2;
            climsim(2) = (nanmax(dataimag)+nanmax(predimag))/2;
            rangeim = climsim(2) - climsim(1);
            if isnan(rangeim)
                error('One of the data sets is entirely NaN for this period.');
            end
            
            separateFigs = 0;
            if nargin > 3
                if ~ischar(fplot) && length(fplot)==6
                    separateFigs = 1;
                end
            end

            % Plot
            cper = sprintf('%03f',info.per(k));
            f=figure('Position',[300,300,1400,1000],...
                'PaperPosition',[1,1,14,10],...
                'PaperOrientation','Portrait');
            % plot the "data"
            %if ~separateFigs; sp=tight_subplot(3,2,[.1 .03],[.1 .1],[.03 .03]); end
            %if ~separateFigs; axes(sp(1)); end
            %if ~separateFigs; sp=tight_subplot(3,2,[.1 .03],[.1 .1],[.03 .03]); end
            if ~separateFigs; subplot(3,2,1); end
            str = ['Re ' comp ' x sqrt(T)'];
            InterpPlot(info.lon,info.lat,datareal,'',str,climsre);
            colormap redblue3;
            if separateFigs
                print('-djpeg95','-r300',fplot{1});
            end
            if ~separateFigs; subplot(3,2,2); end
            str = ['Im ' comp ' x sqrt(T)'];
            InterpPlot(info.lon,info.lat,dataimag,'',str,climsim);
            colormap redblue3;
            if separateFigs
                print('-djpeg95','-r300',fplot{2});
            end
            % plot the "resp"
            if ~separateFigs; subplot(3,2,3); end
            str = ['Re ' comp ' x sqrt(T)'];
            InterpPlot(info.lon,info.lat,predreal,'',str,climsre);
            colormap redblue3;
            if separateFigs
                print('-djpeg95','-r300',fplot{3});
            end
            if ~separateFigs; subplot(3,2,4); end
            str = ['Im ' comp ' x sqrt(T)'];
            InterpPlot(info.lon,info.lat,predimag,'',str,climsim);
            colormap redblue3;
            if separateFigs
                print('-djpeg95','-r300',fplot{4});
            end
            % plot the difference
            if ~separateFigs; subplot(3,2,5); end
            InterpPlot(info.lon,info.lat,diffreal,'','Difference',[-rangere/2 rangere/2]);
            colormap redblue3;
            if separateFigs
                print('-djpeg95','-r300',fplot{5});
            end
            if ~separateFigs; subplot(3,2,6); end
            InterpPlot(info.lon,info.lat,diffimag,'','Difference',[-rangeim/2 rangeim/2]);
            colormap redblue3;
            if separateFigs
                print('-djpeg95','-r300',fplot{6});
            end
            
            % save figure
            if nargin > 3
                print('-djpeg95','-r300',[fplot '_' comp '_' cper '.jpg']);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = plotrms(obj,maxrms,comp,iper,fplot)
            
            % [f] = plotrms(obj,maxrms,comp,iper,fplot)
            %
            % Plots a summary of RMS misfit.
            % If comp and iper are specified, plots that component
            % for one frequency, only.
            % If fplot is specified, saves a jpeg.
            
            %obj = mtdatacompare(obj1,obj2);
            %info = obj.v;
            [rms, info] = misfit(obj);
            nper = length(info.per);
            nsites = length(info.code);
            ncomp = length(info.comp);
            
            icomp = 0;
            if nargin > 2
                [~,icomp]=intersect(info.comp,comp,'rows');
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
            end
            
            if nargin > 3
                if iper<=0 || iper>nper
                    disp(['Unknown frequency # ' iper]);
                    return
                end
            else
                iper = 0;
            end

            % Average over frequency and component bins
            for j = 1:nper
                freq.f(j) = 1/(info.per(j));
                freq.p(j) = info.per(j);
                freq.d(j) = info.per(j)/(24*3600);
                for k = 1:ncomp
                    count = nsites - sum(isnan(info.err(:,j,k)));
                    freq.rms(j,k) = sqrt(nansum(info.res(:,j,k))/count);
                end
            end
            
            % Average over site bins
            if icomp > 0
                j = 1;
                for i = 1:nsites
                    count = nper - sum(isnan(info.res(i,:,icomp)));
                    if count > 0
                        site.codes(j,:) = info.code(i,:);
                        site.lon(j) = info.lon(i);
                        site.lat(j) = info.lat(i);
                        site.resp(j,:) = info.resp(i,:,icomp);
                        site.rms(j,:) = sqrt(info.res(i,:,icomp));
                        j = j+1;
                    end
                end
            else
                j = 1;
                for i = 1:nsites
                    count = ncomp*nper - sum(sum(isnan(info.res(i,:,:))));
                    if count > 0
                        site.codes(j,:) = info.code(i,:);
                        site.lon(j) = info.lon(i);
                        site.lat(j) = info.lat(i);
                        site.rms(j) = sqrt(nansum(nansum(info.res(i,:,:)))/count);
                        j = j+1;
                    end
                end
            end
            
            % Plot
            if iper > 0
                k = iper;
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
                if nargin > 1
                    caxis([1 maxrms]);
                end
                % xlabel('Log(10) T (secs)');
                % ylabel('RMS');
                % legend(h,'Zxx','Zxy','Zyx','Zyy','Tx','Ty','Location','BestOutside');
                if nargin > 4
                    print('-djpeg95','-r300',[fplot '_' comp '_' cper '.jpg']);
                end
            else
                figure('Position',[300,300,600,600],...
                    'PaperPosition',[1,1,6,6],...
                    'PaperOrientation','Portrait');
                subplot(7,1,1:5);
                htitle = ['RMS ' num2str(rms) ' (averaged over ' num2str(nper) ' periods)'];
                InterpPlot(site.lon,site.lat,site.rms,'RMS',htitle);
                colormap(invhot);
                if nargin > 1
                    caxis([1 maxrms]);
                end
                subplot(7,1,6:7);
                h = zeros(1,6);
                for i = 1:ncomp
                    color = [(i-1)/ncomp (i-1)/ncomp (i-1)/ncomp];
                    h(i)=semilogx(freq.p,freq.rms(:,i),'x-','color',color); hold on;
                    set(h(i),'LineWidth',2);
                end
                xlabel('Log(10) T (secs)');
                ylabel('RMS');
                legend(h,info.comp,'Location','BestOutside');
                if nargin > 4
                    print('-djpeg95','-r300',[fplot '_rms.jpg']);
                end
            end

        end
        
    end
end