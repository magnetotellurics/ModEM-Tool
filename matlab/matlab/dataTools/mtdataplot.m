classdef mtdataplot < mtdata
    %   (c) Anna Kelbert, April 2014
    %
    %   Note as of Dec 2016: should also use mttf class for plotting
    
    methods
        
        function [obj] = mtdataplot(varargin)
            %   class constructor
            obj = obj@mtdata;
            
            if nargin == 0                
                return
            end
            
            % create this class from an mtdata object                               
            if isa(varargin{1},'mtdata')                
                mtdataobj = varargin{1};
                obj.header = mtdataobj.header;
                obj.d = mtdataobj.d;
                obj.nPeriods = mtdataobj.nPeriods;
                obj.nComp = mtdataobj.nComp;
                obj.Cmplx = mtdataobj.Cmplx;
                obj.units = mtdataobj.units;
                obj.signConvention = mtdataobj.signConvention;
                obj.orient = mtdataobj.orient;
                obj.compChar = mtdataobj.compChar;
                obj.primaryCoords = mtdataobj.primaryCoords;
                obj.type = mtdataobj.type;
                obj.predicted = mtdataobj.predicted;
                obj.origin = mtdataobj.origin;
                obj = obj.regroup;
            end
                            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = aprespltall(obj,fileName)
            % create and apparent resistivity and phase object from an
            % impedance object; can write the new object to file, or plot
                
            if isempty(obj.apres) || isempty(obj.phase)
                obj.regroup;
            end
            
            for i=1:length(obj.v.lon)
                siteName = char(obj.v.code(i,:));
                if nargin > 2
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
            % create and apparent resistivity and phase object from an
            % impedance object; can write the new object to file, or plot
                
            if isempty(obj.apres) || isempty(obj.phase)
                obj.regroup;
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
            p(1)=loglog(per,obj.apres.xy(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
            if ~obj.predicted
                errorbar(per,obj.apres.xy(i,:),2*obj.apres.xy_se(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5);
            end
            hold on
            p(2)=loglog(per,obj.apres.yx(i,:),'rx','linewidth',1.5); hold on
            if ~obj.predicted
                errorbar(per,obj.apres.yx(i,:),2*obj.apres.yx_se(i,:),'rx','linewidth',1.5);
            end
            hold off
            ymean = (nanmedian(obj.apres.xy(i,:))+nanmedian(obj.apres.yx(i,:)))/2;
            %set(gca,'xlim',[5 10^5],'ylim',[0.1 20000]);
            %set(gca,'xlim',[5 10^5],'ylim',[ymean/20 ymean*20]);
            set(gca,'ylim',[ymean/20 ymean*20]);
            set(gca,'xminortick','on','ytick',[1 10 100 1000 10000]);
            legend(p,'XY','YX');
            ylabel('\rho_a (\Omega m)','fontsize',14,'fontweight','demi');
            %RMS = sqrt(nansum(obj.v.res(i,:))/sum(~isnan(obj.v.err(i,:))));
            %title([char(OBS) ' [ LON = ' lon '; LAT = ' lat '; RMS = ' num2str(RMS)  ' ]'],'fontweight','demi');
            title([char(siteName) ' [ LON = ' lon '; LAT = ' lat ' ]'],'fontweight','demi');
            subplot(3,1,3);
            p(1)=semilogx(per,obj.phase.xy(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
            if ~obj.predicted
                errorbar(per,obj.phase.xy(i,:),2*obj.phase.xy_se(i,:),'o','color',[0. 0.5 0.8],'linewidth',1.5);
            end
            hold on
            p(2)=semilogx(per,obj.phase.yx(i,:),'rx','linewidth',1.5); hold on
            if ~obj.predicted
                errorbar(per,obj.phase.yx(i,:),2*obj.phase.yx_se(i,:),'rx','linewidth',1.5);
            end
            hold off
            ymean = (nanmean(obj.phase.xy(i,:))+nanmean(obj.phase.yx(i,:)))/2;
            %set(gca,'xlim',[5 10^5],'ylim',[-90 90]);
            %set(gca,'xlim',[5 10^5],'ylim',[ymean-45 ymean+45]);
            set(gca,'ylim',[ymean-45 ymean+45]);
            set(gca,'xminortick','on','yminortick','on');
            ylabel('\phi (Degrees)','fontsize',14,'fontweight','demi');

            if nargin > 2
                print(f,'-dpng',[fileName '_' siteName '_apres']);
            end
                        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = impplt(obj,siteName,fileName)
            % plot the transfer functions directly
            
            if isempty(obj.v)
                obj.regroup;
            end
            
            [~,i] = intersect(obj.v.code,siteName,'rows');
            if isempty(i)
                error(['No site in MT data object matching ' char(siteName)]);
            end
            
            lon = num2str(obj.v.lon(i));
            lat = num2str(obj.v.lat(i));
            per = obj.v.per;
            range = 100; %20;
            xlims = [0 1000]; %[5 10^5];
            ncomp = size(obj.d.data,3)/2;
            
            if strcmp(obj.type,'Full_Impedance')
                f=figure('Position',[200,200,1400,900],...
                    'PaperPosition',[1,1,18,10],...
                    'PaperOrientation','Portrait'); clf;
                % Diagonal TFs
                subplot(2,2,1)
                data = obj.v.data(i,:,1).*sqrt(per);
                std2 = 2*obj.v.err(i,:,1).*sqrt(per);
                TFinfo = [char(siteName) ' (' lon '; ' lat ')'];
                p(1)=semilogx(per,real(data),'bx','linewidth',1.5); hold on
                errorbar(per,real(data),std2,'bx','linewidth',1.5); hold on
                p(2)=semilogx(per,imag(data),'o','color',[0.8 0.5 0.],'linewidth',1.5); hold on
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
                std2 = 2*obj.v.err(i,:,4).*sqrt(per);
                TFinfo = [char(OBS) ' (' lon '; ' lat ')'];
                p(1)=semilogx(per,real(data),'bx','linewidth',1.5); hold on
                errorbar(per,real(data),std2,'bx','linewidth',1.5); hold on
                p(2)=semilogx(per,imag(data),'o','color',[0.8 0.5 0.],'linewidth',1.5); hold on
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
                std2 = 2*obj.v.err(i,:,2).*sqrt(per);
                TFinfo = [char(OBS) ' (' lon '; ' lat ')'];
                p(1)=semilogx(per,real(data),'rx','linewidth',1.5); hold on
                errorbar(per,real(data),std2,'rx','linewidth',1.5); hold on
                p(2)=semilogx(per,imag(data),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
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
                std2 = 2*obj.v.err(i,:,3).*sqrt(per);
                TFinfo = [char(OBS) ' (' lon '; ' lat ')'];
                p(1)=semilogx(per,real(data),'rx','linewidth',1.5); hold on
                errorbar(per,real(data),std2,'rx','linewidth',1.5); hold on
                p(2)=semilogx(per,imag(data),'o','color',[0. 0.5 0.8],'linewidth',1.5); hold on
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
                std2 = 2*obj.v.err(i,:,1);
                TFinfo = [char(OBS) ' (' lon '; ' lat ')'];
                p(1)=semilogx(per,real(data),'mx','linewidth',1.5); hold on
                errorbar(per,real(data),std2,'mx','linewidth',1.5); hold on
                p(2)=semilogx(per,imag(data),'go','linewidth',1.5); hold on
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
                std2 = 2*obj.v.err(i,:,2);
                TFinfo = [char(OBS) ' (' lon '; ' lat ')'];
                p(1)=semilogx(per,real(data),'mx','linewidth',1.5); hold on
                errorbar(per,real(data),std2,'mx','linewidth',1.5); hold on
                p(2)=semilogx(per,imag(data),'go','linewidth',1.5); hold on
                errorbar(per,imag(data),std2,'go','linewidth',1.5); hold off
                ymean = 0; %real(nanmedian(data));
                set(gca,'xlim',xlims,'ylim',[ymean-0.5 ymean+0.5]);
                set(gca,'xminortick','on','yminortick','on');
                legend(p,'Re(TY)','Im(TY)');
                ylabel('TY','fontweight','demi');
                xlabel('Period (secs)','fontweight','demi');
                title(TFinfo,'fontweight','demi');
            end
            print(f,'-dpng',[fileName '_fit_' siteName '_TFs']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = plot(obj,comp,iper,fplot)
            
            % [f] = plot(obj,comp,iper,fplot)
            %
            % Plots data as interpolated map.
            % Requires comp and iper.
            % If fplot is specified, saves a jpeg.
            %
            % For data comparison, use plot@mtdatacompare
                         
            if nargin < 2
                disp('Usage: [f] = obj.plot(comp,iper,fplot)');
                return
            end
            
            info = obj.v;
            isCmplx = obj.Cmplx;
            
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
            
            if isCmplx
                datareal = sqrt(info.per(k))*real(info.data(:,k,icomp));
                climsre(1) = (nanmin(rmoutliers(datareal)));
                climsre(2) = (nanmax(rmoutliers(datareal)));
                rangere = climsre(2) - climsre(1);
                strre = ['Re ' comp ' x sqrt(T)'];
                
                dataimag = sqrt(info.per(k))*imag(info.data(:,k,icomp));
                climsim(1) = (nanmin(rmoutliers(dataimag)));
                climsim(2) = (nanmax(rmoutliers(dataimag)));
                rangeim = climsim(2) - climsim(1);
                strim = ['Im ' comp ' x sqrt(T)'];
            else
                datareal = info.data(:,k,icomp);
                climsre(1) = (nanmin(rmoutliers(datareal)));
                climsre(2) = (nanmax(rmoutliers(datareal)));
                rangere = climsre(2) - climsre(1);
                strre = comp;
            end
            
            separateFigs = 1;
            
            plotReal = 1;
            plotImag = 1;

            % Plot
            cper = sprintf('%03f',info.per(k));
            f=figure('Position',[300,300,1400,1000],...
                'PaperPosition',[1,1,14,10],...
                'PaperOrientation','Portrait');
            
            if ~separateFigs && isCmplx && plotReal; subplot(2,1,1); end
            if plotReal
                InterpPlot(info.lon,info.lat,datareal,'',strre,climsre);
                colormap redblue3;
            end
            
            if separateFigs && isCmplx && plotImag
                f=figure('Position',[300,300,1400,1000],...
                    'PaperPosition',[1,1,14,10],...
                    'PaperOrientation','Portrait');
            end
            if ~separateFigs && isCmplx && plotImag; subplot(2,1,2); end
            if plotImag
                InterpPlot(info.lon,info.lat,dataimag,'',strim,climsim);
                colormap redblue3;
            end
            
            % save figure
            if nargin > 3
                print('-djpeg95','-r300',[fplot '_' comp '_' cper '.jpg']);
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = plotap(obj,comp,iper,fplot)
            
            % [f] = plotap(obj,comp,iper,fplot)
            %
            % Plots apparent resistivity & phase as interpolated map.
            % Requires comp and iper.
            % If fplot is specified, saves a jpeg.
            %
            % For data comparison, use plot@mtdatacompare
                         
            if nargin < 2
                disp('Usage: [f] = obj.plot(comp,iper,fplot)');
                return
            end
            
            info = obj.v;
            
            nper = length(info.per);

            if nargin >= 2
                if ~strcmp(lower(comp),'xy') && ~strcmp(lower(comp),'yx')
                    disp(['Unknown component ' comp]);
                    return
                end
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
            
                data1 = log10(obj.apres.(lower(comp))(:,k));
                clims1(1) = (nanmin(rmoutliers(data1)));
                clims1(2) = (nanmax(rmoutliers(data1)));
                range1 = clims1(2) - clims1(1);
                str1 = ['Ap. Res. ' comp];
                
                data2 = obj.phase.(lower(comp))(:,k);
                clims2(1) = (nanmin(rmoutliers(data2)));
                clims2(2) = (nanmax(rmoutliers(data2)));
                range2 = clims2(2) - clims2(1);
                str2 = ['Phase ' comp];
            
            separateFigs = 1;            
            
            plotApres = 1;
            plotPhase = 1;
            firstQuadrant = 0;
            
            if firstQuadrant; clims2(1) = 0; clims2(2) = 90; end
            if abs(clims2(2)-clims2(1))<0.01
                clims2(1) = clims2(1)-0.01;
                clims2(2) = clims2(2)+0.01;
            end

            % Plot
            cper = sprintf('%03f',info.per(k));
            %f=figure('Position',[300,300,2400,2000],...
            %    'PaperPosition',[1,1,24,16],...
            %    'PaperOrientation','Portrait');
            
            if separateFigs && plotApres; f(1)=figure; end
            if ~separateFigs && plotApres; f(1)=figure; subplot(2,1,1); end
            f(1).Position = [300,300,1800,1200];
            f(1).PaperPosition = [1,1,24,16];
            
            if plotApres
                InterpPlot(info.lon,info.lat,data1,'',str1,[0 4],1);
                colormap(flip(redblue3));
            end
            
            if separateFigs && plotPhase; f(2)=figure; end
            f(2).Position = [300,300,1800,1200];
            f(2).PaperPosition = [1,1,24,16];
            if ~separateFigs && plotPhase; subplot(2,1,2); end
            if plotPhase
                InterpPlot(info.lon,info.lat,data2,'',str2,clims2,1);
                colormap redblue3;
            end
            
            % save figure
            if nargin > 3
                print(f(1),'-djpeg95','-r300',[fplot '_' comp '_' cper '_apres.jpg']);
                if length(f)>1
                    print(f(2),'-djpeg95','-r300',[fplot '_' comp '_' cper '_phase.jpg']);
                end
                    
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [f] = pcolor(obj,siteName,fileName)
            % plot the map with colored dots
            % DRAFT ONLY DOES NOT WORK YET
            
            if isempty(obj.v)
                obj.regroup;
            end
            
            cmap = 'jet';
            cmap_dots = 'green';
            size_of_dots = 8;
            site_shape = 'o';
            if nargin >= 4
                if isfield(P,'cmap')
                    cmap = P.cmap;
                end
                if isfield(P,'sitecolor')
                    cmap_dots = P.sitecolor;
                end
                if isfield(P,'sitesize')
                    size_of_dots = P.sitesize;
                end
            end
            site_color = zeros(length(site.codes),3);
            site_round = zeros(length(site.codes),3);
            if (isa(site,'globalemdata') && ~isempty(site.rms)) || isfield(site,'rms')
                maxrms=3.0;
                if ~isempty(strfind(cmap_dots,'grey'))
                    map=cbrewer('seq','Greys',9);[mc,~]=size(map);
                elseif ~isempty(strfind(cmap_dots,'pink'))
                    map=cbrewer('seq','RdPu',9);[mc,~]=size(map);
                else
                    map=cbrewer('seq','YlGn',9);[mc,~]=size(map);
                end
                for k=1:length(site.rms)
                    ik=max(0,round(((site.rms(k)-1)/maxrms)*mc))+1;
                    ik=min(ik,mc);
                    site_color(k,:)=map(ik,:);
                    site_round(k,:) = [0. 0. 0.];
                end
            elseif strfind(cmap,'jet')
                for k=1:length(site.codes)
                    site_color(k,:) = [0.9 0.9 0.9];
                    site_round(k,:) = [0. 0. 0.];
                end
            else
                for k=1:length(site.codes)
                    site_color(k,:) = [0. 1. 0.]; %'g';
                    site_round(k,:) = [0. 0. 0.];
                end
            end
            if nargin >= 4
                if isfield(P,'site_color')
                    for k=1:length(site.codes)
                        site_color(k,:) = P.site_color;
                    end
                end
                if isfield(P,'site_round')
                    for k=1:length(site.codes)
                        site_round(k,:) = P.site_round;
                    end
                end
                if isfield(P,'site_shape')
                    site_shape = P.site_shape;
                end
                if isfield(P,'size_of_dots')
                    size_of_dots = P.size_of_dots;
                end
            end
            
            figure(h); hold on;
            for iobs=1:length(site.codes)
                m_plot(site.lon(iobs),site.lat(iobs),['w' site_shape],...
                    'markerfacecolor',site_color(iobs,:),...
                    'markeredgecolor',site_round(iobs,:),'markersize',size_of_dots);
                m_plot(site.lon(iobs)-360,site.lat(iobs),['w' site_shape],...
                    'markerfacecolor',site_color(iobs,:),...
                    'markeredgecolor',site_round(iobs,:),'markersize',size_of_dots);
            end
            
        end

    end
     
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = read(varargin)
            %
            % Usage:  [obj] = read(cfile,format,newUnits,dataType,predicted)
            %
            %  Reads an mtdata object formerly known as allData array
            
            mtdataobj = mtdata.read(varargin{:});
            
            obj = mtdataplot;
            obj = mtdataobj;
            
        end
    end

end