function [h,hCB] = plotCond(m,ModelGrid,OPTIONS)
% Plots conductivity model
% Usage [h,hCB] = plotCond(m,grid,OPTIONS)
%   m is the 2D conductivity model to plot
%   grid is the model grid
%   OPTIONS is a structure of plotting OPTIONS
%          .nySkip = number of cells to omit from
%                    each end of the profile
%          .nZplot = number of vertical layers to plot
%                    (starting from the earth surface)
%          .siteLoc = site locations to plot (optional)
%          .ncax = color axis limits (vector with 2 elements)
%          .title = plot title

if strcmp(m.paramType,'LOGE')
   m.v = log10(exp(m.v));
else
   m.v = log10(m.v);
end

figure('Position',[100,100,600,400],...
	'PaperPosition',[1,1,6,4])
y = cumsum([0 ModelGrid.Dy])';
z = cumsum([0 ModelGrid.Dz])';
zCenter = z-z(ModelGrid.Nza+1);
%Nmid = fix((length(y)+1)/2);
%yCenter = y-y(Nmid);
zCenter = zCenter(ModelGrid.Nza+1:end)/1000;
yplot = y(OPTIONS.nYskip+1:end-OPTIONS.nYskip)/1000;
yCenter = yplot-mean(yplot);
Cond = real(m.v);
Cond = [Cond Cond(:,end)];
Cond = [Cond ;Cond(end,:)];
Cond = Cond(OPTIONS.nYskip+1:end-OPTIONS.nYskip,1:OPTIONS.nZplot);
h = pcolor(yCenter,zCenter(1:OPTIONS.nZplot),Cond'); ...
	shading flat; axis('ij');
set(gca,'FontWeight','demi','FontSize',13);
caxis(OPTIONS.cax); 
c = colmap;
c = c(end:-1:1,:);
X = [1:17];
XI = [1:.25:17];
c = interp1(X,c,XI);
colormap(c)
XI = [1:.25:17];
hCB = colorbar;

if isfield(OPTIONS,'siteLoc')
    siteLoc = OPTIONS.siteLoc/1000-mean(yplot);
    hold on; plot(siteLoc,0,'kv','markersize',6,'linewidth',2);
end

yt = [floor(OPTIONS.cax(1)):1:ceil(OPTIONS.cax(2))];
ytLabel = 10.^(-yt);

set(hCB,'FontWeight','demi','FontSize',12,...
   'Ytick',yt,'YtickLabel',ytLabel);
ylabel('Depth (km)');
xlabel('km');
title(OPTIONS.title,'interpreter','none');
