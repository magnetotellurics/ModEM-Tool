function [h] = plotPseud(Z,OPTIONS);

rhoAx = [.15,.55,.8,.35];
phiAx = [.15,.1,.8,.35];
c = colmap;
c = c(end:-1:1,:);
X = [1:17];
XI = [1:.25:17];
c = interp1(X,c,XI);
c = c(end:-1:1,:);

[rho,phi,T,y,MODE,rho2,phi2,T2,y2] = mkPseud(Z);

h1 = figure('Position',[100,100,650,900],'PaperPosition',[1,1,6.5,9],...
     'PaperOrientation','Portrait')

%  plot apparent resistivity
h_rho = axes('Position',rhoAx,'FontWeight','demi','FontSize',14);
pcolor(y,log10(T),log10(rho));
axis('ij');
shading interp
colormap(c)
caxis(OPTIONS.rho_cax);
set(h_rho,'FontWeight','demi','FontSize',14);

hCB = colorbar;
yt = [floor(OPTIONS.rho_cax(1)):1:ceil(OPTIONS.rho_cax(2))];
ytLabel = 10.^(yt);
set(hCB,'FontWeight','demi','FontSize',12,...
   'Ytick',yt,'YtickLabel',ytLabel);

ylabel('log_{10} Period (s)');
xlabel('km');

ctitle = [ MODE ' Apparent Resistivity : ' OPTIONS.title];
title(ctitle);

%  plot phase
h_phi = axes('Position',phiAx,'FontWeight','demi','FontSize',14);
pcolor(y,log10(T),phi);
axis('ij')
shading interp
colormap(c)
caxis(OPTIONS.phi_cax);

hCB = colorbar;
yt = [floor(OPTIONS.phi_cax(1)):15:ceil(OPTIONS.phi_cax(2))]
ytLabel = yt;
set(hCB,'FontWeight','demi','FontSize',12,...
   'Ytick',yt,'YtickLabel',ytLabel);
set(h_phi,'FontWeight','demi','FontSize',14);

ylabel('log_{10} Period (s)');
xlabel('km');

ctitle = [ MODE ' Phase' ];
title(ctitle);

if MODE == 'JT'
   axes(h_rho)
   ctitle = ['TE Apparent Resistivity: ' OPTIONS.title ];
   title(ctitle);
   axes(h_phi)
   ctitle = ['TE Phase' ];
   title(ctitle);
   
   %  open a second figure for other mode
   h2 = figure('Position',[100,100,650,900],'PaperPosition',[1,1,6.5,9],...
     'PaperOrientation','Portrait')
   h = [h1 h2];

   h_rho = axes('Position',rhoAx,'FontWeight','demi','FontSize',14);
   pcolor(y2,log10(T2),log10(rho2));
   axis('ij')
   shading interp
   colormap(c)
   caxis(OPTIONS.rho_cax);

   hCB = colorbar;
   yt = [floor(OPTIONS.rho_cax(1)):1:ceil(OPTIONS.rho_cax(2))]
   ytLabel = 10.^(yt);
   set(hCB,'FontWeight','demi','FontSize',12,...
       'Ytick',yt,'YtickLabel',ytLabel);
   set(h_rho,'FontWeight','demi','FontSize',14);

   ylabel('log_{10} Period (s)');
   xlabel('km');

   ctitle = [ 'TM Apparent Resistivity : ' OPTIONS.title];
   title(ctitle);

   %  plot phase
   h_phi = axes('Position',phiAx,'FontWeight','demi','FontSize',14);
   pcolor(y2,log10(T2),phi2);
   axis('ij')
   shading interp
   colormap(c)
   caxis(OPTIONS.phi_cax);

   hCB = colorbar;
   yt = [floor(OPTIONS.phi_cax(1)):15:ceil(OPTIONS.phi_cax(2))]
   ytLabel = yt;
   set(hCB,'FontWeight','demi','FontSize',12,...
      'Ytick',yt,'YtickLabel',ytLabel);
   set(h_phi,'FontWeight','demi','FontSize',14);

   ylabel('log_{10} Period (s)');
   xlabel('km');

   ctitle = [ 'TM Phase' ];
   title(ctitle);
else
   h = h1;
end

