P = get(gcf,'UserData');
[Tx,Ty,X,Y] = getTipper(P);
X = X/1000; Y = Y/1000;
skip = 4;
h0 = figure('Position',[100,100,700,1000],...
        'PaperPosition',[1,1,5,8],...
        'PaperOrientation','Landscape');
axRect1 = [.2,.50,.7,.35];
axRect2 = [.2,.1,.7,.35];

%   script that plots real and imaginary parts of array C
%  plot real part
h1 = axes('Position',axRect1);h2 = axes('Position',axRect2);
axes(h1)
delta = Y(2)-Y(1);

h = quiver(Y(1:skip:end,1:skip:end),...
	X(1:skip:end,1:skip:end),...
	real(Ty(1:skip:end,1:skip:end))*skip*delta, ....
       real(Tx(1:skip:end,1:skip:end))*skip*delta,0);
set(h,'LineWidth',2);
set(gca,'FontWeight','demi','FontSize',16,...
	'Ylim',[40,100],'Xlim',[25,105])
text(.1,.9,'Real','Units','Normalized',...
        'FontWeight','demi','FontSize',16)
ylabel('x   (km)')
xlabel('y   (km)')
axis('xy')

%  plot imag part
axes(h2)
h = quiver(Y(1:skip:end,1:skip:end),...
	X(1:skip:end,1:skip:end),...
	imag(Ty(1:skip:end,1:skip:end))*skip*delta, ....
       imag(Tx(1:skip:end,1:skip:end))*skip*delta,0);
set(gca,'FontWeight','demi','FontSize',16, ...
	'Ylim',[40,100],'Xlim',[25,105])
text(.1,.9,'Imag','Units','Normalized',...
        'FontWeight','demi','FontSize',16)
set(h,'LineWidth',2);
ylabel('x   (km)')
xlabel('y   (km)')
axis('xy')

