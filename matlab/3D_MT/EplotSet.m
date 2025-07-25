function [h0] = EplotSet(E,S,Grid,options)
% Usage EplotSet(E,Grid,options)
%
%  plots electric field components; Grid is a structure containing dx,dy,dz,
%   E a structure containg Ex, Ey Ez, and options gives initial plot options


Nx = size(Grid.dx);
Ny = size(Grid.dy);
Nz = size(Grid.dz);
h0 = figure('Position',[100,100,700,1000],...
	'PaperPosition',[1,1,5,8],...
	'PaperOrientation','Landscape');
axRect1 = [.2,.50,.7,.35];
axRect2 = [.2,.1,.7,.35];
h1 = 0; h2 = 0;

%   get default plotting options from structure "options"
% Comp is E-field component to plot   
Comp = options.Comp
Components = {'x','y','z'};
if(Comp == 'X')
   nComp = 1;
elseif(Comp == 'Y')
   nComp = 2;
else
   nComp = 3;
end
  
% Slice is slice to display (X, Y or Z = constant)
Slice = options.slice;
% Np = number of slice to plot initially
Np = options.Np;
%   T is list of periods
T = options.T;
%   nPer is period number
nPer = options.nPer;
%   mode is mode (Ex or Ey)
jk=0;
nMode = length(options.Modes);
for k = 1:nPer
   for j = 1:nMode
      jk = jk+1;
      PeriodMode{jk} = [options.Modes{j}  ' : T=' num2str(T(k))];
   end
end
fNum = (options.iPer-1)*2 + options.mode;

%  i1:i2, etc. give limits x, y , z coordinates to plot
i1 = options.iXlim(1);
i2 = options.iXlim(2);
j1 = options.iYlim(1);
j2 = options.iYlim(2);
k1 = options.iZlim(1);
k2 = options.iZlim(2);

NzAir = Grid.NzAir;
xEdge = Grid.origin(1)+[0 ; cumsum(Grid.dx)];
yEdge = Grid.origin(2)+[0 ; cumsum(Grid.dy)];
zEdge = Grid.origin(3)+[0 ; cumsum(Grid.dz)];
zSurface = zEdge(NzAir+1);
zEdge = zEdge-zSurface;
del = [ Grid.dx(1); Grid.dx(1:end-1)+Grid.dx(2:end) ]/2;
xCenter = Grid.origin(1)+cumsum(del);
del = [ Grid.dy(1); Grid.dy(1:end-1)+Grid.dy(2:end) ]/2;
yCenter = Grid.origin(2)+cumsum(del);
del = [ Grid.dz(1); Grid.dz(1:end-1)+Grid.dz(2:end) ]/2;
zCenter = Grid.origin(3)+cumsum(del)-zSurface;

H = CurlE(E,Grid);
H.x = H.x/(1i*E.omega);
H.y = H.y/(1i*E.omega);
H.z = H.z/(1i*E.omega);
%   store all information (including data plotted) in data structure P
P = struct('E',E,'H',H,'xEdge',xEdge,'yEdge',yEdge,'zEdge',zEdge,...
	'xCenter',xCenter,'yCenter',yCenter,'zCenter',zCenter,...
	'ax1',h1,'ax2',h2,'Comp',Comp,'T',T,...
	'nPer',nPer,'mode',options.mode,...
	'Slice',Slice,'Np',Np,'i1',i1,'i2',i2,'j1'...
	,j1,'j2',j2,'k1',k1,'k2',k2,'File',options.SolnFile,...
	'CLmode','Auto','CLreal',[],'CLimag',[],'S',S,...
	'PlotCurrents',0,'PlotH',0,'Grid',Grid,'PlotType','E_H');

%  now set up component to plot
[C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
h1 = axes('Position',axRect1);h2 = axes('Position',axRect2);
P.ax1 = h1; P.ax2 = h2;

%   associate  all of the data in structure P with figure
set(h0,'UserData',P)

plotE

for k = 1:length(lev)
   Slices{k} = [ num2str(k) ' : ' num2str(lev(k),'%8.4e')] ;
end

uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.125,.95,.05,.030],...
	'String','-',...
	'Tag','Slice -',...
	'Style','pushbutton',...
	'CallBack','EplotCB')
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.175,.95,.05,.030],...
	'String','+',...
	'Tag','Slice +',...
	'Style','pushbutton',...
	'CallBack','EplotCB')
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.225,.95,.20,.030],...
	'String',Slices,...
	'Tag','Plot Slice',...
	'Style','popupmenu',...
	'CallBack','EplotCB',...
	'Value',Np);
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.425,.95,.10,.03],...
	'String','X',...
	'Tag','X Slice',...
	'Style','radiobutton',...
	'CallBack','EplotCB',...
	'Value',0);
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.525,.95,.10,.03],...
	'String','Y',...
	'Tag','Y Slice',...
	'Style','radiobutton',...
	'CallBack','EplotCB',...
	'Value',0);
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.625,.95,.10,.03],...
	'String','Z',...
	'Tag','Z Slice',...
	'Style','radiobutton',...
	'CallBack','EplotCB',...
	'Value',1);
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.725,.95,.10,.03],...
	'String','Limits',...
	'Tag','Limits',...
	'Style','pushbutton',...
	'CallBack','EplotCB')
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.825,.95,.10,.03],...
	'String',Components,...
	'Tag','Components',...
	'Style','popupmenu',...
	'CallBack','EplotCB',...
	'Value',nComp);
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.500,.91,.25,.030],...
	'String',PeriodMode,...
	'Tag','Mode and Period',...
        'UserData',nMode,...
	'Style','popupmenu',...
        'Value',fNum,...
	'CallBack','EplotCB')
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.350,.91,.15,.030],...
	'String','Mode/Period',...
	'Tag','Period Text',...
	'Style','text')
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.750,.91,.15,.03],...
	'String','Fix Color Scale',...
	'Tag','Fix Scale',...
	'Style','pushbutton',...
	'CallBack','EplotCB',...
	'Value',1);
uicontrol('Units','normalized',...
	'FontWeight','demi',...
	'Position',[.200,.91,.14,.03],...
	'String',{'E-fields','Currents','H-fields'},...
	'Tag','Plot Currents',...
	'Style','popupmenu',...
	'CallBack','EplotCB',...
	'Value',1);
