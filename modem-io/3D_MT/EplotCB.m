action = get(gcbo,'Tag');
h0 = get(gcbo,'Parent');
P = get(h0,'UserData');
NzAir = P.Grid.NzAir;
switch action
case 'Slice -'
  if(P.Np > 1)
    P.Np = P.Np-1;
    set(findobj('Tag','Plot Slice','Parent',h0),'Value',P.Np);
    set(h0,'UserData',P);
    if(P.PlotH)
       [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
    else
       [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
    end
    h1 = P.ax1; h2 = P.ax2;
    plotE
  end
case 'Slice +'
  if(P.PlotH)
     [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
  else
     [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
  end
  Nmax = length(lev);
  if(P.Np < Nmax)
    P.Np = P.Np+1;
    set(findobj('Tag','Plot Slice','Parent',h0),'Value',P.Np);
    set(h0,'UserData',P);
    if(P.PlotH)
       [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
    else
       [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
    end
    h1 = P.ax1; h2 = P.ax2;
    plotE
  end
case 'Plot Slice'
  Np = get(gcbo,'Value');
  P.Np = Np;
  set(h0,'UserData',P);
  if(P.PlotH)
     [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
  else
     [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
  end
  h1 = P.ax1; h2 = P.ax2;
  plotE

case 'X Slice'
  P.Slice = 'X';
  hX = gcbo;
  hP = findobj('Tag','Plot Slice','Parent',h0);
  set(findobj('Tag','Y Slice','Parent',h0),'Value',0);
  set(findobj('Tag','Z Slice','Parent',h0),'Value',0);
  [Nx1,Ny1,Nz1,Np] = CompLims(P);
  P.Np = Np;
  P.i2 = min(P.i2,Nx1);
  P.j2 = min(P.j2,Ny1);
  P.k2 = min(P.k2,Nz1);
  set(h0,'UserData',P)
  if(P.PlotH)
     [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
  else
     [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
  end
  Slices = {}
  for k = 1:Nx1
     Slices{k} = [ num2str(k) ' : ' num2str(lev(k),'%8.4e')] ;
  end
  set(hP,'Value',Np);
  set(hP,'String',Slices)
  h1 = P.ax1; h2 = P.ax2;
  plotE

case 'Y Slice'
  P.Slice = 'Y';
  hY = gcbo;
  hP = findobj('Tag','Plot Slice','Parent',h0);
  set(findobj('Tag','X Slice','Parent',h0),'Value',0);
  set(findobj('Tag','Z Slice','Parent',h0),'Value',0);
  [Nx1,Ny1,Nz1,Np] = CompLims(P);
  P.i2 = min(P.i2,Nx1);
  P.j2 = min(P.j2,Ny1);
  P.k2 = min(P.k2,Nz1);
  set(h0,'UserData',P)
  if(P.PlotH)
     [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
  else
     [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
  end
  Slices = {}
  for k = 1:Ny1
     Slices{k} = [ num2str(k) ' : ' num2str(lev(k),'%8.4e')] ;
  end
  set(hP,'Value',Np);
  set(hP,'String',Slices)
  h1 = P.ax1; h2 = P.ax2;
  plotE
  
case 'Z Slice'
  P.Slice = 'Z';
  hY = gcbo;
  hP = findobj('Tag','Plot Slice','Parent',h0);
  set(findobj('Tag','X Slice','Parent',h0),'Value',0);
  set(findobj('Tag','Y Slice','Parent',h0),'Value',0);
  [Nx1,Ny1,Nz1,Np] = CompLims(P);
  P.Np = Np;
  P.i2 = min(P.i2,Nx1);
  P.j2 = min(P.j2,Ny1);
  P.k2 = min(P.k2,Nz1);
  set(h0,'UserData',P)
  if(P.PlotH)
     [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
  else
     [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
  end
  Slices = {}
  for k = 1:Nz1
     Slices{k} = [ num2str(k) ' : ' num2str(lev(k),'%8.4e')];
  end
  set(hP,'Value',Np);
  set(hP,'String',Slices)
  h1 = P.ax1; h2 = P.ax2;
  plotE

case 'Limits'
  hP = gcbo;
  if(P.PlotH)
     [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
  else
     [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
  end
  [Nx1,Ny1,Nz1,Np] = CompLims(P);
  hLimits = figure('Position',[800,800,300,150], ...
	'MenuBar','none','NumberTitle','off',...
	'Name','Set Grid Limits','UserData',h0);
  
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.1,.75,.5,.15],...
        'String',['X limits : 1-' num2str(Nx1)],...
        'Tag','X limits',...
        'Style','text')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.55,.75,.2,.15],...
        'String',[num2str(P.i1)],...
        'Tag','X limits 1',...
        'Style','edit',...
        'CallBack','EsetLims')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.75,.75,.2,.15],...
        'String',[num2str(P.i2)],...
        'Tag','X limits 2',...
        'Style','edit',...
        'CallBack','EsetLims')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.1,.5,.5,.15],...
        'String',['Y limits : 1-' num2str(Ny1)],...
        'Tag','Y limits',...
        'Style','text')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.55,.5,.2,.15],...
        'String',[num2str(P.j1)],...
        'Tag','Y limits 1',...
        'Style','edit',...
        'CallBack','EsetLims')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.75,.5,.2,.15],...
        'String',[num2str(P.j2)],...
        'Tag','Y limits 2',...
        'Style','edit',...
        'CallBack','EsetLims')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.1,.25,.5,.15],...
        'String',['Z limits : 1-' num2str(Nz1)],...
        'Tag','Z limits',...
        'Style','text')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.55,.25,.2,.15],...
        'String',[num2str(P.k1)],...
        'Tag','Z limits 1',...
        'Style','edit',...
        'CallBack','EsetLims')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.75,.25,.2,.15],...
        'String',[num2str(P.k2)],...
        'Tag','Z limits 2',...
        'Style','edit',...
        'CallBack','EsetLims')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.10,.05,.20,.15],...
        'String',['NzAir=',num2str(NzAir)],...
        'Tag','NzAir',...
        'Style','text');
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.35,.05,.25,.15],...
        'String','OK',...
        'Tag','Limits OK',...
        'Style','pushbutton',...
        'CallBack','EsetLims')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.65,.05,.25,.15],...
        'String','Reset',...
        'Tag','Reset Limits',...
        'Style','pushbutton',...
        'CallBack','EsetLims')

case 'Components'
  COMPS = {'X','Y','Z'};
  P.Comp = COMPS{get(gcbo,'Value')};
  [Nx1,Ny1,Nz1,Np] = CompLims(P);
  P.Np = Np;
  P.i2 = min(P.i2,Nx1);
  P.j2 = min(P.j2,Ny1);
  P.k2 = min(P.k2,Nz1);
  set(h0,'UserData',P);
  if(P.PlotH)
     [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
  else
     [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
  end
  switch P.Slice
     case 'X'
        Nslice = Nx1;
     case 'Y'
        Nslice = Ny1;
     case 'Z'
        Nslice = Nz1;
  end
  Slices = {};
  for k = 1:Nslice
     Slices{k} = [ num2str(k) ' : ' num2str(lev(k),'%8.4e')] ;
  end
  set(findobj('Tag','Plot Slice','Parent',h0),'String',Slices);
  set(findobj('Tag','Plot Slice','Parent',h0),'Value',P.Np);
  h1 = P.ax1; h2 = P.ax2;
  plotE

case 'Mode and Period'
  nModePer = get(gcbo,'Value');
  nMode = get(gcbo,'UserData');
  if(nModePer ~= P.mode)
    %   read in new period/mode
    P.mode = nModePer;
    fNum = floor((nModePer-1)/nMode)+1;
    mNum = nModePer-nMode*(fNum-1);
    [P.E,T,Modes] = rdExyz(P.File, fNum, mNum);
    P.H = CurlE(P.E,P.Grid);
    P.PlotCurrents = 0;
    if(P.PlotH)
       plotWhat = 3
    else
       plotWhat = 1
    end
    set(findobj('Tag','Plot Currents'),'Value',plotWhat);
    set(h0,'UserData',P);
    if(P.PlotH)
       [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
    else
       [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
    end
    h1 = P.ax1; h2 = P.ax2;
    plotE
  end

case 'Fix Scale'
  clReal = get(P.ax1,'CLim');
  clImag = get(P.ax2,'CLim');
  hScale = figure('Position',[700,700,300,150], ...
	'MenuBar','none','NumberTitle','off',...
	'Name','Set Color Scale','UserData',h0);
  
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.1,.75,.5,.15],...
        'String',['Real '],...
        'Tag','Real limits',...
        'Style','text')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.55,.75,.2,.15],...
        'String',[num2str(clReal(1))],...
        'Tag','Real limits 1',...
        'Style','edit')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.75,.75,.2,.15],...
        'String',[num2str(clReal(2))],...
        'Tag','Real limits 2',...
        'Style','edit')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.1,.5,.5,.15],...
        'String',['Imaginary'],...
        'Tag','Imag limits',...
        'Style','text')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.55,.5,.2,.15],...
        'String',[num2str(clImag(1))],...
        'Tag','Imag limits 1',...
        'Style','edit')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.75,.5,.2,.15],...
        'String',[num2str(clImag(2))],...
        'Tag','Imag limits 2',...
        'Style','edit')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.25,.05,.25,.15],...
        'String','OK',...
        'Tag','Clim OK',...
        'Style','pushbutton',...
        'CallBack','EsetLims')
  uicontrol('Units','normalized',...
        'FontWeight','demi',...
        'Position',[.55,.05,.25,.15],...
        'String','Auto',...
        'Tag','Auto',...
        'Style','pushbutton',...
        'CallBack','EsetLims')
case 'Plot Currents'
  plotWhat = get(gcbo,'Value');
  PlotE = plotWhat==1;
  PlotH = plotWhat==3;
  P.PlotH = PlotH;
  [Nx1,Ny1,Nz1,Np] = CompLims(P);
  P.Np = Np;
  P.i2 = min(P.i2,Nx1);
  P.j2 = min(P.j2,Ny1);
  P.k2 = min(P.k2,Nz1);
  set(h0,'UserData',P);
  if(plotWhat == 2)
  %  currents
     if(~P.S.CondRead)
       %  if conductivity has not yet been read in,
       %    first need to read grid file
        Nchar = length(P.File);
        for k = Nchar:-1:1
          if(P.File(k) == '/')
             break
          end
        end
        filt = [P.File(1:k) '*'];
        [filename, pathname] = uigetfile(filt, 'Grid File');
        cfile = [pathname filename];
        [ModelGrid,rho] = rdModelRM(cfile);
        dzAir = 10*(3.^[1:ModelGrid.NzAir]);
        ModelGrid.dz = [dzAir(end:-1:1)' ; ModelGrid.dz];
        P.S = EdgeCond(ModelGrid,rho);
     end
     if(P.PlotCurrents == 0)
        P.PlotCurrents = 1;
        P.E.x = P.E.x.*P.S.x;
        P.E.y = P.E.y.*P.S.y;
        P.E.z = P.E.z.*P.S.z;
        set(h0,'UserData',P)
     end 
     [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
     switch P.Slice
        case 'X'
           Nslice = Nx1;
        case 'Y'
           Nslice = Ny1;
        case 'Z'
           Nslice = Nz1;
     end
     Slices = {}
     for k = 1:Nslice
        Slices{k} = [ num2str(k) ' : ' num2str(lev(k),'%8.4e')] ;
     end
     set(findobj('Tag','Plot Slice','Parent',h0),'String',Slices);
     set(findobj('Tag','Plot Slice','Parent',h0),'Value',P.Np);
     h1 = P.ax1; h2 = P.ax2;
     plotE
  elseif(plotWhat==1)
  %  E-fields
     if(P.PlotCurrents == 1)
        P.PlotCurrents = 0;
        P.E.x = P.E.x./P.S.x;
        P.E.y = P.E.y./P.S.y;
        P.E.z = P.E.z./P.S.z;
        set(h0,'UserData',P)
     end
     [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
     switch P.Slice
        case 'X'
           Nslice = Nx1;
        case 'Y'
           Nslice = Ny1;
        case 'Z'
           Nslice = Nz1;
     end
     Slices = {}
     for k = 1:Nslice
        Slices{k} = [ num2str(k) ' : ' num2str(lev(k),'%8.4e')] ;
     end
     set(findobj('Tag','Plot Slice','Parent',h0),'String',Slices);
     set(findobj('Tag','Plot Slice','Parent',h0),'Value',P.Np);
     h1 = P.ax1; h2 = P.ax2;
     plotE
  elseif(plotWhat == 3)
  %  H-Fields
     [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
     switch P.Slice
        case 'X'
           Nslice = Nx1;
        case 'Y'
           Nslice = Ny1;
        case 'Z'
           Nslice = Nz1;
     end
     Slices = {}
     for k = 1:Nslice
        Slices{k} = [ num2str(k) ' : ' num2str(lev(k),'%8.4e')] ;
     end
     set(findobj('Tag','Plot Slice','Parent',h0),'String',Slices);
     set(findobj('Tag','Plot Slice','Parent',h0),'Value',P.Np);
     h1 = P.ax1; h2 = P.ax2;
     plotE
  end
end
