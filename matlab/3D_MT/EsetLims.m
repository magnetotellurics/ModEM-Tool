action = get(gcbo,'Tag');
%  Handle for main figure which generated this subwindow
hWindow = get(gcbo,'Parent');
h0 = get(hWindow,'UserData');
P = get(h0,'UserData');
switch action
case 'X limits 1'
   P.i1 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'X limits 2'
   P.i2 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Y limits 1'
   P.j1 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Y limits 2'
   P.j2 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Z limits 1'
   P.k1 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Z limits 2'
   P.k2 = str2num(get(gcbo,'String'));
   set(h0,'UserData',P);
case 'Limits OK'
   if(P.PlotType == 'E_H')
      if(P.PlotH)
         [C,xLab,yLab,X,Y,AxMode,lev] = getHcomp(P);
      else
         [C,xLab,yLab,X,Y,AxMode,lev] = getEcomp(P);
      end
      h1 = P.ax1; h2 = P.ax2;
      plotE
   elseif(P.PlotType == 'TFs')
      [C,xLab,yLab,X,Y,AxMode,lev] = getZcomp(P);
      h1 = P.ax1; h2 = P.ax2;
      plotAmpPh
   end
   close(hWindow)
case 'Reset Limits'
  if(P.PlotType == 'E_H')
     [Nx1,Ny1,Nz1,Np] = CompLims(P);
  else
     Nx1 = P.nX;
     Ny1 = P.nY;
     Nz1 = P.nPer;
  end
  P.i1 = 1;
  P.i2 = Nx1; 
  P.j1 = 1;
  P.j2 = Ny1; 
  P.k1 = 1;
  P.k2 = Nz1; 
  set(h0,'UserData',P);
  set(findobj('Tag','X limits 1'),'String',num2str(P.i1));
  set(findobj('Tag','X limits 2'),'String',num2str(P.i2));
  set(findobj('Tag','Y limits 1'),'String',num2str(P.j1));
  set(findobj('Tag','Y limits 2'),'String',num2str(P.j2));
  set(findobj('Tag','Z limits 1'),'String',num2str(P.k1));
  set(findobj('Tag','Z limits 2'),'String',num2str(P.k2));
case 'Cancel'
   close(hWindow)
%case 'Real limts 1'
%   clReal(1) = str2num(get(gcbo,'String'));
%case 'Real limts 2'
%   clReal(2) = str2num(get(gcbo,'String'));
%case 'Imag limts 1'
%   clImag(1) = str2num(get(gcbo,'String'));
%case 'Imag limts 2'
%   clImag(2) = str2num(get(gcbo,'String'));
case 'Clim OK'
   clReal(1) = str2num(get(findobj('Tag','Real limits 1',...
	'Parent',hWindow),'String'));
   clReal(2) = str2num(get(findobj('Tag','Real limits 2',...
	'Parent',hWindow),'String'));
   clImag(1) = str2num(get(findobj('Tag','Imag limits 1',...
	'Parent',hWindow),'String'));
   clImag(2) = str2num(get(findobj('Tag','Imag limits 2',...
	'Parent',hWindow),'String'));
   P.CLmode = 'manu';
   P.CLreal = clReal;
   P.CLimag = clImag;
   set(h0,'UserData',P);
   set(P.ax1,'Clim',clReal);
   axes(P.ax1); hcb = colorbar;
   set(hcb,'FontWeight','demi')
   set(P.ax2,'Clim',clImag);
   axes(P.ax2); hcb = colorbar;
   set(hcb,'FontWeight','demi')
   close(hWindow)
case 'Auto'
   P.CLmode = 'auto';
   set(h0,'UserData',P);
   close(hWindow)
case 'Save CL'
   clReal(1) = str2num(get(findobj('Tag','Real limits 1',...
	'Parent',hWindow),'String'));
   clReal(2) = str2num(get(findobj('Tag','Real limits 2',...
	'Parent',hWindow),'String'));
   clImag(1) = str2num(get(findobj('Tag','Imag limits 1',...
	'Parent',hWindow),'String'));
   clImag(2) = str2num(get(findobj('Tag','Imag limits 2',...
	'Parent',hWindow),'String'));
   [file,path] = uiputfile('*.clm');
   cfile = [path file];
   eval(['save ' cfile ' clReal clImag']);
case 'Load CL'
   [file,path] = uigetfile('*.clm');
   eval(['load ' path file ' -MAT']);
   set(findobj('Tag','Real limits 1', ...
		'Parent',hWindow),'String',num2str(clReal(1)));
   set(findobj('Tag','Real limits 2', ...
		'Parent',hWindow),'String',num2str(clReal(2)));
   set(findobj('Tag','Imag limits 1', ...
		'Parent',hWindow),'String',num2str(clImag(1)));
   set(findobj('Tag','Imag limits 2', ...
		'Parent',hWindow),'String',num2str(clImag(2)));
end
