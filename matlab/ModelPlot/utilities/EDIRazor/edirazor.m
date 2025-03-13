function varargout = edirazor(varargin)
% edirazor - A GUI tool to remove bad data in EDI file.
%
%  NOTES: 
%  
%
%  See also MTDataRazor.
%  --------------------
%  Bo Yang, 2014.
%  Institute of Geophysics and Geomatics,
%  China University of Geosciences, Wuhan, China.
%  Comments, bug reports and questions, please send to:
%  yangbo.cug@163.com.
%  Copyright 2014-2018, Bo Yang, IGG, CUG.
%  $Revision: 1.0 $ $Date: 2014/02/18 $
%  Last changed: 2014/03/28 12:45:11.

%  Revision log:
%  2014/02/18 : Version 1.0 released.
%  2014/03/26 : Version 1.1 released.
%  2014/03/27 : Version 1.2 released. Modified dispedi.

global edi
clear edi
%
%  define the GUI.
%
hMainGUIFig = figure(...       % the main GUI figure.	'Toolbar','none', ...	'MenuBar','none', ...
	                'HandleVisibility','callback', ...
	                'Name', 'EDIRazor 1.0', ...
	                'NumberTitle','off', ...
	                'units','normalized',...
	                'tag','tMainGUIFig',...
	                'position',[0.01,0.06,0.98,0.80],...
	                'Color', get(0, 'defaultuicontrolbackgroundcolor'));
%---------------------------define the main GUIs------------------------
%
%  define the btn size.
%
	btnr = 5;
	btnw = 0.1;
	btnh = 0.1;
	axew = 0.7;
	lisw = 0.14;
% 	hAxesDisp   =    axes(...         % the axes for displaying data.
% 		                'Parent', hMainGUIFig, ...
% 		                'Units', 'normalized', ...
% 		                'DrawMode','fast',...
% 		                'Visible','on',...
% 		                'HandleVisibility','callback', ...
% 		                'tag','tAxesDisp',...
% 		                'Position',[0.02 0.05 axew 0.9]);
	h0   =    axes(...         % the axes for displaying data.
		                'Parent', hMainGUIFig, ...
		                'Units', 'normalized', ...
		                'DrawMode','fast',...
		                'Visible','off',...
		                'tag','tAxesh0',...
		                'Position',[0.03,0.05,0.7,0.90]);
	h1   =    axes(...         % the axes for displaying data.
		                'Parent', hMainGUIFig, ...
		                'Units', 'normalized', ...
		                'DrawMode','fast',...
		                'Visible','off',...
		                'tag','tAxesh1',...
		                'Position',[0.03,0.45,0.7,0.50]);
	h2   =    axes(...         % the axes for displaying data.
		                'Parent', hMainGUIFig, ...
		                'Units', 'normalized', ...
		                'DrawMode','fast',...
		                'Visible','off',...
		                'tag','tAxesh2',...
		                'Position',[0.03,0.05,0.7,0.40]);
	hEdiList  = uicontrol(...
		                'Parent', hMainGUIFig, ...
		                'Units', 'normalized', ...
		                'Style','listbox',...
		                'tag','tEdiList',...
		                'callback',@tEdiList_Callback,...
		                'Position',[0.76 0.05 lisw 0.9]);
	colname     = {'','Comp','Visible','Selectable'};
	colfmt      = {'char','char','logical','logical'};
	coleditable = [false false true true];
	colwidth    = {10,40,40,40};
	comptab     = {'','Rhoxy',true, true;'','Rhoyx',true, true;...
		           '','Phsxy',true,true;'','Phsyx',true,true;...
		           '','ZXXR',false,false;'','ZXXI',false,false;...
		           '','ZXYR',false,false;'','ZXYI',false,false;...
		           '','ZYXR',false,false;'','ZYXI',false,false;...
		           '','ZYYR',false,false;'','ZYYI',false,false;...
		           '','TXR', false,false;'','TXI', false,false;...
		           '','TYR', false,false;'','TYI', false,false;...
		           };
	hTabComp   = uitable(... % the table of properties of component.
		                'Parent',hMainGUIFig,...
		                'Units','Normalized',...
		                'tag','tTabComp',...
		                'Position',[0.90,0.25,btnw,0.70],...
		                'ColumnName',colname,...
		                'ColumnFormat',colfmt,...
		                'ColumnEditable',coleditable,...
		                'ColumnWidth',colwidth,...
		                'Data',comptab,...
		                'CellEditCallback',@tTabComp_CellEditCallback,...
		                'CellSelectionCallback',@tTabComp_CellSelectionCallback,...
		                'RowName',[]);
	hBtnLoadList = uicontrol(... % button for loading edi list.
		                'Parent',hMainGUIFig,...
		                'Units','Normalized',...
		                'HandleVisibility','callback',...
		                'Position',[0.90 0.15 btnw btnh],...
		                'FontSize',10,...
		                'tag','tBtnLoadList',...
		                'Enable','on',...
		                'String','Load list',...
		                'Callback',@tBtnLoadList_Callback);
	hBtnSaveAs = uicontrol(... % button for saving edi file.
		                'Parent',hMainGUIFig,...
		                'Units','Normalized',...
		                'HandleVisibility','callback',...
		                'Position',[0.90 0.05 btnw btnh],...
		                'FontSize',10,...
		                'tag','tBtnSaveAs',...
		                'Enable','on',...
		                'String','Save FSEL',...
		                'Callback',@tBtnSaveAs_Callback);
%
%  set all msg response functions.
%
	set(hMainGUIFig,'WindowButtonDownFcn',{@ButtonDown_Callback});
%
%  save gui data.
%
	handles = guihandles(hMainGUIFig);
	guidata(hMainGUIFig,handles);
%
% End of the function.
%
%----------------------------------------------------------------------%
%                              btn callbacks                           %
%----------------------------------------------------------------------%

%----------------------------------------------------------------------%
	function tBtnLoadList_Callback(hObject,eventdata,handles)
		global edilist edi dispcomp selcomp
		handles = guidata(gcbo);
		[file,path] = uigetfile({'*.list';'*.*'},'Open a edi list file');
		% read the list.
		if (file ~= 0)
			fid = fopen([path file]);
			k = 0;
			while 1
				tline = fgetl(fid);
				if ischar(tline)
					k = k + 1;
					edilist{k} = tline;
				else
					break;
				end % check file name.
			end % while.
			fclose(fid);
			set(handles.tEdiList,'String',edilist);
			edi = loadedi(edilist);
			tab = get(handles.tTabComp,'Data');
			dispcomp = tab(:,3);
			selcomp  = tab(:,4);
			% initialize the axis with plotyy.
			dispedi(handles,1);
		else
			return;
		end % check list file name.

%----------------------------------------------------------------------%
	function tBtnSaveAs_Callback(hObject,eventdata,handles)
		global edi
		for k = 1:length(edi);
% 			[path,name,ext] = fileparts(edi(k).filename);
% 			filename = [path,'/',name,'_fsel',ext];
% 			edi(k).write(filename);
			edi(k).writefsel();
			disp(['FSEL file for: ',edi(k).filename,' saved!']);
		end % k

%----------------------------------------------------------------------%
%                            msg callbacks                             %
%----------------------------------------------------------------------%
	function tEdiList_Callback(hObject,eventdata,handles)
		handles = guidata(gcbo);
		iSel = get(handles.tEdiList,'Value');
		dispedi(handles,iSel);

%----------------------------------------------------------------------%
	function tTabComp_CellEditCallback(hObject,eventdata,handles)
		global dispcomp selcomp
		handles = guidata(gcbo);
		tab = get(handles.tTabComp,'Data');
		dispcomp = tab(:,3);
		selcomp  = tab(:,4);
		iSel = get(handles.tEdiList,'Value');
		dispedi(handles,iSel);

%----------------------------------------------------------------------%
	function tTabComp_CellSelectionCallback(hObject,eventdata,handles)
	
%----------------------------------------------------------------------%
	function ButtonDown_Callback(h, eventdata, handles, varargin)
		global edi
		handles = guidata(gcbo);
		iSel = get(handles.tEdiList,'Value');
		bIsInAxes = IsInAxes;
		if bIsInAxes ~= 1
			return
		else
			CP = get(gca,'CurrentPoint');
			x = CP(1,1);
			freq = edi(iSel).freq;
			kf = getindex(x,freq);
			ms = get(get(gca,'parent'),'SelectionType');
			idx = getselidx;
			if strcmp(ms,'normal') %LBtn down
				edi(iSel).asel(kf,idx) = false;
			elseif strcmp(ms,'alt') %RBtn down
				edi(iSel).asel(kf,idx) = true;
			end
		end
		% draw edi.
		dispedi(handles,iSel);

%----------------------------------------------------------------------%
%                          general functions                           %
%----------------------------------------------------------------------%
	function idx = getselidx()
		global selcomp
		idx(1) = selcomp{1}; % rhoxy
		idx(2) = selcomp{2}; % rhoyx
		idx(3) = selcomp{3}; % phsxy
		idx(4) = selcomp{4}; % phsyx
		idx(5) = selcomp{5 } | selcomp{6 }; % zxx
		idx(6) = selcomp{7 } | selcomp{8 }; % zxy
		idx(7) = selcomp{9 } | selcomp{10}; % zyx
		idx(8) = selcomp{11} | selcomp{12}; % zyy
		idx(9) = selcomp{13} | selcomp{14}; % tzx
		idx(10)= selcomp{15} | selcomp{16}; % tzx

%----------------------------------------------------------------------%
	function dispedi(handles,k)
		global dispcomp edi
		edi(k) = edi(k).asel2csel;
		vdisp = cell2mat(dispcomp);
		%
		% loglog subfigure.
		%
		if (sum(vdisp(1:2))+sum(vdisp(5:12))) * (sum(vdisp(3:4))+sum(vdisp(13:16))) > eps
			axes(handles.tAxesh1);
			set(handles.tAxesh1,'Visible','on');
			set(handles.tAxesh2,'Visible','on');
			set(handles.tAxesh0,'Visible','off');
		else
			axes(handles.tAxesh0);
			set(handles.tAxesh1,'Visible','off');
			set(handles.tAxesh2,'Visible','off');
			set(handles.tAxesh0,'Visible','on');
		end % if disp loglog & semilogx comp.
		if dispcomp{1}
			sel = edi(k).rhoxys > 0;
			loglog(edi(k).freq,edi(k).rhoxy,'b-');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).rhoxy(sel),edi(k).drhoxy(sel),'bo');
			hold on;
		end
		if dispcomp{2}
			sel = edi(k).rhoyxs > 0;
			loglog(edi(k).freq,edi(k).rhoyx,'b--');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).rhoyx(sel),edi(k).drhoyx(sel),'bs');
			hold on;
		end
		if dispcomp{5}
			sel = edi(k).zxxs > 0;
			loglog(edi(k).freq,edi(k).zxxr,'r-.');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zxxr(sel),edi(k).dzxx(sel),'r^');
			hold on;
		end
		if dispcomp{6}
			sel = edi(k).zxxs > 0;
			loglog(edi(k).freq,edi(k).zxxi,'r-.');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zxxi(sel),edi(k).dzxx(sel),'rv');
			hold on;
		end
		if dispcomp{7}
			sel = edi(k).zxys > 0;
			loglog(edi(k).freq,edi(k).zxyr,'r-');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zxyr(sel),edi(k).dzxy(sel),'r^');
			hold on;
		end
		if dispcomp{8}
			sel = edi(k).zxys > 0;
			loglog(edi(k).freq,edi(k).zxyi,'r-');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zxyi(sel),edi(k).dzxy(sel),'rv');
			hold on;
		end
		if dispcomp{9}
			sel = edi(k).zyxs > 0;
			loglog(edi(k).freq,edi(k).zyxr,'r--');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zyxr(sel),edi(k).dzyx(sel),'r^');
			hold on;
		end
		if dispcomp{10}
			sel = edi(k).zyxs > 0;
			loglog(edi(k).freq,edi(k).zyxi,'r--');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zyxi(sel),edi(k).dzyx(sel),'rv');
			hold on;
		end
		if dispcomp{11}
			sel = edi(k).zyys > 0;
			loglog(edi(k).freq,edi(k).zyyr,'r:');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zyyr(sel),edi(k).dzyy(sel),'r^');
			hold on;
		end
		if dispcomp{12}
			sel = edi(k).zyys > 0;
			loglog(edi(k).freq,edi(k).zyyi,'r:');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zyyi(sel),edi(k).dzyy(sel),'rv');
			hold on;
		end
		grid on
		set(gca,'GridLineStyle','-','XMinorGrid','off','YMinorGrid','off');
		set(gca,'Xcolor',[0.5 0.5 0.5]);
		set(gca,'Ycolor',[0.5 0.5 0.5]);
		set(gca,'xdir','reverse');
		if (sum(vdisp(1:2))+sum(vdisp(5:12))) * (sum(vdisp(3:4))+sum(vdisp(13:16))) > eps
			set(gca,'xticklabel',[]);
		end % if disp loglog & semilogx comp.
		hold off
		%
		% semilogx subfigure.
		%
		if (sum(vdisp(1:2))+sum(vdisp(5:12))) * (sum(vdisp(3:4))+sum(vdisp(13:16))) > eps
			axes(handles.tAxesh2);
		end % if disp loglog & semilogx comp.
		if dispcomp{3}
			sel = edi(k).phsxys > 0;
			semilogx(edi(k).freq,edi(k).phsxy,'b-');
			hold on;
			errorbare('vlogx',edi(k).freq(sel),edi(k).phsxy(sel),edi(k).dphsxy(sel),'bo');
			hold on;
		end
		if dispcomp{4}
			sel = edi(k).phsyxs > 0;
			semilogx(edi(k).freq,edi(k).phsyx,'b--');
			hold on;
			errorbare('vlogx',edi(k).freq(sel),edi(k).phsyx(sel),edi(k).dphsyx(sel),'bs');
			hold on;
		end
		if dispcomp{13}
			sel = edi(k).tzxs > 0;
			semilogx(edi(k).freq,edi(k).tzxr,'m-');
			hold on;
			errorbare('vlogx',edi(k).freq(sel),edi(k).tzxr(sel),edi(k).dtzx(sel),'m^');
			hold on;
		end
		if dispcomp{14}
			sel = edi(k).tzxs > 0;
			semilogx(edi(k).freq,edi(k).tzxi,'m-');
			hold on;
			errorbare('vlogx',edi(k).freq(sel),edi(k).tzxi(sel),edi(k).dtzx(sel),'mv');
			hold on;
		end
		if dispcomp{15}
			sel = edi(k).tzys > 0;
			semilogx(edi(k).freq,edi(k).tzyr,'m--');
			hold on;
			errorbare('vlogx',edi(k).freq(sel),edi(k).tzyr(sel),edi(k).dtzy(sel),'m^');
			hold on;
		end
		if dispcomp{16}
			sel = edi(k).tzys > 0;
			semilogx(edi(k).freq,edi(k).tzyi,'m--');
			hold on;
			errorbare('vlogx',edi(k).freq(sel),edi(k).tzyi(sel),edi(k).dtzy(sel),'mv');
			hold on;
		end
		grid on
		set(gca,'GridLineStyle','-','XMinorGrid','off','YMinorGrid','off');
		set(gca,'Xcolor',[0.5 0.5 0.5]);
		set(gca,'Ycolor',[0.5 0.5 0.5]);
		set(gca,'xdir','reverse');
		hold off;

%----------------------------------------------------------------------%
	function dispedi_old(handles,k)
		global dispcomp edi
		legstr = cell(1);
		kk = 0;
		%sel = edi(k).fsel > 0;
		edi(k) = edi(k).asel2csel;
		if dispcomp{1}
			sel = edi(k).rhoxys > 0;
			loglog(edi(k).freq,edi(k).rhoxy,'b');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).rhoxy(sel),edi(k).drhoxy(sel),'bo');
			hold on;
			kk = kk + 1; legstr{kk} = 'rhoxy';
		end
		if dispcomp{2}
			sel = edi(k).rhoyxs > 0;
			loglog(edi(k).freq,edi(k).rhoyx,'g');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).rhoyx(sel),edi(k).drhoyx(sel),'gs');
			hold on;
			kk = kk + 1; legstr{kk} = 'rhoyx';
		end
		if dispcomp{5}
			sel = edi(k).zxxs > 0;
			loglog(edi(k).freq,edi(k).zxxr,'y');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zxxr(sel),edi(k).dzxx(sel),'yv');
			hold on;
			kk = kk + 1; legstr{kk} = 'zxxr';
		end
		if dispcomp{6}
			sel = edi(k).zxxs > 0;
			loglog(edi(k).freq,edi(k).zxxi,'k');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zxxi(sel),edi(k).dzxx(sel),'k^');
			hold on;
			kk = kk + 1; legstr{kk} = 'zxxi';
		end
		if dispcomp{7}
			sel = edi(k).zxys > 0;
			loglog(edi(k).freq,edi(k).zxyr,'m');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zxyr(sel),edi(k).dzxy(sel),'m>');
			hold on;
			kk = kk + 1; legstr{kk} = 'zxyr';
		end
		if dispcomp{8}
			sel = edi(k).zxys > 0;
			loglog(edi(k).freq,edi(k).zxyi,'c');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zxyi(sel),edi(k).dzxy(sel),'c<');
			hold on;
			kk = kk + 1; legstr{kk} = 'zxyi';
		end
		if dispcomp{9}
			sel = edi(k).zyxs > 0;
			loglog(edi(k).freq,edi(k).zyxr,'y');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zyxr(sel),edi(k).dzyx(sel),'y>');
			hold on;
			kk = kk + 1; legstr{kk} = 'zyxr';
		end
		if dispcomp{10}
			sel = edi(k).zyxs > 0;
			loglog(edi(k).freq,edi(k).zyxi,'k');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zyxi(sel),edi(k).dzyx(sel),'k<');
			hold on;
			kk = kk + 1; legstr{kk} = 'zyxi';
		end
		if dispcomp{11}
			sel = edi(k).zyys > 0;
			loglog(edi(k).freq,edi(k).zyyr,'m');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zyyr(sel),edi(k).dzyy(sel),'mv');
			hold on;
			kk = kk + 1; legstr{kk} = 'zyyr';
		end
		if dispcomp{12}
			sel = edi(k).zyys > 0;
			loglog(edi(k).freq,edi(k).zyyi,'c');
			hold on;
			errorbare('vlogd',edi(k).freq(sel),edi(k).zyyi(sel),edi(k).dzyy(sel),'c^');
			hold on;
			kk = kk + 1; legstr{kk} = 'zyyi';
		end
		grid on
		set(gca,'xdir','reverse');
		legend(legstr);
		hold off;

%----------------------------------------------------------------------%
	function dispedi_2(handles,k)
		global dispcomp edi
		axes(handles.tAxesDisp);
		% generate a bi-vertical-axis first.
		[AX,H1,H2] = plotyy(edi(k).freq,edi(k).rhoxy,edi(k).freq,edi(k).phsxy,@loglog,@semilogx);
		xlabel('Frequency(Hz)');
		title(edi(k).filename);
		set(H1,'Visible','off');
		set(H2,'Visible','off');
		set(AX,'xdir','reverse');
		set(AX(1),'YLim',[1e-4 1e4]);
		set(AX(1),'YTick',[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4]);
		set(AX(2),'YLim',[-180 180]);
		set(AX(2),'YTick',[-180 -90 0 90 180]);
		set(AX,'YMinorTick','off');
% 		set(AX,'XGrid','on');
% 		set(AX(2),'YGrid','on');
		set(AX,'nextplot','add');
		drawnow
		% draw components.
		% AX(1) for left vertical-axis.
		% AX(2) for right vertical-axis.
		edi(k) = edi(k).asel2csel;
		if dispcomp{1}
			sel = edi(k).rhoxys > 0;
			axes(AX(1));
			loglog(edi(k).freq,edi(k).rhoxy,'b');
			errorbare('vlogd',edi(k).freq(sel),edi(k).rhoxy(sel),edi(k).drhoxy(sel),'bo');
		end
		if dispcomp{2}
			sel = edi(k).rhoyxs > 0;
			axes(AX(1));
			loglog(edi(k).freq,edi(k).rhoyx,'g');
			errorbare('vlogd',edi(k).freq(sel),edi(k).rhoyx(sel),edi(k).drhoyx(sel),'gs');
		end
		if dispcomp{3}
			sel = edi(k).phsxys > 0;
			axes(AX(2));
			semilogx(edi(k).freq,edi(k).phsxy,'g');
			errorbare('vlogx',edi(k).freq(sel),edi(k).phsxy(sel),edi(k).dphsxy(sel),'b+');
		end
		if dispcomp{4}
			sel = edi(k).phsyxs > 0;
			axes(AX(2));
			semilogx(edi(k).freq,edi(k).phsyx,'g');
			errorbare('vlogx',edi(k).freq(sel),edi(k).phsyx(sel),edi(k).dphsyx(sel),'gx');
		end
		%if dispcomp{5}
			%sel = edi(k).zxxs > 0;
			%loglog(edi(k).freq,edi(k).zxxr,'y');
			%hold on;
			%errorbare('vlogd',edi(k).freq(sel),edi(k).zxxr(sel),edi(k).dzxx(sel),'yv');
			%hold on;
		%end
		%if dispcomp{6}
			%sel = edi(k).zxxs > 0;
			%loglog(edi(k).freq,edi(k).zxxi,'k');
			%hold on;
			%errorbare('vlogd',edi(k).freq(sel),edi(k).zxxi(sel),edi(k).dzxx(sel),'k^');
			%hold on;
		%end
		%if dispcomp{7}
			%sel = edi(k).zxys > 0;
			%loglog(edi(k).freq,edi(k).zxyr,'m');
			%hold on;
			%errorbare('vlogd',edi(k).freq(sel),edi(k).zxyr(sel),edi(k).dzxy(sel),'m>');
			%hold on;
		%end
		%if dispcomp{8}
			%sel = edi(k).zxys > 0;
			%loglog(edi(k).freq,edi(k).zxyi,'c');
			%hold on;
			%errorbare('vlogd',edi(k).freq(sel),edi(k).zxyi(sel),edi(k).dzxy(sel),'c<');
			%hold on;
		%end
		%if dispcomp{9}
			%sel = edi(k).zyxs > 0;
			%loglog(edi(k).freq,edi(k).zyxr,'y');
			%hold on;
			%errorbare('vlogd',edi(k).freq(sel),edi(k).zyxr(sel),edi(k).dzyx(sel),'y>');
			%hold on;
		%end
		%if dispcomp{10}
			%sel = edi(k).zyxs > 0;
			%loglog(edi(k).freq,edi(k).zyxi,'k');
			%hold on;
			%errorbare('vlogd',edi(k).freq(sel),edi(k).zyxi(sel),edi(k).dzyx(sel),'k<');
			%hold on;
		%end
		%if dispcomp{11}
			%sel = edi(k).zyys > 0;
			%loglog(edi(k).freq,edi(k).zyyr,'m');
			%hold on;
			%errorbare('vlogd',edi(k).freq(sel),edi(k).zyyr(sel),edi(k).dzyy(sel),'mv');
			%hold on;
		%end
		%if dispcomp{12}
			%sel = edi(k).zyys > 0;
			%loglog(edi(k).freq,edi(k).zyyi,'c');
			%hold on;
			%errorbare('vlogd',edi(k).freq(sel),edi(k).zyyi(sel),edi(k).dzyy(sel),'c^');
			%hold on;
		%end
% 		grid on
% 		set(gca,'xdir','reverse');
% 		legend(legstr);
		hold off;

%----------------------------------------------------------------------%
	function edi = loadedi(edilist)
		[nl,nc] = size(edilist);
		for k = 1:nc
			edi(k) = classEDI(edilist{k});
			edi(k) = edi(k).complete;
			% convert Z to Z*sqrt(T);
			sT = sqrt(1./edi(k).freq);
			edi(k).zxxr = edi(k).zxxr .* sT;
			edi(k).zxxi = edi(k).zxxi .* sT;
			edi(k).zxyr = edi(k).zxyr .* sT;
			edi(k).zxyi = edi(k).zxyi .* sT;
			edi(k).zyxr = edi(k).zyxr .* sT;
			edi(k).zyxi = edi(k).zyxi .* sT;
			edi(k).zyyr = edi(k).zyyr .* sT;
			edi(k).zyyi = edi(k).zyyi .* sT;
		end

%----------------------------------------------------------------------%
	function rtv = getindex(x,freq)
		if (freq(1)>freq(end))
			if (x>freq(1)) rtv = 1; return;end
			if (x<freq(end)) rtv = length(freq); return;end
		else
			if (x<freq(1)) rtv = 1; return;end
			if (x>freq(end)) rtv = length(freq); return;end
		end % if out of freq range.
		for ki = 1:length(freq)-1
			if (x>min(freq(ki),freq(ki+1)) && x<max(freq(ki),freq(ki+1)))
				rtv(1) = ki;
				break;
			end % if
		end % ki

%----------------------------------------------------------------------%
	function bIsInAxes = IsInAxes(hObject, eventdata)
		minx = min(get(gca,'xlim'));
		maxx = max(get(gca,'xlim'));
		miny = min(get(gca,'ylim'));
		maxy = max(get(gca,'ylim'));
		Point=mean(get(gca,'currentpoint'));
		if((Point(1)-minx)*(maxx-Point(1))>0 && (Point(2)-miny)*(maxy-Point(2))>0 );
			bIsInAxes = 1;
		else
			bIsInAxes = 0;
		end

