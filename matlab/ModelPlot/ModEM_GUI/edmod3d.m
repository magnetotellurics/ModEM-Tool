function varargout = edmod3d(varargin)
% edmod3d - A GUI tool for editing and displaying 3D model.
%
%  --------------------
%  Bo Yang, 2013.
%  Institute of Geophysics and Geomatics,
%  China University of Geosciences, Wuhan, China.
%  College of Earth, Ocean, and Atomspheric Sciences,
%  Oregon State University, Corvallis, USA.
%  Comments, bug reports and questions, please send to:
%  yangbo@cug.edu.cn or yang@coas.oregonstate.edu.
%  Copyright 2013-2017 Bo Yang, IGG, CUG.
%  $Revision: 1.0 $ $Date: 2013/10/29 $
%  Last change: 2013/11/08 10:26:51.

%  Revision log:
%  2013/10/29 : Created.

%
%  define and initialize the global arguments.
%
	global model

%
%  define the GUI.
%
	hEdmod3dMainFig = figure(... % the main GUI figure of edmod3d.
	                  'ToolBar','none',...
	                  'MenuBar','none',...
	                  'HandleVisibility','callback',...
	                  'Name','edmod3d',...
	                  'NumberTitle','off',...
	                  'Units','Normalized',...
	                  'tag','tEdmod3dMainFig',...
	                  'Position',[0.01,0.45,0.3,0.35],...
	                  'Color', get(0, 'defaultuicontrolbackgroundcolor'));
%--------------------------------menu items-----------------------------
	mFile       =   uimenu(...   % menu File.
	                'Parent',hEdmod3dMainFig,...
	                'Label','File');
		mNew    =   uimenu(...  % menu new model
	                'Parent', mFile,...
	                'Label','New Model ...',...
	                'Enable','off',...
	                'Callback',@mNew_Callback);
		mOpen   =   uimenu(...  % menu open model
	                'Parent', mFile,...
	                'Label','Open Model ...',...
	                'Callback',@mOpen_Callback);
		mSave   =   uimenu(...  % menu save model
	                'Parent', mFile,...
	                'Label','Save Model',...
	                'Enable','off',...
	                'Callback',@mSave_Callback);
		mSaveAs =   uimenu(...  % menu save model as
	                'Parent', mFile,...
	                'Label','Save Model As ...',...
	                'Enable','off',...
	                'Callback',@mSaveAs_Callback);
		mExport =   uimenu(...  % menu export model
	                'Parent', mFile,...
	                'Label','Export ...',...
	                'Enable','off',...
	                'Callback',@mExport_Callback);
		mQuit =   uimenu(...  % menu exit
	                'Parent', mFile,...
	                'Label','Quit',...
	                'Separator','on',...
	                'Accelerator','Q',...
	                'Callback','close all');
	mView       =   uimenu(...   % menu View.
	                'Parent',hEdmod3dMainFig,...
	                'Label','View');
		mGrid =   uimenu(...  % menu grid on/off
	                'Parent', mView,...
	                'Label','Grid',...
	                'Checked','on',...
	                'Callback',@mGrid_Callback);
		mLog =   uimenu(...  % menu log on/off
	                'Parent', mView,...
	                'Label','Log',...
	                'Checked','off',...
	                'Callback',@mLog_Callback);
	mHelp       =   uimenu(...   % menu Help.
	                'Parent',hEdmod3dMainFig,...
	                'Label','Help');
		mLoadHelp=   uimenu(...  % menu load the help doc.
	                'Parent', mHelp,...
	                'Label','Help',...
	                'Enable','off',...
	                'Callback',@mLoadHelp_Callback);
		mAbout   =   uimenu(...  % menu load about dlg.
	                'Parent', mHelp,...
	                'Label','About',...
	                'Enable','off',...
	                'Callback',@mAbout_Callback);
%-------------------------define the main GUIs--------------------------
	% define the btn size.
	btnr = 2;           % number of btn rows.
	btnc = 4;           % number of btn column.
	btnw = 1.0/btnc;    % btn width.
	btnh = 0.3/btnr;    % btn height.
	tabh = 0.6;         % table height.
	sldh = 0.08;        % slider height.
	colname     = {'','# slice','Del','Visible','XYZ','# block','Position'};
	colfmt      = {'char','numeric','logical','logical',{'X' 'Y' 'Z'},'numeric','numeric'};
	coleditable = [false false true true true true true];
	colwidth    = {16,44,30,50,60,80,120};
	hTabSlice        = uitable(... % the table of properties of slices.
	                  'Parent',hEdmod3dMainFig,...
	                  'Units','Normalized',...
	                  'tag','tTabSlice',...
	                  'Position',[0.01,0.99-tabh,0.98,tabh],...
	                  'ColumnName',colname,...
	                  'ColumnFormat',colfmt,...
	                  'ColumnEditable',coleditable,...
	                  'ColumnWidth',colwidth,...
	                  'CellEditCallback',@tTabSlice_CellEditCallback,...
	                  'CellSelectionCallback',@tTabSlice_CellSelectionCallback,...
	                  'RowName',[]);
	hSliderSlice     = uicontrol(... % the slider to move the selected slice.
	                  'Style','slider',...
	                  'Parent',hEdmod3dMainFig,...
	                  'Units','Normalized',...
	                  'tag','tSliderSlice',...
	                  'Position',[0.01,0.99-tabh-0.01-sldh,0.98,sldh],...
	                  'Callback',@tSliderSlice_Callbak);
	hBtnAddSlice = uicontrol(... % button for adding slice.
	                  'Parent',hEdmod3dMainFig,...
	                  'Units','Normalized',...
	                  'HandleVisibility','callback',...
	                  'Position',[0.0,(btnr-1)*btnh,btnw,btnh],...
	                  'FontSize',10,...
	                  'tag','tBtnAddSlice',...
	                  'Enable','on',...
	                  'String','Add Slice',...
	                  'Callback',@tBtnAddSlice_Callback);
	hBtnDelSlice = uicontrol(... % button for deleting slice.
	                  'Parent',hEdmod3dMainFig,...
	                  'Units','Normalized',...
	                  'HandleVisibility','callback',...
	                  'Position',[0.0,(btnr-2)*btnh,btnw,btnh],...
	                  'FontSize',10,...
	                  'tag','tBtnDelSlice',...
	                  'Enable','on',...
	                  'String','Delete Slice',...
	                  'Callback',@tBtnDelSlice_Callback);
	hBtnLoadModel = uicontrol(... % button for loading model.
	                  'Parent',hEdmod3dMainFig,...
	                  'Units','Normalized',...
	                  'HandleVisibility','callback',...
	                  'Position',[1*btnw,(btnr-1)*btnh,btnw,btnh],...
	                  'FontSize',10,...
	                  'tag','tBtnLoadModel',...
	                  'Enable','on',...
	                  'String','Load Model',...
	                  'Callback',@tBtnLoadModel_Callback);
	hBtnSaveModel = uicontrol(... % button for saving model.
	                  'Parent',hEdmod3dMainFig,...
	                  'Units','Normalized',...
	                  'HandleVisibility','callback',...
	                  'Position',[1*btnw,(btnr-2)*btnh,btnw,btnh],...
	                  'FontSize',10,...
	                  'tag','tBtnSaveModel',...
	                  'Enable','on',...
	                  'String','Save Model',...
	                  'Callback',@tBtnSaveModel_Callback);
	hBtnClose = uicontrol(... % button for Closing.
	                  'Parent',hEdmod3dMainFig,...
	                  'Units','Normalized',...
	                  'HandleVisibility','callback',...
	                  'Position',[2*btnw,(btnr-2)*btnh,btnw,btnh],...
	                  'FontSize',10,...
	                  'tag','tBtnClose',...
	                  'Enable','on',...
	                  'String','Close',...
	                  'Callback',@tBtnClose_Callback);
	hBtnEditSlice= uicontrol(... % button for Editing selected slice.
	                  'Parent',hEdmod3dMainFig,...
	                  'Units','Normalized',...
	                  'HandleVisibility','callback',...
	                  'Position',[3*btnw,(btnr-2)*btnh,btnw,2*btnh],...
	                  'FontSize',10,...
	                  'tag','tBtnEditSlice',...
	                  'Enable','on',...
	                  'String','Edit Slice',...
	                  'Callback',@tBtnEditSlice_Callback);
%
%  set msg response functions.
%
	set(hEdmod3dMainFig,'CloseRequestFcn',{@CloseRequestFcn_Callback});
%
%  initialize a figure tag for display the 3D model slices.
%
	tDispFig = figure('Name','Display slices','visible','off');
	tAxSlice = axes;
	handles = guihandles(hEdmod3dMainFig);
	handles.tDispFig = tDispFig;
	handles.tAxSlice = tAxSlice;
	guidata(hEdmod3dMainFig,handles);
%
%  deal with input arguments.
%
	if (length(varargin) > 0)
		arg = varargin{1};
		if (ischar(arg))
			if (exist(arg,'file'))
				model.MESHFILE = arg;
			else
				str = strcat('ERROR: File',arg,'does not exist!');
				error(str);
			end % file exist.
		elseif (isstruct(arg))
			model.xn = arg.xn;
			model.yn = arg.yn;
			model.zn = arg.zn;
			model.v  = arg.v;
			if (isfield(arg,'slicetab'));
				model.slicetab = arg.slicetab;
			else
				nxc = round(length(model.xn)/2);
				nyc = round(length(model.yn)/2);
				nzc = round(length(model.zn)/2);
				xnc = model.xn(nxc);
				ync = model.yn(nyc);
				znc = model.zn(nzc);
				model.slicetab = {'',1,false,true,'X',nxc,xnc;...
				'',2,false,true,'Y',nyc,ync;'',3,false,true,'Z',nzc,znc};
			end % if difined the slices.
			model.iSelected = 0;
		end % meshfile or model struct inputed.
		set(hTabSlice,'Data',model.slicetab);
		Display(handles);
	end % called with model struct.
%
%  deal with gui output.
%
%	uiwait(hEdmod3dMainFig);
	if (nargout > 0)
		varargout(1) = {model.v};
	end % nargout.
%
%  End of main function.
%

%----------------------------------------------------------------------%
%                              btn callbacks                           %
%----------------------------------------------------------------------%
	function tBtnAddSlice_Callback(hObject,eventdata,handles)
		global model
		handles = guidata(gcbo);
		[ns,temp] = size(model.slicetab);
		for k = 1:ns
			temp(k) = model.slicetab{k,2};
		end % k
		temp = max(temp);
		model.slicetab(ns+1,:) = {'',temp+1,false,true,'X',1,1};
		set(handles.tTabSlice,'Data',model.slicetab);
		Display(handles);

%----------------------------------------------------------------------%
	function tBtnDelSlice_Callback(hObject,eventdata,handles)
		global model
		handles = guidata(gcbo);
		tabdata = get(handles.tTabSlice,'Data');
		[ns,temp] = size(tabdata);
		model.slicetab = {};
		k2 = 0;
		for k1 = 1:ns
			if ~tabdata{k1,3} 
				k2 = k2 + 1;
				model.slicetab(k2,:) = tabdata(k1,:);
			end % if not del.
		end % k1
		set(handles.tTabSlice,'Data',model.slicetab);
		Display(handles);

%----------------------------------------------------------------------%
	function tBtnLoadModel_Callback(hObject,eventdata,handles)
		msgbox('not implemented yet!');

%----------------------------------------------------------------------%
	function tBtnSaveModel_Callback(hObject,eventdata,handles)
		msgbox('not implemented yet!');

%----------------------------------------------------------------------%
	function tBtnClose_Callback(hObject,eventdata,handles)
		msgbox('not implemented yet!');

%----------------------------------------------------------------------%
	function tBtnEditSlice_Callback(hObject,eventdata,handles)
		global model
		handles = guidata(gcbo);
		kk = model.slicetab{model.iSelected,6};
		switch model.slicetab{model.iSelected,5}
		case 'X'
			m2d.xn  = model.yn;
			m2d.zn  = model.zn;
			pp      = squeeze(model.v(:,kk,:));
			m2d.pp  = pp';
			pp      = edmod2d(m2d);
			model.v(:,kk,:) = pp';
			set(handles.tTabSlice,'Data',model.slicetab);
			Display(handles);
		case 'Y'
			m2d.xn  = model.xn;
			m2d.zn  = model.zn;
			pp      = squeeze(model.v(kk,:,:));
			m2d.pp  = pp';
			pp      = edmod2d(m2d);
			model.v(kk,:,:) = pp';
			set(handles.tTabSlice,'Data',model.slicetab);
			Display(handles);
		case 'Z'
			m2d.xn  = model.xn;
			m2d.zn  = model.yn;
			m2d.pp  = squeeze(model.v(:,:,kk));
			model.v(:,:,kk) = edmod2d(m2d);
			set(handles.tTabSlice,'Data',model.slicetab);
			Display(handles);
		end % switch.

%----------------------------------------------------------------------%
	function tSliderSlice_Callbak(hObject,eventdata,handles)
		global model
		handles = guidata(gcbo);
		if (model.iSelected < 1)
			msgbox('Select one slice by clicking the # slice!');
		else
			model.slicetab{model.iSelected,7} = get(handles.tSliderSlice,'Value');
			set(handles.tTabSlice,'Data',model.slicetab);
			updatetab(7);
			Display(handles);
		end % not selected.

%----------------------------------------------------------------------%
%                          menu callbacks.                             %
%----------------------------------------------------------------------%
	function mOpen_Callback(hObject, eventdata)

%----------------------------------------------------------------------%
%                            msg callbacks                             %
%----------------------------------------------------------------------%
	function varargout = CloseRequestFcn_Callback(h,eventdata,handles)
		uiresume(h);
		delete(h);

%----------------------------------------------------------------------%
	function tTabSlice_CellEditCallback(hObject,eventdata,handles)
		global model
		icol = eventdata.Indices(2);
		updatetab(icol);
		handles = guidata(gcbo);
		model.slicetab = get(handles.tTabSlice,'Data');
		Display(handles);

%----------------------------------------------------------------------%
	function tTabSlice_CellSelectionCallback(hObject,eventdata,handles)
		global model
		ei = eventdata.Indices;
		[temp1,temp2] = size(ei);
		if (temp1*temp2 > eps)
			model.iSelected = ei(1,1);
		else
			return;
		end % if error ei.
		h = guidata(gcbo);
		model.slicetab(:,1) = {''};
		model.slicetab{model.iSelected,1} = '=>';
		set(h.tTabSlice,'Data',model.slicetab);
		switch model.slicetab{model.iSelected,5}
		case 'X'
			set(h.tSliderSlice,'Min',model.xn(1  ));
			set(h.tSliderSlice,'Max',model.xn(end));
			set(h.tSliderSlice,'Value',model.slicetab{model.iSelected,7});
		case 'Y'
			set(h.tSliderSlice,'Min',model.yn(1  ));
			set(h.tSliderSlice,'Max',model.yn(end));
			set(h.tSliderSlice,'Value',model.slicetab{model.iSelected,7});
		case 'Z'
			set(h.tSliderSlice,'Min',model.zn(1  ));
			set(h.tSliderSlice,'Max',model.zn(end));
			set(h.tSliderSlice,'Value',model.slicetab{model.iSelected,7});
		end % switch.

%----------------------------------------------------------------------%
%                          general functions                           %
%----------------------------------------------------------------------%
	function updatetab(icol)
		global model
		h = guidata(gcbo);
		%h = guidata(gcf);
		model.slicetab = get(h.tTabSlice,'Data');
		if icol == 6
			kk = model.slicetab{model.iSelected,6};
			switch model.slicetab{model.iSelected,5}
			case 'X'
				model.slicetab(model.iSelected,7) = {model.xn(kk)};
			case 'Y'
				model.slicetab(model.iSelected,7) = {model.yn(kk)};
			case 'Z'
				model.slicetab(model.iSelected,7) = {model.zn(kk)};
			end % switch
		elseif icol == 7
			kc = model.slicetab{model.iSelected,7};
			switch model.slicetab{model.iSelected,5}
			case 'X'
				kk = find(model.xn > kc);
				model.slicetab(model.iSelected,6) = {kk(1)};
			case 'Y'
				kk = find(model.yn > kc);
				model.slicetab(model.iSelected,6) = {kk(1)};
			case 'Z'
				kk = find(model.zn > kc);
				model.slicetab(model.iSelected,6) = {kk(1)};
			end % switch
		end % icol
		set(h.tTabSlice,'Data',model.slicetab);

%----------------------------------------------------------------------%
	function Display(handles)
		global model
		set(handles.tDispFig,'visible','on');
		tabdata = model.slicetab;
		[ns,temp] = size(tabdata);
		xs = []; ys = []; zs = [];
		for k = 1:ns
			if tabdata{k,4}
				switch tabdata{k,5}
				case 'X'
					xs = [xs tabdata{k,7}];
				case 'Y'
					ys = [ys tabdata{k,7}];
				case 'Z'
					zs = [zs tabdata{k,7}];
				end % switch.
			end % if visible.
		end % k
		slice(handles.tAxSlice,model.xn,model.yn,model.zn,model.v,xs,ys,zs);
		%colormap hsv;
		xlabel('X');ylabel('Y');zlabel('Z');
		set(gca,'Zdir','reverse');
		shading flat;

