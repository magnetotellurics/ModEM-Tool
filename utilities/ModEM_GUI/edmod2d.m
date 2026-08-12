function varargout = edmod2d(varargin)
% edmod2d - A GUI editor for 2D meshed model.
%
%  See also EDMOD3D
%  --------------------
%  Bo Yang, 2012.
%  Institute of Geophysics and Geomatics,
%  China University of Geosciences, Wuhan, China.
%  Comments, bug reports and questions, please send to:
%  yangbo.cug@163.com.
%  Copyright 2012-2016 Bo Yang, IGG, CUG.
%  $Revision: 1.0 $ $Date: 2012/12/25 $
%  Last change: 2013/11/07 15:03:52.

%  Revision log:
%  2012/12/25 : Version 1.0 released.
%  2013/11/05 : Version 1.1 released.

%
%  define and initialize the global arguments.
%
	global m
	global bButton pOrg pCur
	global bButton bDrawed bGrid bColorbar bAir bLog
	bButton   = 0;
	bDrawed   = 0;
	pOrg      = zeros(2,3);
	pCur      = pOrg;

	bGrid    = 1;
	bColorbar= 0;
	bAir     = 0;
	bLog     = 0;

	if (nargin < 1)
		m.MESHFILE = '';
		m.xn  = 0;
		m.zn  = 0;
		m.pp  = 0;
	end % nargin.
%
%  define the GUI.
%
	hMainGUIFig   = figure(...       % the main GUI figure 
		                'Toolbar','none', ...
		                'MenuBar','none', ...
		                'HandleVisibility','callback', ...
		                'Name', 'EDMOD2D 1.0', ...
		                'NumberTitle','off', ...
		                'units','normalized',...
		                'tag','tMainFig',...
		                'position',[0.01,0.06,0.98,0.86],...
		                'Color', get(0, 'defaultuicontrolbackgroundcolor'));
%--------------------------------menu items--------------------------------
	mFile       =   uimenu(...   % menu File.
		                'Parent',hMainGUIFig,...
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
		                'Callback',@mSave_Callback);
		mSaveAs =   uimenu(...  % menu save model as
		                'Parent', mFile,...
		                'Label','Save Model As ...',...
		                'Callback',@mSaveAs_Callback);
		mExport =   uimenu(...  % menu export model
		                'Parent', mFile,...
		                'Label','Export ...',...
		                'Enable','off',...
		                'Callback',@mExport_Callback);
		mLoadTXRX =   uimenu(...  % menu load txrx file.
		                'Parent', mFile,...
		                'Label','Load TXRX ...',...
		                'Enable','off',...
		                'Callback',@mLoadTXRX_Callback);
		mQuit =   uimenu(...  % menu exit
		                'Parent', mFile,...
		                'Label','Quit',...
		                'Separator','on',...
		                'Accelerator','Q',...
		                'Callback','close all');
	mEdit       =   uimenu(...   % menu Edit.
		                'Parent',hMainGUIFig,...
		                'Label','Edit');
		mXMesh =   uimenu(...  % menu edit x-mesh
		                'Parent', mEdit,...
		                'Label','X-Mesh ...',...
		                'Enable','off',...
		                'Callback',@mEditXMesh_Callback);
		mZMesh =   uimenu(...  % menu edit z-mesh
		                'Parent', mEdit,...
		                'Label','Z-Mesh ...',...
		                'Enable','off',...
		                'Callback',@mEditZMesh_Callback);
		mSetBG =   uimenu(...  % menu set background
		                'Parent', mEdit,...
		                'Label','Set Backgound ...',...
		                'Enable','off',...
		                'Callback',@mSetBG_Callback);
		mAddAir =   uimenu(...  % menu add air layers
		                'Parent', mEdit,...
		                'Label','Add air layer ...',...
		                'Enable','off',...
		                'Callback',@mAddAir_Callback);
	mView       =   uimenu(...   % menu View.
		                'Parent',hMainGUIFig,...
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
		mColorBar =   uimenu(...  % menu color bar
		                'Parent', mView,...
		                'Label','Colorbar',...
		                'Checked','off',...
		                'Callback',@mColorBar_Callback);
		mAirLayers =   uimenu(...  % menu air layers
		                'Parent', mView,...
		                'Label','Air layers',...
		                'Checked','off',...
		                'Enable','off',...
		                'Callback',@mAirLayers_Callback);
		mTXRX   =   uimenu(...  % menu view txrx
		                'Parent', mView,...
		                'Label','TXRX',...
		                'Enable','off',...
		                'Checked','off',...
		                'Callback',@mTXRX_Callback);
	mHelp       =   uimenu(...   % menu Help.
		                'Parent',hMainGUIFig,...
		                'Label','Help');
		mLoadHelp=   uimenu(...  % menu load the help doc.
		                'Parent', mHelp,...
		                'Label','Help',...
		                'Enable','off',...
		                'Callback',@mLoadHelp_Callback);
		mAbout   =   uimenu(...  % menu load about dlg.
		                'Parent', mHelp,...
		                'Label','About',...
		                'Callback',@mAbout_Callback);
%---------------------------define the main GUIs------------------------
	hAxesModel   =    axes(...         % the axes for edit model
		                'Parent', hMainGUIFig, ...
		                'Units', 'normalized', ...
		                'DrawMode','fast',...
		                'Visible','off',...
		                'HandleVisibility','callback', ...
		                'tag','tAxesModel',...
		                'Position',[0.05 0.05 0.80 0.9]);
	hTextLeftBrush  =  uicontrol(...    % text "Left Brush"
		                'Parent', hMainGUIFig, ...
		                'Style','Text',...
		                'Units','normalized',...
		                'Position',[0.91 0.92 0.08 0.06],...
		                'FontSize',11,...
		                'String','Left Brush');
	hEditLeftBrush  =  uicontrol(...    % edit "LeftBrush"
		                'Parent', hMainGUIFig, ...
		                'Style','edit',...
		                'Units','normalized',...
		                'BackgroundColor',([1 1 1]),...
		                'FontSize',10,...
		                'String','0.001',...
		                'tag','tEditLeftBrush',...
		                'Position',[0.92 0.90 0.07 0.04]);
	hTextRightBrush =  uicontrol(...    % text "Right Brush"
		                'Parent', hMainGUIFig, ...
		                'Style','Text',...
		                'Units','normalized',...
		                'Position',[0.91 0.83 0.08 0.06],...
		                'FontSize',11,...
		                'String','Right Brush');
	hEditRightBrush  =  uicontrol(...    % edit "RightBrush"
		                'Parent', hMainGUIFig, ...
		                'Style','edit',...
		                'Units','normalized',...
		                'BackgroundColor',([1 1 1]),...
		                'FontSize',10,...
		                'String','0.01',...
		                'tag','tEditRightBrush',...
		                'Position',[0.92 0.81 0.07 0.04]);
%
%  set all msg response functions.
%
	structData = struct('bButton',[0],'pOrg',[],'pCur',[],'pp',zeros(1,1));
	set(hMainGUIFig,'WindowButtonMotionFcn',{@ButtonMove_Callback}); 
	set(hMainGUIFig,'WindowButtonDownFcn',{@ButtonDown_Callback}); 
	set(hMainGUIFig,'WindowButtonUpFcn',{@ButtonUp_Callback});
	set(hMainGUIFig,'CloseRequestFcn',{@CloseRequestFcn_Callback});
%
%  save the gui data.
%
	handles = guihandles(hMainGUIFig);
	guidata(hMainGUIFig,handles);
%
%  deal with the input arguments.
%
	if (length(varargin) > 0)
		arg = varargin{1};
		if (ischar(arg))
			if (exist(arg,'file'))
				m.MESHFILE = arg;
				[m.nx,m.nz,m.nza,m.xn,m.zn,m.pp] = ReadMesh2D(m.MESHFILE);
			else
				str = strcat('ERROR: File',arg,'does not exist!');
				error(str);
			end % file exist.
		elseif (isstruct(arg))
			m.xn  = arg.xn;
			m.zn  = arg.zn;
			m.pp  = arg.pp;
		end % meshfile or model struct.
		drawmod(handles);
	end % called with arguments.
%
%  deal with the gui output.
%
	uiwait(hMainGUIFig);
	if (nargout > 0)
		varargout(1) = {m.pp};
	end % nargout.
%
%  End of main GUI function.
%
%----------------------------------------------------------------------%
%                             menu callbacks                           %
%----------------------------------------------------------------------%
	function mAbout_Callback(hObject, eventdata)
		msg{1} = 'simuGUI';
		msg{2} = 'A GUI tool for editing 2D mesh model ';
		msg{3} = 'Version 0.1';
		msg{4} = 'Dec 25, 2012';
		msg{5} = 'Author: Bo Yang';
		msg{6} = 'E-mail: yangbo.cug@163.com';
		msg{7} = 'Institute of Geophysics & Geomatics';
		msg{8} = 'China Univ. of Geosciences';
		msg{9} = 'Copyright: 2012-2016';
		h = msgbox(msg,'About EDMOD2D');

%----------------------------------------------------------------------%
	function mOpen_Callback(hObject, eventdata)
		global m
		%handles = guihandles(gcf);
		handles = guidata(gcbo);
		[FILE,PATH,FINDEX] = uigetfile('*.mesh', 'Open mesh file.');
		m.MESHFILE = strcat(PATH,FILE);
		[m.nx,m.nz,m.nza,m.xn,m.zn,m.pp] = ReadMesh2D(m.MESHFILE);
		drawmod(handles);

%----------------------------------------------------------------------%
	function mSaveAs_Callback(hObject, eventdata)
		global m
		[FILE,PATH,FINDEX] = uiputfile('*.mesh', 'Save as.');
		m.MESHFILE = strcat(PATH,FILE);
		WriteMesh2D(m.MESHFILE,m.nx,m.nz,m.nza,m.xn,m.zn,m.pp);

%----------------------------------------------------------------------%
	function mSave_Callback(hObject, eventdata)
		global m
		WriteMesh2D(m.MESHFILE,m.nx,m.nz,m.nza,m.xn,m.zn,m.pp);

%----------------------------------------------------------------------%
	function mGrid_Callback(hObject, eventdata)
		global bGrid
		%handles = guihandles(gcf);
		handles = guidata(gcbo);
		if strcmp(get(gcbo,'Checked'),'on')
			set(gcbo,'Checked','off');
			bGrid = 0;
			drawmod(handles);
		else
			set(gcbo,'Checked','on');
			bGrid = 1;
			drawmod(handles);
		end

%----------------------------------------------------------------------%
	function mLog_Callback(hObject, eventdata)
		global bLog
		%handles = guihandles(gcf);
		handles = guidata(gcbo);
		if strcmp(get(gcbo,'Checked'),'on')
			set(gcbo,'Checked','off');
			bLog = 0;
			drawmod(handles);
		else
			set(gcbo,'Checked','on');
			bLog = 1;
			drawmod(handles);
		end

%----------------------------------------------------------------------%
	function mColorBar_Callback(hObject, eventdata)
		global bColorbar
		%handles = guihandles(gcf);
		handles = guidata(gcbo);
		if strcmp(get(gcbo,'Checked'),'on')
			set(gcbo,'Checked','off');
			bColorbar = 0;
			drawmod(handles);
		else
			set(gcbo,'Checked','on');
			bColorbar = 1;
			drawmod(handles);
		end


%----------------------------------------------------------------------%
%                              msg callbacks                           %
%----------------------------------------------------------------------%
	function varargout = CloseRequestFcn_Callback(h,eventdata,handles)
		uiresume(h);
		delete(h);

%----------------------------------------------------------------------%
	function ButtonDown_Callback(h, eventdata, handles, varargin)
		global pOrg pCur bButton
		global bDrawed
		if (bDrawed < 1)
			return
		end
		bIsInAxes = IsInAxes;
		if bIsInAxes ~= 1
			return
		else
			pOrg = get(gca,'CurrentPoint');
			pCur = pOrg;
			ms = get(get(gca,'parent'),'SelectionType');
			bButton = 0;
			if strcmp(ms,'normal')
				bButton = -1;  %LBtn down
			else
				bButton = 1;   %RBtn down
			end
		end

%----------------------------------------------------------------------%
	function ButtonMove_Callback(h, eventdata, handles, varargin)
		global pOrg pCur bButton
		global bDrawed
		if (bDrawed < 1)
			return
		end
		bIsInAxes = IsInAxes;
		if bIsInAxes ~= 1
			return
		else
%         ShowPointData;
			%cp  = get(gca,'CurrentPoint');
			%cpx = cp(1,1); cpy = cp(1,2);
			%str = ['(',num2str(cpx,'%0.3g'),',',num2str(cpy,'%0.3g'),')'];
			%text(cpx,cpy,str,'VerticalAlignment','bottom');drawnow;
			if bButton ~= 0     %Btn down
				DynV_Sel_Box = line(...
				              'Color','r',...
				              'EraseMode','Xor',...
				              'XData',[pOrg(1,1) pCur(1,1) pCur(1,1) pOrg(1,1) pOrg(1,1)],...
				              'YData',[pOrg(1,2) pOrg(1,2) pCur(1,2) pCur(1,2) pOrg(1,2)]);
				drawnow;
				pCur = get(gca,'CurrentPoint');
				DynV_Sel_Box = line(...
				              'Color','r',...
				              'EraseMode','Xor',...
				              'XData',[pOrg(1,1) pCur(1,1) pCur(1,1) pOrg(1,1) pOrg(1,1)],...
				              'YData',[pOrg(1,2) pOrg(1,2) pCur(1,2) pCur(1,2) pOrg(1,2)]);
				drawnow;
			end
		end

%----------------------------------------------------------------------%
	function ButtonUp_Callback(h, eventdata, handles, varargin)    
		global pOrg pCur bButton
		global m
		global bDrawed
		if (bDrawed < 1)
			return
		end
		bIsInAxes = IsInAxes;
		if bIsInAxes ~= 1
			return
		end
%     handles = guidata(gcf);
		lhandles = guihandles(gcf);
		handles = guidata(gcf);
		lval = eval(get(lhandles.tEditLeftBrush,'String'));
		rval = eval(get(lhandles.tEditRightBrush,'String'));
		XMesh = m.xn;
		ZMesh = m.zn;
		ticx = XMesh(2) - XMesh(1);
		ticy = ZMesh(2) - ZMesh(1);
		t2 = ZMesh;t1 = XMesh;
    %set range
		XS = pOrg(1,1); YS = pOrg(1,2);
		XE = pCur(1,1); YE = pCur(1,2);
		if XS > XE
			Temp = XS; XS = XE; XE = Temp;
		end
		if YS > YE
			Temp = YS; YS = YE; YE = Temp;
		end
    %set value in the range
		for jj = 1:length(t2)
			if t2(jj) >= YS-ticy & t2(jj) <= YE
				for ii = 1:length(t1)
					if t1(ii) >= XS-ticx & t1(ii) <= XE
						if bButton == -1  %LBtn down, set up.
							m.pp(jj,ii) = lval;
						elseif bButton == 1 %RBtn down, delete.
							m.pp(jj,ii) = rval;
						else
						%no change
						end
					end
				end
			end
		end
    %update
		drawmod(handles);
		bButton = 0;
		guidata(gcf,handles);

%----------------------------------------------------------------------%
%                           general functions                          %
%----------------------------------------------------------------------%
	function [nx,nz,nza,xn,zn,pp] = ReadMesh2D(MESHFILE)
		fid = fopen(MESHFILE);
		A   = fscanf(fid,'%d %d %d',3);
		nx  = A(1);
		nz  = A(2);
		nza = A(3);
		tmp = fgetl(fid); % skip the comments.
		xn  = fscanf(fid,'%f',[nx 1]);
		zn  = fscanf(fid,'%f',[nz 1]);
		pp  = fscanf(fid,'%f',[1 (nx-1)*(nz-1)]);
		pp  = reshape(pp,nx-1,nz-1);
		pp  = pp';
		fclose(fid);

%----------------------------------------------------------------------%
	function [] = WriteMesh2D(MESHFILE,nx,nz,nza,xn,zn,pp)
		fid = fopen(MESHFILE,'w');
		fprintf(fid,'%d %d %d %s\n',nx,nz,nza,'!ny nz nza');
		fprintf(fid,'%f  ',xn);fprintf(fid,'\n');
		fprintf(fid,'%f  ',zn);fprintf(fid,'\n');
		for k = 1:nz-1
			for j = 1:nx-1
				fprintf(fid,'%e\n',pp(k,j));
			end % j
		end % k
		fclose(fid);

%----------------------------------------------------------------------%
	function drawmod(handles)
		global m
		global bGrid bColorbar bAir bLog
		global bDrawed
		temp = size(m.pp);
		if ((length(m.xn)-temp(2)) == 1)
			sig = wextend('2d','sym',m.pp,1,'d');
		else
			sig = m.pp;
		end % size of sig.
		if (bLog)
			pcolor(handles.tAxesModel,m.xn,m.zn,log10(sig));
		else
			pcolor(handles.tAxesModel,m.xn,m.zn,sig);
		end % bLog.
		set(handles.tAxesModel,'ydir','reverse');
		if (bGrid > 0)
			shading faceted;
		else
			shading flat;
		end % bGrid
		if (bColorbar > 0)
			colorbar;
		else
			colorbar('off');
		end % bColorbar
		bDrawed = 1;

%----------------------------------------------------------------------%
	function bIsInAxes = IsInAxes(hObject, eventdata)
		minx = min(get(gca,'xlim'));
		maxx = max(get(gca,'xlim'));
		miny = min(get(gca,'ylim'));
		maxy = max(get(gca,'ylim'));
		Point=mean(get(gca,'currentpoint'));
		if((Point(1)-minx)*(maxx-Point(1))>0 && (Point(2)-miny)*(maxy-Point(2))>0 )
			bIsInAxes = 1;
		else
			bIsInAxes = 0;
		end

