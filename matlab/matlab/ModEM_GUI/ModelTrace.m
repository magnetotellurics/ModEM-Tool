function varargout = ModelTrace(varargin)
% ModelTrace - A GUI tool to extract a trace from a 3D model by clicking
%              the mouse.
%
%  --------------------
%  Bo Yang, 2013.
%  Institute of Geophysics and Geomatics,
%  China University of Geosciences, Wuhan, China.
%  Comments, bug reports and questions, please send to:
%  yangbo.cug@163.com.
%  Copyright 2012-2016 Bo Yang, IGG, CUG.
%  $Revision: 1.0 $ $Date: 2013/12/04 $
%  Last changed: 2013/12/04 14:55:20.

%  Revision log:
%  2013/12/04 : Created.

%
%  define and initialize the global variables.
%
	global cond
	global direc posx posy posz trace
	global bDrawed
	direc = 'XY';
	posx   = 1;
	posy   = 1;
	posz   = 1;
	trace = [];
	bDrawed = 0;
%
%  define the GUI.
%
	hMainGUIFig  = figure(...       % the main GUI figure 
		                'Toolbar','none', ...
		                'MenuBar','none', ...
		                'HandleVisibility','callback', ...
		                'Name', 'ModelTrace', ...
		                'NumberTitle','off', ...
		                'units','normalized',...
		                'tag','tMainFig',...
		                'position',[0.01,0.06,0.98,0.86],...
		                'Color', get(0, 'defaultuicontrolbackgroundcolor'));
%--------------------------------menu items--------------------------------
	mFile       =   uimenu(...   % menu File.
		                'Parent',hMainGUIFig,...
		                'Label','File');
		mOpen   =   uimenu(...  % menu open model
		                'Parent', mFile,...
		                'Label','Open Model ...',...
		                'Callback',@mOpen_Callback);
		mSave   =   uimenu(...  % menu save trace
		                'Parent', mFile,...
		                'Label','Save Trace ...',...
		                'Callback',@mSave_Callback);
		mQuit =   uimenu(...  % menu exit
		                'Parent', mFile,...
		                'Label','Quit',...
		                'Separator','on',...
		                'Accelerator','Q',...
		                'Callback','close all');
	mView       =   uimenu(...   % menu View.
		                'Parent',hMainGUIFig,...
		                'Label','View');
		mSlice   =   uimenu(...  % menu slice direction.
		                'Parent', mView,...
		                'Label','Slice Direction');
			mSliceXY =   uimenu(...  % menu slice direction xy.
		                'Parent', mSlice,...
		                'Label','X-Y',...
		                'Checked','on',...
		                'tag','tSliceXY',...
		                'Callback',@mSliceXY_Callback);
			mSliceXZ =   uimenu(...  % menu slice direction xz.
		                'Parent', mSlice,...
		                'Label','X-Z',...
		                'Checked','off',...
		                'tag','tSliceXZ',...
		                'Callback',@mSliceXZ_Callback);
			mSliceYZ =   uimenu(...  % menu slice direction yz.
		                'Parent', mSlice,...
		                'Label','Y-z',...
		                'Checked','off',...
		                'tag','tSliceYZ',...
		                'Callback',@mSliceYZ_Callback);
		mSlicePos =   uimenu(...  % menu slice position.
		                'Parent', mView,...
		                'Label','Slice Position ...',...
		                'Callback',@mSlicePos_Callback);
		mSliceNext =   uimenu(...  % menu next slice in current direction.
		                'Parent', mView,...
		                'Label','Next Slice',...
		                'Callback',@mSliceNext_Callback);
		mSlicePrev =   uimenu(...  % menu previous slice in current direction.
		                'Parent', mView,...
		                'Label','Previous Slice',...
		                'Callback',@mSlicePrev_Callback);
	hAxesModel   =    axes(...         % the axes for edit model
		                'Parent', hMainGUIFig, ...
		                'Units', 'normalized', ...
		                'DrawMode','fast',...
		                'Visible','off',...
		                'HandleVisibility','callback', ...
		                'tag','tAxesModel',...
		                'Position',[0.05 0.05 0.80 0.9]);
%
%  set all msg response functions.
%
	set(hMainGUIFig,'WindowButtonDownFcn',{@ButtonDown_Callback});
%
%  save the gui data.
%
	handles = guihandles(hMainGUIFig);
	guidata(hMainGUIFig,handles);
%
% End of main GUI function.
%
%----------------------------------------------------------------------%
%                             menu callbacks                           %
%----------------------------------------------------------------------%
	function mOpen_Callback(hObject, eventdata)
		global cond
		% get the model file name and path.
		[modelfile,modelpath] = uigetfile('*.rho','Choose a ModEM model file');
		% read the model file.
		if (modelfile ~= 0)
			cond = readCond_3D([modelpath modelfile],2);
		else
			return;
		end % check model file name
		nx = cond.grid.Nx;
		ny = cond.grid.Ny;
		nze= cond.grid.NzEarth;
		nza= cond.grid.NzAir;
		ori= cond.grid.origin;
		dx = cond.grid.dx;
		dy = cond.grid.dy;
		dz = cond.grid.dz;
		sig= cond.v;
		x0 = ori(1);
		y0 = ori(2);
		z0 = ori(2);
		% compute xn, yn, zn.
		xn(1) = x0;
		for k = 2:nx+1
			xn(k) = xn(k-1) + dx(k-1);
		end % k
		yn(1) = y0;
		for k = 2:ny+1
			yn(k) = yn(k-1) + dy(k-1);
		end % k
		zn(1) = z0;
		for k = 2:nze+1
			zn(k) = zn(k-1) + dz(k-1);
		end % k
		sige = zeros(nx+1,ny+1,nze+1);
		sige(1:nx,1:ny,1:nze) = sig;
		cond.xn = xn; cond.yn = yn; cond.zn = zn;
		cond.sig= sige;
		% draw slice.
		handles = guidata(gcf);
		drawslice(handles);

%----------------------------------------------------------------------%
	function mSave_Callback(hObject, eventdata)
		global trace
		global cond
		[FILE,PATH,FINDEX] = uiputfile('*.txt', 'Save.');
		fid = fopen([PATH FILE],'w');
		out = trace;
		ss = size(trace);
		kt = ss(1);
		for k = 1:kt
			out(k,4) = (cond.xn(trace(k,1)+1) + cond.xn(trace(k,1))) / 2.0;
			out(k,5) = (cond.yn(trace(k,2)+1) + cond.yn(trace(k,2))) / 2.0;
			out(k,6) = (cond.zn(trace(k,3)+1) + cond.zn(trace(k,3))) / 2.0;
			out(k,7) = cond.sig(trace(k,1),trace(k,2),trace(k,3));
			fprintf(fid,'%d %d %d %f %f %f %f\n',out(k,1),out(k,2),out(k,3),...
			        out(k,4),out(k,5),out(k,6),out(k,7));
		end % for
		fclose(fid);

%----------------------------------------------------------------------%
	function mSlicePos_Callback(hObject, eventdata)
		global posx posy posz
		prompt = {'# of X slice';'# of Y slice';'# of Z slice';};
		dlg_title = 'Input # of slice';
		num_lines = 1;
		def = {num2str(posx),num2str(posy),num2str(posz)};
		answer = inputdlg(prompt,dlg_title,num_lines,def);
		posx = str2num(answer{1});
		posy = str2num(answer{2});
		posz = str2num(answer{3});
		% draw slice.
		handles = guidata(gcf);
		drawslice(handles);

%----------------------------------------------------------------------%
	function mSliceXY_Callback(hObject, eventdata)
		global direc
		direc = 'XY';
		handles = guidata(gcf);
		set(handles.tSliceXY,'Checked','on')
		set(handles.tSliceXZ,'Checked','off')
		set(handles.tSliceYZ,'Checked','off')
		% draw slice.
		drawslice(handles);

%----------------------------------------------------------------------%
	function mSliceXZ_Callback(hObject, eventdata)
		global direc
		direc = 'XZ';
		handles = guidata(gcf);
		set(handles.tSliceXY,'Checked','off')
		set(handles.tSliceXZ,'Checked','on')
		set(handles.tSliceYZ,'Checked','off')
		% draw slice.
		drawslice(handles);

%----------------------------------------------------------------------%
	function mSliceYZ_Callback(hObject, eventdata)
		global direc
		direc = 'YZ';
		handles = guidata(gcf);
		set(handles.tSliceXY,'Checked','off')
		set(handles.tSliceXZ,'Checked','off')
		set(handles.tSliceYZ,'Checked','on')
		% draw slice.
		drawslice(handles);

%----------------------------------------------------------------------%
%                              msg callbacks                           %
%----------------------------------------------------------------------%
	function ButtonDown_Callback(h, eventdata, handles, varargin)
		global bDrawed trace
		if (bDrawed < 1)
			return
		end
		bIsInAxes = IsInAxes;
		if bIsInAxes ~= 1
			return
		else
			CP = get(gca,'CurrentPoint');
			ms = get(get(gca,'parent'),'SelectionType');
			bButton = 0;
			if strcmp(ms,'normal')
				bButton = -1;  %LBtn down
				rtv = getindex(CP);
				trace(end+1,1:3) = rtv;
				trace = unique(trace,'rows');
			else
				bButton = 1;   %RBtn down
				rtv = getindex(CP);
				trace(end+1,1:3) = rtv;
				[trace,km,kn] = unique(trace,'rows');
				kn = sort(kn);
				kn = kn(2:end) - kn(1:end-1);
				kdel = find(kn==0);
				trace(kdel,:) = [];
			end
		end
		% draw slice.
		handles = guidata(gcf);
		drawslice(handles);

%----------------------------------------------------------------------%
%                           general functions                          %
%----------------------------------------------------------------------%
	function drawslice(handles)
		global cond
		global direc posx posy posz trace
		global bDrawed
		ss = size(trace);
		kt = ss(1);
		switch direc
		case 'XY'
			sig2d  = squeeze(cond.sig(:,:,posz));
			pcolor(handles.tAxesModel,cond.yn,cond.xn,sig2d);
			hold on
			for kk = 1:kt
				if (trace(kk,3) == posz)
					xr = cond.yn(trace(kk,2));
					yr = cond.xn(trace(kk,1));
					wr = cond.yn(trace(kk,2)+1) - cond.yn(trace(kk,2));
					hr = cond.xn(trace(kk,1)+1) - cond.xn(trace(kk,1));
					rectangle('Position',[xr,yr,wr,hr],'LineWidth',2,'EdgeColor','r')
				end % if
			end % kk
			hold off
			bDrawed = 1;
		case 'XZ'
			sig2d  = squeeze(cond.sig(:,posy,:));
			pcolor(handles.tAxesModel,cond.xn,cond.zn,sig2d);
			set(gca,'YDir','reverse');
			hold on
			for kk = 1:kt
				if (trace(kk,2) == posy)
					xr = cond.xn(trace(kk,1));
					yr = cond.zn(trace(kk,3));
					wr = cond.xn(trace(kk,1)+1) - cond.xn(trace(kk,1));
					hr = cond.zn(trace(kk,3)+1) - cond.zn(trace(kk,3));
					rectangle('Position',[xr,yr,wr,hr],'LineWidth',2,'EdgeColor','r')
				end % if
			end % kk
			hold off
			bDrawed = 1;
		case 'YZ'
			sig2d  = squeeze(cond.sig(posx,:,:));
			pcolor(handles.tAxesModel,cond.yn,cond.zn,sig2d');
			set(gca,'YDir','reverse');
			hold on
			for kk = 1:kt
				if (trace(kk,1) == posx)
					xr = cond.yn(trace(kk,2));
					yr = cond.zn(trace(kk,3));
					wr = cond.yn(trace(kk,2)+1) - cond.yn(trace(kk,2));
					hr = cond.zn(trace(kk,3)+1) - cond.zn(trace(kk,3));
					rectangle('Position',[xr,yr,wr,hr],'LineWidth',2,'EdgeColor','r')
				end % if
			end % kk
			hold off
			bDrawed = 1;
		otherwise
		end % switch direc.

%----------------------------------------------------------------------%
	function rtv = getindex(cp)
		global cond
		global direc posx posy posz
		switch direc
		case 'XY'
			for ki = 1:length(cond.xn)-1
				if (cp(1,2)>cond.xn(ki) && cp(1,2)<cond.xn(ki+1))
					rtv(1) = ki;
					break;
				end % if
			end % ki
			for kj = 1:length(cond.yn)-1
				if (cp(1,1)>cond.yn(kj) && cp(1,1)<cond.yn(kj+1))
					rtv(2) = kj;
					break;
				end % if
			end % kj
			rtv(3) = posz;
		case 'XZ'
			for ki = 1:length(cond.xn)-1
				if (cp(1,1)>cond.xn(ki) && cp(1,1)<cond.xn(ki+1))
					rtv(1) = ki;
					break;
				end % if
			end % ki
			for kj = 1:length(cond.zn)-1
				if (cp(1,2)>cond.zn(kj) && cp(1,2)<cond.zn(kj+1))
					rtv(3) = kj;
					break;
				end % if
			end % kj
			rtv(2) = posy;
		case 'YZ'
			for ki = 1:length(cond.yn)-1
				if (cp(1,1)>cond.yn(ki) && cp(1,1)<cond.yn(ki+1))
					rtv(2) = ki;
					break;
				end % if
			end % ki
			for kj = 1:length(cond.zn)-1
				if (cp(1,2)>cond.zn(kj) && cp(1,2)<cond.zn(kj+1))
					rtv(3) = kj;
					break;
				end % if
			end % kj
			rtv(1) = posx;
		otherwise
		end % switch direc.

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

