function varargout = ModTracer(varargin)
% ModTracer - A GUI tool to trace specified features of 3D model.
%
%  Usage:
%    use X/Y/Z key to switch the slice direction.
%    use RIGHT/UP arrow key to go to the next slice.
%    use LEFT/DOWN arrow key to go to the previous slice.
%    use LEFT mouse button to select.
%    use RIGHT mouse button to unselect.
%    Use L/R mouse button down and move as brush to select or unselect.
%    Use <F2> to load an exised selection file.
%    Use <F3> to save current selection to a txt file.
%    Use <F4> to set the selection range. xrange=1 means one more block
%      will be selected at both the left and right sides of current selection.
%    Use <Shift+LEFT click> to active the displaying of current position.
%    Use standard figure toolbar to Zoom in and Zoom out.
%
%  NOTES: 
%  
%
%  See also PickModCells.
%  --------------------
%  Bo Yang, 2014.
%  Institute of Geophysics and Geomatics,
%  China University of Geosciences, Wuhan, China.
%  Comments, bug reports and questions, please send to:
%  yangbo.cug@163.com.
%  Copyright 2014-2018, Bo Yang, IGG, CUG.
%  $Revision: 1.0 $ $Date: 2014/02/26 $
%  Last changed: 2014/03/07 11:07:17.

%  Revision log:
%  2014/02/26 : Version 1.0 released.
%
%  initialize global arguments.
%
	global xn yn trace hrect xrange yrange bBtn
	global nx ny nz
	global kp direc
	xn = []; yn = []; trace = []; hrect = []; xrange = 0; yrange = 0; bBtn = 0;
	kp = 1;
	direc = 'Z';
%
%  get the mesh size info first.
%
	nx = length(ls('figs/X*.fig'));
	ny = length(ls('figs/Y*.fig'));
	nz = length(ls('figs/Z*.fig'));
%
%  deal with input arguments.
%
	if nargin == 1;
		trace = varargin{1};
	end % nargin.
%
%  define the GUI.
%
	hMainGUIFig  = figure(...       % the main GUI figure 'Toolbar','none', ...'MenuBar','none', ...
		                'HandleVisibility','callback', ...
		                'Name', 'ModTracer', ...
		                'NumberTitle','off', ...
		                'units','normalized',...
		                'tag','tMainFig',...
		                'position',[0.04,0.06,0.92,0.82],...
		                'Color', get(0, 'defaultuicontrolbackgroundcolor'));
%
%  set all msg response functions.
%
	set(hMainGUIFig,'WindowKeyPressFcn',{@WindowKeyPressFcn_Callback});
	set(hMainGUIFig,'WindowButtonDownFcn',{@ButtonDown_Callback});
	set(hMainGUIFig,'WindowButtonMotionFcn',{@ButtonMove_Callback}); 
	set(hMainGUIFig,'WindowButtonUpFcn',{@ButtonUp_Callback});
%
%  save the gui data.
%
	handles = guihandles(hMainGUIFig);
	guidata(hMainGUIFig,handles);
	
	FIGFILE = ['figs/','Z',sprintf('%3.3d',kp),'.fig'];
	hrt = loadfig(FIGFILE,hMainGUIFig);
	h  = get(hrt(2),'Children');
	axes(hrt(2));
	xn = get(h,'Xdata');
	yn = get(h,'Ydata');
%
% End of main function.
%
%----------------------------------------------------------------------%
%                              msg callbacks                           %
%----------------------------------------------------------------------%
	function WindowKeyPressFcn_Callback(h,evt)
	global kp direc trace xrange yrange nx ny nz
	switch direc
	case 'X'
		maxkp = nx;
	case 'Y'
		maxkp = ny;
	case 'Z'
		maxkp = nz;
	otherwise
	end % switch.
	switch evt.Key
	case {'rightarrow','uparrow'}
		kp = min(kp + 1,maxkp);
		dispslice(direc,kp)
	case {'leftarrow','downarrow'}
		kp = max(kp - 1,1);
		dispslice(direc,kp)
	case 'x'
		direc = 'X';
		kp = min(kp,nx);
		dispslice(direc,kp)
	case 'y'
		direc = 'Y';
		kp = min(kp,ny);
		dispslice(direc,kp)
	case 'z'
		direc = 'Z';
		kp = min(kp,nz);
		dispslice(direc,kp)
	case 'f2' % load existed selection.
		% get the existed selection.
		[file,path] = uigetfile({'*.txt';'*.*'},'Choose a file');
		if (file ~= 0)
			trace = load([path file]);
		else
			return;
		end % check model file name
		% draw slice.
		handles = guidata(gcf);
		drawpicks(handles);
	case 'f3' % save selection to *.txt file.
		[file, path] = uiputfile('*.txt', 'Save selection as');
		if (file ~= 0)
			fid = fopen([path file],'w');
			for k = 1:length(trace)
				fprintf(fid,'%d %d %d\n',trace(k,1),trace(k,2),trace(k,3));
			end % k
			fclose(fid);
		else
			return;
		end % check model file name
	case 'f4' % set xrange and yrange.
		prompt = {'H range:','V range:'};
		dlg_title = 'Select range';
		num_lines = 1;
		def = {'1','1'};
		answer = inputdlg(prompt,dlg_title,num_lines,def);
		xrange = str2num(answer{1});
		yrange = str2num(answer{2});
		disp(['selection range: V=',sprintf('%d',xrange),' H=',sprintf('%d',yrange)]);
	otherwise
	end % switch key.

%----------------------------------------------------------------------%
	function ButtonDown_Callback(h, eventdata, handles, varargin)
		global trace MAP_PROJECTION xn yn xrange yrange bBtn
		bIsInAxes = IsInAxes;
		if bIsInAxes ~= 1
			return
		else
			CP = get(gca,'CurrentPoint');
			% if using m_map, m_proj has been set up
% 			if exist('MAP_PROJECTION','var')
% 				[x,y] = m_xy2ll(CP(1,1),CP(1,2));
% 			else
				x = CP(1,1); y = CP(1,2);
% 			end
% 			str = sprintf('Selected [%f,%f]',x,y);
% 			disp(str);
			updatetrace(x,y);
		end
		% draw slice.
		handles = guidata(gcf);
		drawpicks(handles);

%----------------------------------------------------------------------%
	function ButtonUp_Callback(h, eventdata, handles, varargin)
		global bBtn
		bBtn = 0;

%----------------------------------------------------------------------%
	function ButtonMove_Callback(h, eventdata, handles, varargin)
		global bBtn
		ms = get(get(gca,'parent'),'SelectionType');
		if strcmp(ms,'extend') %SHFIT + LBtn down
			CP = get(gca,'CurrentPoint');
			x = CP(1,1); y = CP(1,2);
			str = ['(',num2str(x,'%0.3g'),',',num2str(y,'%0.3g'),')'];
			if (isempty(findobj('tag','htext')))
				text(x,y,str,'tag','htext','VerticalAlignment','bottom','EraseMode','normal');
			else
				ht = findobj('tag','htext');
				set(ht,'position',[x,y],'string',str);
			end % if
		end
		if bBtn ~= 0 % L/R btn down.
			CP = get(gca,'CurrentPoint');
			x = CP(1,1); y = CP(1,2);
			updatetrace(x,y);
			handles = guidata(gcf);
			drawpicks(handles);
		end % L/R btn.

%----------------------------------------------------------------------%
%                           general functions                          %
%----------------------------------------------------------------------%
	function dispslice(direc,k)
		global hrect xn yn
		switch direc
		case 'X'
			FIGFILE = ['figs/','X',sprintf('%3.3d',k),'.fig'];
		case 'Y'
			FIGFILE = ['figs/','Y',sprintf('%3.3d',k),'.fig'];
		case 'Z'
			FIGFILE = ['figs/','Z',sprintf('%3.3d',k),'.fig'];
		otherwise
		end % switch.
		clf
		hrt=loadfig(FIGFILE,gcf);
		h  = get(hrt(2),'Children');
		xn = get(h,'Xdata');
		yn = get(h,'Ydata');
		% change the current axis to pcolor axis.
		axes(hrt(2));
		hrect = [];
		handles = guidata(gcf);
		drawpicks(handles);

%----------------------------------------------------------------------%
	function updatetrace(x,y)
		global trace xn yn xrange yrange bBtn nx ny nz kp direc
		rtv = getindex(x,y);
		% deal with the fast moving of mouse might return 0 index.
		if isnan(sum(rtv))
			return
		end % if.
		sel = rtv;
		for kx = sel(1,1)-xrange:sel(1,1)+xrange
			for ky = sel(1,2)-yrange:sel(1,2)+yrange
					if (kx >= 1 && kx < length(xn) && ...
						ky >= 1 && ky < length(yn))
						sel(end+1,1:2) = [kx ky];
					end % if
			end % ky
		end % kx
		sel = unique(sel,'rows');
		[nsel temp] = size(sel);
		% extend the selection indeces to 3D.
		switch direc
		case 'X'
			sel = [ones(nsel,1)*kp sel];
		case 'Y'
			sel = [sel(:,1) ones(nsel,1)*kp sel(:,2)];
		case 'Z'
			sel = [sel ones(nsel,1)*kp];
		otherwise
		end % switch.
		% add or remove current selection.
		ms = get(get(gca,'parent'),'SelectionType');
		if strcmp(ms,'normal') %LBtn down
			if ~isempty(trace)
				trace = union(trace,sel,'rows');
			else
				trace = sel;
			end % if
			bBtn = 1;
		elseif strcmp(ms,'alt') %RBtn down
			if ~isempty(trace)
				trace = setdiff(trace,sel,'rows');
			else
				return
			end % if
			bBtn = -1;
		end

%----------------------------------------------------------------------%
	function drawpicks(handles)
		global trace xn yn hrect MAP_PROJECTION direc kp
		if ~isempty(hrect)
			delete(hrect);
			hrect = [];
		end
		if isempty(trace)
			return
		end % if
		switch direc
		case 'X'
			ind = (trace(:,1)==kp);
			trace2d = trace(ind,:);
			trace2d = trace2d(:,2:3);
		case 'Y'
			ind = (trace(:,2)==kp);
			trace2d = trace(ind,:);
			trace2d = [trace2d(:,1) trace2d(:,3)];
		case 'Z'
			ind = (trace(:,3)==kp);
			trace2d = trace(ind,:);
			trace2d = trace2d(:,1:2);
		otherwise
		end % switch.
		if isempty(trace2d)
			return;
		end % if.
		ss = size(trace2d);
		kt = ss(1);
		%axes(gca);
		for kk = 1:kt
% 			if exist('MAP_PROJECTION','var')
% 				[xr,yr] = m_ll2xy(xn(trace(kk,1)),yn(trace(kk,2)));
% 				[xr2,yr2] = m_ll2xy(xn(trace(kk,1)+1),yn(trace(kk,2)+1));
% 				wr = xr2-xr;
% 				hr = yr2-yr;
% 			else
				xr = xn(trace2d(kk,1));
				yr = yn(trace2d(kk,2));
				wr = xn(trace2d(kk,1)+1) - xn(trace2d(kk,1));
				hr = yn(trace2d(kk,2)+1) - yn(trace2d(kk,2));
% 			end
			hrect(kk) = rectangle('Position',[xr,yr,wr,hr],'LineWidth',2,'EdgeColor','r');
		end % kk

%----------------------------------------------------------------------%
	function rtv = getindex(x,y)
		global xn yn
		rtv = [NaN NaN];
		for ki = 1:length(xn)-1
			if (x>xn(ki) && x<xn(ki+1))
				rtv(1) = ki;
				break;
			end % if
		end % ki
		for kj = 1:length(yn)-1
			if (y>yn(kj) && y<yn(kj+1))
				rtv(2) = kj;
				break;
			end % if
		end % kj

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

