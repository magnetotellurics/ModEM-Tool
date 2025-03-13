function varargout = mvemod3d(action,varargin)
% mvemod3d - A GUI tool for making, viewing and editing the 3D mesh model.
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
%  $Revision: 1.0 $ $Date: 2013/11/05 $
%  Last change: 2013/11/08 11:31:41.

%  Revision log:
%  2013/11/05 : Created.

	global model

	if (nargin < 1)
		action = 'initialize';
		model.file = '';
		model.path = '';
		model.ext  = '';
	end % nargin.

	if (length(varargin) > 0)
		arg = varargin{1};
		if isstruct(arg)
			model.file = arg.modelfile;
			model.path = arg.modelpath;
			model.ext  = arg.modelext;
		elseif ischar(arg)
			[model.path,model.file,model.ext] = fileparts(arg);
		else
			error('Unsupported arguments!');
		end % isstruct.
	end % called with model file name or a model struct.

	switch action
%----------------------------------------------------------------------%
	case 'initialize'
		% define the btn size.
		btnr = 9;      % number of btn rows.
		btnw = 0.98;   % btn width.
		btnh = 1/btnr; % btn height.
		hMvemod3dMainFig = figure(... % the main GUI figure of mvemod3d.
		                  'ToolBar','none',...
		                  'MenuBar','none',...
		                  'HandleVisibility','callback',...
		                  'Name','mvemod3d',...
		                  'NumberTitle','off',...
		                  'Units','Normalized',...
		                  'tag','tMainFig',...
		                  'Position',[0.01,0.25,0.16,0.7],...
		                  'Color', get(0, 'defaultuicontrolbackgroundcolor'));
		callbackCmd = 'mvemod3d(''open'')';
		hBtnOpenModel = uicontrol(... % button for opening the model file.
		                  'Parent',hMvemod3dMainFig,...
		                  'Units','Normalized',...
		                  'HandleVisibility','callback',...
		                  'Position',[0.01,1-btnh*1,btnw,btnh],...
		                  'FontSize',10,...
		                  'tag','tBtnOpenModel',...
		                  'Enable','on',...
		                  'String','Open Model file',...
		                  'Callback',callbackCmd);
		callbackCmd = 'mvemod3d(''new'');';
		hBtnNewModel = uicontrol(... % button for creating new the model file.
		                  'Parent',hMvemod3dMainFig,...
		                  'Units','Normalized',...
		                  'HandleVisibility','callback',...
		                  'Position',[0.01,1-btnh*2,btnw,btnh],...
		                  'FontSize',10,...
		                  'tag','tBtnNewModel',...
		                  'Enable','on',...
		                  'String','New Model file',...
		                  'Callback',callbackCmd);
		callbackCmd = 'mvemod3d(''loadmat'');';
		hBtnLoadMat = uicontrol(... % button for loading the model struct mat file.
		                  'Parent',hMvemod3dMainFig,...
		                  'Units','Normalized',...
		                  'HandleVisibility','callback',...
		                  'Position',[0.01,1-btnh*3,btnw,btnh],...
		                  'FontSize',10,...
		                  'tag','tBtnLoadMat',...
		                  'Enable','on',...
		                  'String','Load Model mat',...
		                  'Callback',callbackCmd);
		callbackCmd = 'mvemod3d(''savemat'');';
		hBtnSaveMat = uicontrol(... % button for saving the model struct mat file.
		                  'Parent',hMvemod3dMainFig,...
		                  'Units','Normalized',...
		                  'HandleVisibility','callback',...
		                  'Position',[0.01,1-btnh*4,btnw,btnh],...
		                  'FontSize',10,...
		                  'tag','tBtnSaveMat',...
		                  'Enable','on',...
		                  'String','Save Model mat',...
		                  'Callback',callbackCmd);
		callbackCmd = 'mvemod3d(''setrange'');';
		hBtnSetRange = uicontrol(... % button for setting the range of displaying.
		                  'Parent',hMvemod3dMainFig,...
		                  'Units','Normalized',...
		                  'HandleVisibility','callback',...
		                  'Position',[0.01,1-btnh*5,btnw,btnh],...
		                  'FontSize',10,...
		                  'tag','tBtnSetRange',...
		                  'Enable','on',...
		                  'String','Set Range',...
		                  'Callback',callbackCmd);
		callbackCmd = 'mvemod3d(''display'');';
		hBtnDispEdit = uicontrol(... % button for displaying and editing.
		                  'Parent',hMvemod3dMainFig,...
		                  'Units','Normalized',...
		                  'HandleVisibility','callback',...
		                  'Position',[0.01,1-btnh*6,btnw,btnh],...
		                  'FontSize',10,...
		                  'tag','tBtnDispEdit',...
		                  'Enable','on',...
		                  'String','Display & Edit',...
		                  'Callback',callbackCmd);
		callbackCmd = 'mvemod3d(''callmeshtools3d'');';
		hBtnCallMeshTools = uicontrol(... % button for calling meshtools3d.
		                  'Parent',hMvemod3dMainFig,...
		                  'Units','Normalized',...
		                  'HandleVisibility','callback',...
		                  'Position',[0.01,1-btnh*7,btnw,btnh],...
		                  'FontSize',10,...
		                  'tag','tBtnCallMeshTools',...
		                  'Enable','on',...
		                  'String','Call Meshtools3d',...
		                  'Callback',callbackCmd);
		callbackCmd = 'mvemod3d(''save'');';
		hBtnSave = uicontrol(... % button for saving the model file.
		                  'Parent',hMvemod3dMainFig,...
		                  'Units','Normalized',...
		                  'HandleVisibility','callback',...
		                  'Position',[0.01,1-btnh*8,btnw,btnh],...
		                  'FontSize',10,...
		                  'tag','tBtnSave',...
		                  'Enable','on',...
		                  'String','Save Model File',...
		                  'Callback',callbackCmd);
		callbackCmd = 'mvemod3d(''close'');';
		hBtnClose = uicontrol(... % button for closing.
		                  'Parent',hMvemod3dMainFig,...
		                  'Units','Normalized',...
		                  'HandleVisibility','callback',...
		                  'Position',[0.01,0,btnw,btnh],...
		                  'FontSize',10,...
		                  'tag','tBtnClose',...
		                  'Enable','on',...
		                  'String','Close',...
		                  'Callback',callbackCmd);
%----------------------------------------------------------------------%
	case 'open'
		% get the model file name and path.
		[model.file,model.path] = uigetfile({'*.rho';'*.*';},'Open a ModEM model file.');
		% read the model file.
		if (model.file ~= 0)
			cond = readCond_3D([model.path,model.file],2);
			model.nx  = cond.grid.Nx;
			model.ny  = cond.grid.Ny;
			model.nze = cond.grid.NzEarth;
			model.nza = cond.grid.NzAir;
			model.ori = cond.grid.origin;
			model.dx  = cond.grid.dx;
			model.dy  = cond.grid.dy;
			model.dz  = cond.grid.dz;
			model.v   = cond.v;
			x0 = model.ori(1); y0 = model.ori(2); z0 = model.ori(2);
			% compute xn, yn, zn.
			xn(1) = x0;
			for k = 2:model.nx+1
				xn(k) = xn(k-1) + model.dx(k-1);
			end % k
			yn(1) = y0;
			for k = 2:model.ny+1
				yn(k) = yn(k-1) + model.dy(k-1);
			end % k
			zn(1) = z0;
			for k = 2:model.nze+1
				zn(k) = zn(k-1) + model.dz(k-1);
			end % k
			model.xn = xn;
			model.yn = yn;
			model.zn = zn;
			model.disp.xn = model.xn;
			model.disp.yn = model.yn;
			model.disp.zn = model.zn;
			nx = length(model.xn); ny = length(model.yn); nz = length(model.zn);
			model.disp.v  = zeros(ny,nx,nz);
			for ky = 1:ny-1
				for kx = 1:nx-1
					for kz = 1:nz-1
						model.disp.v(ky,kx,kz) = model.v(kx,ky,kz);
					end % kz
				end % kx
			end % ky
			disp('Model file loaded.');
		else
			return;
		end % check model file name.
%----------------------------------------------------------------------%
	case 'new'
%----------------------------------------------------------------------%
		[model.file,model.path] = uiputfile({'*.rho';'*.*';},'Specify the name of the new ModEM model file.');
%----------------------------------------------------------------------%
	case 'setrange'
		msgbox('not implemented yet!');
%----------------------------------------------------------------------%
	case 'loadmat'
		[cfile,cpath] = uigetfile({'*.mat';},'Load MAT file.');
		model = load([cpath,cfile]);
		disp('model struct loaded.');
%----------------------------------------------------------------------%
	case 'savemat'
		[cfile,cpath] = uiputfile({'*.mat';},'Save to MAT file.');
		save([cpath,cfile],'-struct','model');
		disp('model struct saved.');
%----------------------------------------------------------------------%
	case 'display'
		out = edmod3d(model.disp);
%----------------------------------------------------------------------%
	case 'callmeshtools3d'
		[cpath,cfile,cext] = fileparts(model.file);
		% write out the mesh file for meshtools3d.
		cmesh = [model.path,cfile,'.mesh'];
		fid = fopen(cmesh,'w');
		fprintf(fid,'%d %d %d\n',model.nx,model.ny,model.nze);
		fprintf(fid,'%f %f %f\n',model.xn(1),model.yn(1),-model.zn(end));
		fprintf(fid,'%f ',model.dx);fprintf(fid,'\n');
		fprintf(fid,'%f ',model.dy);fprintf(fid,'\n');
		fprintf(fid,'%f ',model.dz);fprintf(fid,'\n');
		fclose(fid);
		% write out the model file for meshtools3d.
		cmodel = [model.path,cfile,'.model'];
		h = waitbar(0,'Writing model file...');
		fid = fopen(cmodel,'w');
		for kj = 1:model.ny
			waitbar(kj/model.ny);
			for ki = 1:model.nx
				for kk = 1:model.nze
					fprintf(fid,'%e\n',model.v(ki,kj,kk));
				end % kk
			end % ki
		end % kj
		fclose(fid);
		close(h);
		cmdstr = ['MeshTools3d',' "',cmesh,'" ','"',cmodel,'"'];
		system(cmdstr);
%----------------------------------------------------------------------%
	case 'save'
		msgbox('not implemented yet!');
%----------------------------------------------------------------------%
	case 'close'
		close all;
%----------------------------------------------------------------------%
	otherwise
		error('Unsupported action! You must specify an action.');
	end % switch action.

%
% End of the function.
%
