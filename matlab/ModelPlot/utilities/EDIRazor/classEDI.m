classdef classEDI
% classEDI - A class for IO of MT data in edi/plt format.
%
%  NOTES: 
%  
%
%  See also XXX.
%  --------------------
%  Bo Yang, 2014.
%  Institute of Geophysics and Geomatics,
%  China University of Geosciences, Wuhan, China.
%  Comments, bug reports and questions, please send to:
%  yangbo.cug@163.com.
%  Copyright 2014-2018 Bo Yang, IGG, CUG.
%  $Revision: 1.0 $ $Date: 2014/02/14 $
%  Last changed: 2014/10/27 18:56:14.

%  Revision log:
%  2014/02/14 : Version 1.0 released.
%  2014/03/26 : Version 1.1 released. Added 14 keywords for all
%               components selection.
%  2014/03/27 : Version 1.2 released. Separated the FSEL keywords as
%               a sole file from the original EDI file.
%  2014/04/28 : I) Modified the drho & dphs computation, treated the
%               varince as standard deviation.
%               II) shift phsyx +180 degree.
%  2015/01/22 : Fixed the bug in write function when using the default
%               filename.

%----------------------------------------------------------------------%
%  list the properties.                                                %
%----------------------------------------------------------------------%
	properties
		filename= '';
		dataid = '';
		lat    = 0.0;
		lon    = 0.0;
		elev   = 0.0;
		nfreq  = 0;
		asel   = [];
		fsel   = [];
		freq   = [];
		zxxr   = [];
		zxxi   = [];
		zxxs   = [];
		zxyr   = [];
		zxyi   = [];
		zxys   = [];
		zyxr   = [];
		zyxi   = [];
		zyxs   = [];
		zyyr   = [];
		zyyi   = [];
		zyys   = [];
		dzxx   = [];
		dzxy   = [];
		dzyx   = [];
		dzyy   = [];
		rhoxx  = [];
		rhoxxs = [];
		rhoxy  = [];
		rhoxys = [];
		rhoyx  = [];
		rhoyxs = [];
		rhoyy  = [];
		rhoyys = [];
		phsxx  = [];
		phsxxs = [];
		phsxy  = [];
		phsxys = [];
		phsyx  = [];
		phsyxs = [];
		phsyy  = [];
		phsyys = [];
		drhoxx = [];
		drhoxy = [];
		drhoyx = [];
		drhoyy = [];
		dphsxx = [];
		dphsxy = [];
		dphsyx = [];
		dphsyy = [];
		tzxr   = [];
		tzxi   = [];
		tzxs   = [];
		tzyr   = [];
		tzyi   = [];
		tzys   = [];
		dtzx   = [];
		dtzy   = [];
	end % properties.
	properties(Constant)
	end
	properties(SetAccess = public)
	end
	properties(SetAccess = protected)
	end
%----------------------------------------------------------------------%
%  define the methods.                                                 %
%----------------------------------------------------------------------%
	methods
%
%  constructor.
%
	function [obj] = classEDI(varargin)
		if nargin < 1
			return
		else
			if ischar(varargin{1})
				obj.filename = varargin{1};
				obj = obj.read;
			end
		end % empty constructor.
	end % constructor.
%
%  read edi.
%
	function obj = read(obj)
		if (isempty(obj.filename))
			error('No EDI file name information!');
		end % if empty filename.
		try
			[obj.dataid,obj.nfreq,obj.lat,obj.lon,obj.elev] = readedi(obj.filename,'>HEAD');
		catch err
			error(['No header information in EDI file: ',obj.filename,' OR wrong edi format!']);
		end % try read header info.
		% read all data components.
		try [obj.freq] = readedi(obj.filename,'>FREQ'); end
		try [obj.zxxr] = readedi(obj.filename,'>ZXXR'); end
		try [obj.zxxi] = readedi(obj.filename,'>ZXXI'); end
		try [obj.zxyr] = readedi(obj.filename,'>ZXYR'); end
		try [obj.zxyi] = readedi(obj.filename,'>ZXYI'); end
		try [obj.zyxr] = readedi(obj.filename,'>ZYXR'); end
		try [obj.zyxi] = readedi(obj.filename,'>ZYXI'); end
		try [obj.zyyr] = readedi(obj.filename,'>ZYYR'); end
		try [obj.zyyi] = readedi(obj.filename,'>ZYYI'); end
		try [obj.dzxx] = readedi(obj.filename,'>ZXX.VAR'); end
		try [obj.dzxy] = readedi(obj.filename,'>ZXY.VAR'); end
		try [obj.dzyx] = readedi(obj.filename,'>ZYX.VAR'); end
		try [obj.dzyy] = readedi(obj.filename,'>ZYY.VAR'); end
		try [obj.tzxr] = readedi(obj.filename,'>TXR.EXP'); end
		try [obj.tzxi] = readedi(obj.filename,'>TXI.EXP'); end
		try [obj.tzyr] = readedi(obj.filename,'>TYR.EXP'); end
		try [obj.tzyi] = readedi(obj.filename,'>TYI.EXP'); end
		try [obj.dtzx] = readedi(obj.filename,'>TXVAR.EXP'); end
		try [obj.dtzy] = readedi(obj.filename,'>TYVAR.EXP'); end
		try [obj.rhoxx ] = readedi(obj.filename,'>RHOXX'); end
		try [obj.rhoxy ] = readedi(obj.filename,'>RHOXY'); end
		try [obj.rhoyx ] = readedi(obj.filename,'>RHOYX'); end
		try [obj.rhoyy ] = readedi(obj.filename,'>RHOYY'); end
		try [obj.phsxx ] = readedi(obj.filename,'>PHSXX'); end
		try [obj.phsxy ] = readedi(obj.filename,'>PHSXY'); end
		try [obj.phsyx ] = readedi(obj.filename,'>PHSYX'); end
		try [obj.phsyy ] = readedi(obj.filename,'>PHSYY'); end
		try [obj.drhoxx] = readedi(obj.filename,'>RHOXX.ERR'); end
		try [obj.drhoxy] = readedi(obj.filename,'>RHOXY.ERR'); end
		try [obj.drhoyx] = readedi(obj.filename,'>RHOYX.ERR'); end
		try [obj.drhoyy] = readedi(obj.filename,'>RHOYY.ERR'); end
		try [obj.dphsxx] = readedi(obj.filename,'>PHSXX.ERR'); end
		try [obj.dphsxy] = readedi(obj.filename,'>PHSXY.ERR'); end
		try [obj.dphsyx] = readedi(obj.filename,'>PHSYX.ERR'); end
		try [obj.dphsyy] = readedi(obj.filename,'>PHSYY.ERR'); end
		% read frequency selection info.
		FSELFILE = [obj.filename '.fsel'];
		if exist(FSELFILE,'file')
			try [obj.fsel] = readedi(FSELFILE,'>FSEL'); end
			try [obj.zxxs] = readedi(FSELFILE,'>ZXXS'); end
			try [obj.zxys] = readedi(FSELFILE,'>ZXYS'); end
			try [obj.zyxs] = readedi(FSELFILE,'>ZYXS'); end
			try [obj.zyys] = readedi(FSELFILE,'>ZYYS'); end
			try [obj.tzxs] = readedi(FSELFILE,'>TZXS'); end
			try [obj.tzys] = readedi(FSELFILE,'>TZYS'); end
			try [obj.rhoxxs] = readedi(FSELFILE,'>RHOXXS'); end
			try [obj.rhoxys] = readedi(FSELFILE,'>RHOXYS'); end
			try [obj.rhoyxs] = readedi(FSELFILE,'>RHOYXS'); end
			try [obj.rhoyys] = readedi(FSELFILE,'>RHOYYS'); end
			try [obj.phsxxs] = readedi(FSELFILE,'>PHSXXS'); end
			try [obj.phsxys] = readedi(FSELFILE,'>PHSXYS'); end
			try [obj.phsyxs] = readedi(FSELFILE,'>PHSYXS'); end
			try [obj.phsyys] = readedi(FSELFILE,'>PHSYYS'); end
		end % if.
		obj = obj.csel2asel;
	end % read.
%
%  write edi.
%
	function write(obj,varargin)
		if (length(varargin))
			OUTFILE = varargin{1};
		else
			if isempty(obj.filename)
				obj.filename = [obj.dataid,'.edi'];
			end % if.
			OUTFILE = obj.filename;
		end % if.
		writeedi(OUTFILE,'>HEAD',obj.dataid,obj.nfreq,obj.lat,obj.lon,obj.elev);
		if (~isempty(obj.freq)) writeedi(OUTFILE,'>FREQ',obj.freq); end
		if (~isempty(obj.zxxr)) writeedi(OUTFILE,'>ZXXR',obj.zxxr); end
		if (~isempty(obj.zxxi)) writeedi(OUTFILE,'>ZXXI',obj.zxxi); end
		if (~isempty(obj.zxyr)) writeedi(OUTFILE,'>ZXYR',obj.zxyr); end
		if (~isempty(obj.zxyi)) writeedi(OUTFILE,'>ZXYI',obj.zxyi); end
		if (~isempty(obj.zyxr)) writeedi(OUTFILE,'>ZYXR',obj.zyxr); end
		if (~isempty(obj.zyxi)) writeedi(OUTFILE,'>ZYXI',obj.zyxi); end
		if (~isempty(obj.zyyr)) writeedi(OUTFILE,'>ZYYR',obj.zyyr); end
		if (~isempty(obj.zyyi)) writeedi(OUTFILE,'>ZYYI',obj.zyyi); end
		if (~isempty(obj.dzxx)) writeedi(OUTFILE,'>ZXX.VAR',obj.dzxx); end
		if (~isempty(obj.dzxy)) writeedi(OUTFILE,'>ZXY.VAR',obj.dzxy); end
		if (~isempty(obj.dzyx)) writeedi(OUTFILE,'>ZYX.VAR',obj.dzyx); end
		if (~isempty(obj.dzyy)) writeedi(OUTFILE,'>ZYY.VAR',obj.dzyy); end
		if (~isempty(obj.tzxr)) writeedi(OUTFILE,'>TXR.EXP',obj.tzxr); end
		if (~isempty(obj.tzxi)) writeedi(OUTFILE,'>TXI.EXP',obj.tzxi); end
		if (~isempty(obj.tzyr)) writeedi(OUTFILE,'>TYR.EXP',obj.tzyr); end
		if (~isempty(obj.tzyi)) writeedi(OUTFILE,'>TYI.EXP',obj.tzyi); end
		if (~isempty(obj.dtzx)) writeedi(OUTFILE,'>TXVAR.EXP',obj.dtzx); end
		if (~isempty(obj.dtzy)) writeedi(OUTFILE,'>TYVAR.EXP',obj.dtzy); end
		if (~isempty(obj.rhoxx )) writeedi(OUTFILE,'>RHOXX',obj.rhoxx ); end
		if (~isempty(obj.rhoxy )) writeedi(OUTFILE,'>RHOXY',obj.rhoxy ); end
		if (~isempty(obj.rhoyx )) writeedi(OUTFILE,'>RHOYX',obj.rhoyx ); end
		if (~isempty(obj.rhoyy )) writeedi(OUTFILE,'>RHOYY',obj.rhoyy ); end
		if (~isempty(obj.phsxx )) writeedi(OUTFILE,'>PHSXX',obj.phsxx ); end
		if (~isempty(obj.phsxy )) writeedi(OUTFILE,'>PHSXY',obj.phsxy ); end
		if (~isempty(obj.phsyx )) writeedi(OUTFILE,'>PHSYX',obj.phsyx ); end
		if (~isempty(obj.phsyy )) writeedi(OUTFILE,'>PHSYY',obj.phsyy ); end
		if (~isempty(obj.drhoxx)) writeedi(OUTFILE,'>RHOXX.ERR',obj.drhoxx); end
		if (~isempty(obj.drhoxy)) writeedi(OUTFILE,'>RHOXY.ERR',obj.drhoxy); end
		if (~isempty(obj.drhoyx)) writeedi(OUTFILE,'>RHOYX.ERR',obj.drhoyx); end
		if (~isempty(obj.drhoyy)) writeedi(OUTFILE,'>RHOYY.ERR',obj.drhoyy); end
		if (~isempty(obj.dphsxx)) writeedi(OUTFILE,'>PHSXX.ERR',obj.dphsxx); end
		if (~isempty(obj.dphsxy)) writeedi(OUTFILE,'>PHSXY.ERR',obj.dphsxy); end
		if (~isempty(obj.dphsyx)) writeedi(OUTFILE,'>PHSYX.ERR',obj.dphsyx); end
		if (~isempty(obj.dphsyy)) writeedi(OUTFILE,'>PHSYY.ERR',obj.dphsyy); end
	end % write.
%
%  write edi.
%
	function writefsel(obj,varargin)
		if (length(varargin) > 0)
			FSELFILE = varargin{1};
		else
			FSELFILE = [obj.filename '.fsel'];
		end % if.
		writeedi(FSELFILE,'>HEAD',obj.dataid,obj.nfreq,obj.lat,obj.lon,obj.elev);
		if (~isempty(obj.fsel)) writeedi(FSELFILE,'>FSEL',obj.fsel); end
		if (~isempty(obj.zxxs)) writeedi(FSELFILE,'>ZXXS',obj.zxxs); end
		if (~isempty(obj.zxys)) writeedi(FSELFILE,'>ZXYS',obj.zxys); end
		if (~isempty(obj.zyxs)) writeedi(FSELFILE,'>ZYXS',obj.zyxs); end
		if (~isempty(obj.zyys)) writeedi(FSELFILE,'>ZYYS',obj.zyys); end
		if (~isempty(obj.tzxs)) writeedi(FSELFILE,'>TZXS',obj.tzxs); end
		if (~isempty(obj.tzys)) writeedi(FSELFILE,'>TZYS',obj.tzys); end
		if (~isempty(obj.rhoxxs)) writeedi(FSELFILE,'>RHOXXS',obj.rhoxxs); end
		if (~isempty(obj.rhoxys)) writeedi(FSELFILE,'>RHOXYS',obj.rhoxys); end
		if (~isempty(obj.rhoyxs)) writeedi(FSELFILE,'>RHOYXS',obj.rhoyxs); end
		if (~isempty(obj.rhoyys)) writeedi(FSELFILE,'>RHOYYS',obj.rhoyys); end
		if (~isempty(obj.phsxxs)) writeedi(FSELFILE,'>PHSXXS',obj.phsxxs); end
		if (~isempty(obj.phsxys)) writeedi(FSELFILE,'>PHSXYS',obj.phsxys); end
		if (~isempty(obj.phsyxs)) writeedi(FSELFILE,'>PHSYXS',obj.phsyxs); end
		if (~isempty(obj.phsyys)) writeedi(FSELFILE,'>PHSYYS',obj.phsyys); end
	end % write.
%
%  complete edi.
%
	function obj = complete(obj)
		if (isempty(obj.asel)) obj.asel = ones(obj.nfreq,10);end
		if (isempty(obj.fsel)) obj.fsel = ones(obj.nfreq,1); end
		if (isempty(obj.zxxs)) obj.zxxs = ones(obj.nfreq,1); end
		if (isempty(obj.zxys)) obj.zxys = ones(obj.nfreq,1); end
		if (isempty(obj.zyxs)) obj.zyxs = ones(obj.nfreq,1); end
		if (isempty(obj.zyys)) obj.zyys = ones(obj.nfreq,1); end
		if (isempty(obj.tzxs)) obj.tzxs = ones(obj.nfreq,1); end
		if (isempty(obj.tzys)) obj.tzys = ones(obj.nfreq,1); end
		if (isempty(obj.rhoxys)) obj.rhoxys = ones(obj.nfreq,1); end
		if (isempty(obj.rhoyxs)) obj.rhoyxs = ones(obj.nfreq,1); end
		if (isempty(obj.phsxys)) obj.phsxys = ones(obj.nfreq,1); end
		if (isempty(obj.phsyxs)) obj.phsyxs = ones(obj.nfreq,1); end
		if (isempty(obj.rhoxx)) obj.rhoxx = classEDI.comprho(complex(obj.zxxr,obj.zxxi),obj.freq); end
		if (isempty(obj.rhoxy)) obj.rhoxy = classEDI.comprho(complex(obj.zxyr,obj.zxyi),obj.freq); end
		if (isempty(obj.rhoyx)) obj.rhoyx = classEDI.comprho(complex(obj.zyxr,obj.zyxi),obj.freq); end
		if (isempty(obj.rhoyy)) obj.rhoyy = classEDI.comprho(complex(obj.zyyr,obj.zyyi),obj.freq); end
		if (isempty(obj.phsxx)) obj.phsxx = classEDI.compphs(complex(obj.zxxr,obj.zxxi)); end
		if (isempty(obj.phsxy)) obj.phsxy = classEDI.compphs(complex(obj.zxyr,obj.zxyi)); end
		if (isempty(obj.phsyx)) obj.phsyx = classEDI.compphs(complex(obj.zyxr,obj.zyxi)); obj.phsyx = obj.phsyx + 180; end
		if (isempty(obj.phsyy)) obj.phsyy = classEDI.compphs(complex(obj.zyyr,obj.zyyi)); end
		if (isempty(obj.drhoxx)) obj.drhoxx = classEDI.compdrho(complex(obj.zxxr,obj.zxxi),obj.dzxx,obj.freq); end
		if (isempty(obj.drhoxy)) obj.drhoxy = classEDI.compdrho(complex(obj.zxyr,obj.zxyi),obj.dzxy,obj.freq); end
		if (isempty(obj.drhoyx)) obj.drhoyx = classEDI.compdrho(complex(obj.zyxr,obj.zyxi),obj.dzyx,obj.freq); end
		if (isempty(obj.drhoyy)) obj.drhoyy = classEDI.compdrho(complex(obj.zyyr,obj.zyyi),obj.dzyy,obj.freq); end
		if (isempty(obj.dphsxx)) obj.dphsxx = classEDI.compdphs(complex(obj.zxxr,obj.zxxi),obj.dzxx); end
		if (isempty(obj.dphsxy)) obj.dphsxy = classEDI.compdphs(complex(obj.zxyr,obj.zxyi),obj.dzxy); end
		if (isempty(obj.dphsyx)) obj.dphsyx = classEDI.compdphs(complex(obj.zyxr,obj.zyxi),obj.dzyx); end
		if (isempty(obj.dphsyy)) obj.dphsyy = classEDI.compdphs(complex(obj.zyyr,obj.zyyi),obj.dzyy); end
	end % complete.
%
%  convert sel matrix to sel components.
%
	function obj = asel2csel(obj)
		if (~isempty(obj.asel))
			obj.rhoxys = obj.asel(:,1);
			obj.rhoyxs = obj.asel(:,2);
			obj.phsxys = obj.asel(:,3);
			obj.phsyxs = obj.asel(:,4);
			obj.zxxs   = obj.asel(:,5);
			obj.zxys   = obj.asel(:,6);
			obj.zyxs   = obj.asel(:,7);
			obj.zyys   = obj.asel(:,8);
			obj.tzxs   = obj.asel(:,9);
			obj.tzys   = obj.asel(:,10);
		end % if.
	end % asel2csel.
%
%  convert sel components to sel matrix.
%
	function obj = csel2asel(obj)
		if ~isempty(obj.rhoxys) obj.asel(:,1) = obj.rhoxys; end
		if ~isempty(obj.rhoyxs) obj.asel(:,2) = obj.rhoyxs; end
		if ~isempty(obj.phsxys) obj.asel(:,3) = obj.phsxys; end
		if ~isempty(obj.phsyxs) obj.asel(:,4) = obj.phsyxs; end
		if ~isempty(obj.zxxs) obj.asel(:,5) = obj.zxxs; end
		if ~isempty(obj.zxys) obj.asel(:,6) = obj.zxys; end
		if ~isempty(obj.zyxs) obj.asel(:,7) = obj.zyxs; end
		if ~isempty(obj.zyys) obj.asel(:,8) = obj.zyys; end
		if ~isempty(obj.tzxs) obj.asel(:,9) = obj.tzxs; end
		if ~isempty(obj.tzys) obj.asel(:,10)= obj.tzys; end
	end % csel2asel.
%
%  plot.
%
	function plot(obj,TYPE)
		switch TYPE
		case 'rhophi'
			hf = figure('Position',[10,10,500,800],'PaperPosition',[2 2 5 8]); clf;
			subplot(211);
			loglog(obj.freq,obj.rhoxy,'o','color',[0 .5 .8],'linewidth',1.5);
			hold on
			errorbar(obj.freq,obj.rhoxy,2*obj.drhoxy,'color',[0 .5 .8],'linewidth',1.5);
			loglog(obj.freq,obj.rhoyx,'rs','linewidth',1.5);
			errorbar(obj.freq,obj.rhoyx,2*obj.drhoyx,'rs','linewidth',1.5);
			hold off
			subplot(212)
			semilogx(obj.freq,obj.phsxy,'o','color',[0 .5 .8],'linewidth',1.5);
			hold on
			errorbar(obj.freq,obj.phsxy,2*obj.dphsxy,'color',[0 .5 .8],'linewidth',1.5);
			semilog(obj.freq,obj.phsyx,'rs','linewidth',1.5);
			errorbar(obj.freq,obj.phsyx,2*obj.dphsyx,'rs','linewidth',1.5);
			hold off
		case 'impedance'
		case 'vtf'
		otherwise
		end % type.
	end % plot


%----------------------------------------------------------------------%
%  end of methods definition.                                          %
%----------------------------------------------------------------------%
	end % methods.
%----------------------------------------------------------------------%
%  static methods.                                                     %
%----------------------------------------------------------------------%
	methods (Static)
%
%  compute rho.
%
	function rho = comprho(z,freq)
		% covert units from [mV/km]/[nT] to [V/A].
		z = z * 4.e-4 * pi;
		% rho = abs(z) / (2*pi*f*mu);
		rho = abs(z).^2 ./ freq / 8.e-7 / pi / pi;
	end % comprho.
%
%  compute drho.
%
	function drho = compdrho(z,dz,freq)
		% drho = 2*mu/(2*pi*f) * abs(z) * sqrt(dz) * 1e6;
		%drho = 0.4 ./ freq .* abs(z) .* sqrt(dz);
		drho = 0.4 ./ freq .* abs(z) .* (dz);
	end % compdrho.
%
%  compute phs.
%
	function phs = compphs(z)
		% phs = atan2(imag(z),real(z));
		phs = atan2(imag(z),real(z)) * 180.0 / pi;
	end % compphs.
%
%  compute dphs.
%
	function dphs = compdphs(z,dz)
		% dphs = sqrt(dz) / abs(z) * 180 / pi;
		%dphs = sqrt(dz) ./ abs(z) * 180.0 / pi;
		dphs = (dz) ./ abs(z) * 180.0 / pi;
	end % compdphs.
	end % static methods.
%
% End of the class.
%
end % classEDI.

