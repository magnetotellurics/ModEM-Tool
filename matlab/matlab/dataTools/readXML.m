function [allData,units,isign,info] = readXML(cfile,newUnits,dataType)
%  Usage: [allData,units,isign,info] = readXML(cfile,newUnits,dataType);
%   read contents of cell array allData from direction [cfile] of all
%   XML files. There is one cell per period; each cell contains all
%   informaiton necessary to define data (location, values, error
%   standard dev) for each period (transmitter).
%
%  Written by Bo Yang, 2013.
%  Last change: 2013/11/12 11:21:52.

files = findfiles('xml',cfile);
nfile = length(files);
for i = 1:nfile
	disp(['Read file ',files{i},'...']);
	t = xmltree(files{i});
	% get number of period.
	uid = find(t,'//Data');
	nper= attributes(t,'get',uid(1),1);
	nper= str2num(nper.val);
	if (i == 1)
		% get units.
		uid   = find(t,'//DataType');
		units = attributes(t,'get',uid(1),5);
		units = units.val;
		% get isign.
		signstr = get(t,children(t,find(t,'//SignConvention')),'value');
		if findstr(signstr,'-')
			isign = -1;
		else
			isign = +1;
		end % isign.
		for k = 1:nper
			allData{k}.Cmplx          = 1;
			allData{k}.units          = units;
			allData{k}.nSite          = [];
			allData{k}.signConvention = isign;
			allData{k}.orient         = 0;
			switch strtrim(dataType)
			case 'Full_Impedance'
				allData{k}.type           = 'Full_Impedance';
				allData{k}.compChar       = ['ZXX';'ZXY';'ZYX';'ZYY'];
				allData{k}.nComp          = 8;
			case 'Off_Diagonal_Impedance'
				allData{k}.type           = 'Off_Diagonal_Impedance';
				allData{k}.compChar       = ['ZXY';'ZYX'];
				allData{k}.nComp          = 4;
			case 'Full_Vertical_Components'
				allData{k}.type           = 'Full_Vertical_Components';
				allData{k}.compChar       = ['TX';'TY'];
				allData{k}.nComp          = 4;
			otherwise
				disp('Unknown data type!');
				break;
			end % switch dataType.
			allData{k}.primaryCoords  = 'latlon';
			allData{k}.origin         = [];
			allData{k}.topography     = [];
			allData{k}.dataErrors     = [];
			allData{k}.EarthRad       = 6378.1;
			allData{k}.EarthEcc2      = 0.0067;
		end % for k.
		% assign info.
		info{1}.ncomp = allData{1}.nComp;
		info{1}.comp  = allData{1}.compChar;
	end % if i==1.

	% get site location.
	lat = get(t,children(t,find(t,'/EM_TF/Site/Location/Latitude')),'value');
	lat = str2num(lat);
	lon = get(t,children(t,find(t,'/EM_TF/Site/Location/Longitude')),'value');
	lon = str2num(lon);
	elev= get(t,children(t,find(t,'/EM_TF/Site/Location/Elevation')),'value');
	elev= str2num(elev);
	% get site name.
	Id  = get(t,children(t,find(t,'/EM_TF/Site/Id')),'value');
	% assign info.
	info{1}.lat(i)   = lat;
	info{1}.lon(i)   = lon;
	info{1}.loc(i,:) = [NaN NaN NaN];
	info{1}.code(i,:)= Id;

	uid = find(t,'//Period');
	for k = 1:nper
		% get current period.
		sper = attributes(t,'get',uid(k),1);
		allData{k}.T      = str2num(sper.val);
		% set site name.
		allData{k}.siteChar(i,:)  = Id;
		% get the data and error.
		switch strtrim(dataType)
		case 'Full_Impedance'
			allData{k}.Z(i,:)         = readXMLperiod(t,k,'Z');
			allData{k}.Zerr(i,:)      = readXMLperiod(t,k,'Z.VAR');
		case 'Off_Diagonal_Impedance'
			Zfull    = readXMLperiod(t,k,'Z');
			Zfullerr = readXMLperiod(t,k,'Z.VAR');
			allData{k}.Z(i,:)         = Zfull(2:3);
			allData{k}.Zerr(i,:)      = Zfullerr(2:3);
		case 'Full_Vertical_Components'
			allData{k}.Z(i,:)         = readXMLperiod(t,k,'T');
			allData{k}.Zerr(i,:)      = readXMLperiod(t,k,'T.VAR');
		otherwise
			disp(['Unknown data type: ',dataType,'!']);
			break;
		end % switch dataType.
		% set lat and lon.
		allData{k}.lat(i)         = lat;
		allData{k}.lon(i)         = lon;
		allData{k}.siteLoc(i,:)   = [NaN NaN NaN];
		% assign info.
		info{1}.data(i,k,:) = allData{k}.Z(i,:);
		info{1}.err(i,k,:)  = allData{k}.Zerr(i,:);
	end % k

	if (i == 1)
		% assign info.
		for k = 1:nper
			info{1}.per(k) = allData{k}.T;
		end % for k
	end % if i==1.
end % i
