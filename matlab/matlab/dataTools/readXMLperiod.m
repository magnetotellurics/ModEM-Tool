function varargout = readXMLperiod(t,np,type)
	periods = find(t,'/EM_TF/Data/Period');
	switch type
		case 'Z'
			uid = find(t,periods(np),'name','Z');
			val = find(t,uid,'name','value');
			Zxx = getZTvalue(t,val,'Zxx');
			Zxy = getZTvalue(t,val,'Zxy');
			Zyx = getZTvalue(t,val,'Zyx');
			Zyy = getZTvalue(t,val,'Zyy');
			varargout{1} = [Zxx Zxy Zyx Zyy];
			return;
		case 'Z.VAR'
			uid = find(t,periods(np),'name','Z.VAR');
			val = find(t,uid,'name','value');
			Zxx = getZTvalue(t,val,'Zxx');
			Zxy = getZTvalue(t,val,'Zxy');
			Zyx = getZTvalue(t,val,'Zyx');
			Zyy = getZTvalue(t,val,'Zyy');
			varargout{1} = [Zxx Zxy Zyx Zyy];
			return;
		case 'T'
			uid = find(t,periods(np),'name','T');
			val = find(t,uid,'name','value');
			Tx = getZTvalue(t,val,'Tx');
			Ty = getZTvalue(t,val,'Ty');
			varargout{1} = [Tx Ty];
			return;
		case 'T.VAR'
			uid = find(t,periods(np),'name','T.VAR');
			val = find(t,uid,'name','value');
			Tx = getZTvalue(t,val,'Tx');
			Ty = getZTvalue(t,val,'Ty');
			varargout{1} = [Tx Ty];
			return;
		otherwise
			disp(['unsupported type: ',type]);
	end % switch type.

	function val = getZTvalue(t,uid,kw)
		n = length(uid);
		for k = 1:n
			vname = attributes(t,'get',uid(k),1);
			vname = vname.val;
			if (vname == kw)
				str = get(t,uid(k));
				str = get(t,str.contents);
				str = str.value;
				num = str2num(str);
				if (length(num) == 2)
					val = num(1) + 1i*num(2);
				elseif (length(num) == 1)
					val = num;
				else
					error(['ERROR: the value of ',kw,'is either complex nor double!']);
				end % if lenght(num).
				return;
			end % if
		end % k
		error(['ERROR: Can not find the value of ',kw,'!']);
	end % function getZTvalue.
%
% end of function readXMLperiod.
%
end