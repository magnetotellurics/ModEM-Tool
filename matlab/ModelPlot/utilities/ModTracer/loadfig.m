function [hr] = loadfig(FIGFILE,hdl)
% loadfig - function to load *.fig file into figure handle.
%
%  NOTES: 
%  
%
%  See also 
%  --------------------
%  Bo Yang, 2014.
%  Institute of Geophysics and Geomatics,
%  China University of Geosciences, Wuhan, China.
%  Comments, bug reports and questions, please send to:
%  yangbo.cug@163.com.
%  Copyright 2014-2018, Bo Yang, IGG, CUG.
%  $Revision: 1.0 $ $Date: 2014/02/26 $
%  Last changed: 2014/02/26 14:17:02.

%  Revision log:
%  2014/02/26 : Version 1.0 released.
% clf
hf = hgload(FIGFILE,struct('visible','off'));
ha = get(hf,'Children'); 
hr = copyobj(ha,hdl);
delete(hf);
%
% End of the function.
%
