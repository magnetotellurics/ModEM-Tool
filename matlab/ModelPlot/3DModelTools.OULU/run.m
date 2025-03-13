
% this is netcdf fiel format
% available for all SMAP layear

fname = 'data/seawater.grd'  %NETCDF file
%fname = 'data/ETOPO1_Ice_g_gdal.grd';

% shape fiels can be used instead but 
% might need to be converted to a proper format


% TMaxwell object should probably always be called
% to be created from file, so input argument
% would be a file
% however we need it for multi grid.....

m = TBGModel(fname)
% load lat.long file
% project and make a grid
% should be compibet eith TGrid3D object, just an example
%m.LoadFromXYZFile(fname);
%m.PlotPlainView(1);




% load from .grd file
m.LoadFromNETCDFFile(fname);
m.PlotPlainView(1);



