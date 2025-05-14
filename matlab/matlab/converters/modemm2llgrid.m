function objOut = modemm2llgrid(objIn,lat0,lon0)
% objOut = modemm2llgrid(objIn,lat0,lon0)
%
% converter from ModEMM TGrid3D to ModelPlot llgrid
% lat0,lon0 are optional for when the input is a cartesian grid

if ~isa(objIn,'TGrid3D')
    error('Input to modemm2llgrid should be of ModEMM TGrid3D class');
end

if strcmp(objIn.type,'xy')
    disp('Input to modemm2llgrid is in cartesian coordinates. Converting to llgrid...');
    obj = modemm2xygrid(objIn);
    objOut = llgrid(obj,lat0,lon0);
    return
end

% RELOCATE THE GRID TO BE CENTERED AROUND LAT0 & LON0
if nargin >= 3
    lowerleftcorner = [lat0 lon0 0] - [sum(objIn.Dx)/2 sum(objIn.Dy)/2 0];
else
    lowerleftcorner = objIn.origin;
end    

% compute zAir from objIn.Dz(1:objIn.Nza)
zAir(1:objIn.Nza+1,1) = NaN;
zAir(1) = 0;
zAir(2:end) = cumsum(objIn.Dz(objIn.Nza:-1:1));

% define grid limits - note a potential problem: origin always defines the
% lower left corner, otherwise can't get lat/lon from dlat/dlon - but
% that's not a good location for the origin in coordinate conversions!
lims.latmin = lowerleftcorner(1);
lims.latmax = lowerleftcorner(1) + sum(objIn.Dx);
lims.lonmin = lowerleftcorner(2);
lims.lonmax = lowerleftcorner(2) + sum(objIn.Dy);
lims.depthmin = lowerleftcorner(3);
lims.depthmax = lowerleftcorner(3) + sum(objIn.Dz(objIn.Nza+1:end));

% call llgrid constructor
objOut = llgrid(lims);

% set the array sizes and the air
objOut.nlat = objIn.Nx;
objOut.nlon = objIn.Ny;
objOut.nzEarth = objIn.NzEarth;
objOut.nzAir = objIn.Nza;
objOut.zAir = zAir;

% set some metadata
objOut.padding = '';
objOut.logdz = [];

% set the values
objOut.dlat = objIn.Dx;
objOut.dlon = objIn.Dy;
objOut.dz = objIn.Dz(objIn.Nza+1:end);
objOut.lat = [0; cumsum(objOut.dlat)] + lowerleftcorner(1);
objOut.lon = [0; cumsum(objOut.dlon)] + lowerleftcorner(2);
objOut.depth = [0; cumsum(objOut.dz)] + lowerleftcorner(3);

% convert from meters to km
objOut.units = objIn.units;
if strcmp(objOut.units,'m')
    objOut.dz = objOut.dz * 1e-3; % km
    objOut.depth = objOut.depth * 1e-3; % km
    objOut.zAir = objOut.zAir * 1e-3; % km
    objOut.units = 'km';
end

% quick final check
if isa(objOut,'llgrid')
    disp('Output grid object is in the format ModelPlot llgrid.');
else
    error('Error converting ModEMM TGrid3D to ModelPlot llgrid');
end
