function objOut = llgrid2modemm(objIn)
% objOut = llgrid2modemm(objIn)
%
% converter from ModelPlot llgrid to ModEMM TGrid3D (in meters)

if ~isa(objIn,'llgrid')
    error('Input to llgrid2modemm should be of ModelPlot llgrid class');
end

Dx = objIn.dlat; % degrees
Dy = objIn.dlon; % degrees

if strcmp(objIn.units,'km')
    Dz = 1000*objIn.dz; % meters
    zorigin = 1000*objIn.depth(1);
else
    Dz = objIn.dz; % meters
    zorigin = objIn.depth(1);
end

% now update Dz with air layers if they exist
if ~isempty(objIn.zAir)
    Dz = [- diff(sort(objIn.zAir,1,'descend')); Dz];
end
Nza = max(length(objIn.zAir)-1,0);
objOut = TGrid3D(Dx,Dy,Dz,Nza);
objOut.units = 'm';
objOut.type = 'latlon';

% critically important: set the correct "origin" to lower left corner
objOut.origin(1) = objIn.lat(1);
objOut.origin(2) = objIn.lon(1);
objOut.origin(3) = zorigin;

% sometimes zAir is empty but nzAir just indicates how many air
% layers we want... here we accommodate this case
if Nza == 0
    Nza = objIn.nzAir;
    objOut = setAirLayers(objOut,'maxheight',1e+6,...
        'nlayers',Nza);%,'method','mirror','nlayers',10);
end

% quick final check
if isa(objOut,'TGrid3D')
    disp('Output grid object is in the format ModEMM TGrid3D spherical.');
else
    error('Error converting ModelPlot llgrid to ModEMM TGrid3D');
end
