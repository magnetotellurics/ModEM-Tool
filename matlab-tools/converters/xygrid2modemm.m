function objOut = xygrid2modemm(objIn)
% objOut = xygrid2modemm(objIn)
%
% converter from ModelPlot xygrid to ModEMM TGrid3D (in meters)
% need to be careful with zAir in meters vs km - review this later

if ~isa(objIn,'xygrid')
    error('Input to xygrid2modemm should be of ModelPlot xygrid class');
end

if strcmp(objIn.units,'km')
    Dx = 1000*objIn.dx; % meters
    Dy = 1000*objIn.dy; % meters
    Dz = 1000*objIn.dz; % meters
else
    Dx = objIn.dx; % meters
    Dy = objIn.dy; % meters
    Dz = objIn.dz; % meters
end

% now update Dz with air layers if they exist
if ~isempty(objIn.zAir)
    Dz = [- diff(sort(objIn.zAir,1,'descend')); Dz];
end
Nza = max(length(objIn.zAir)-1,0);
objOut = TGrid3D(Dx,Dy,Dz,Nza);
objOut.units = 'm';
objOut.type = 'xy';

% sometimes zAir is empty but nzAir just indicates how many air
% layers we want... here we accommodate this case
if Nza == 0
    Nza = objIn.nzAir;
    objOut = setAirLayers(objOut,'maxheight',1e+6,...
        'nlayers',Nza);%,'method','mirror','nlayers',10);
end

% quick final check
if isa(objOut,'TGrid3D')
    disp('Output grid object is in the format ModEMM TGrid3D cartesian.');
else
    error('Error converting ModelPlot xygrid to ModEMM TGrid3D');
end
