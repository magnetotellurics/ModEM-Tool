function objOut = modemm2xygrid(objIn)
% objOut = modemm2xygrid(objIn)
%
% converter from ModEMM TGrid3D to ModelPlot xygrid

if ~isa(objIn,'TGrid3D')
    error('Input to modemm2xygrid should be of ModEMM TGrid3D class');
end

if strcmp(objIn.type,'latlon')
    disp('Input to modemm2xygrid is in spherical coordinates. Converting to xygrid...');
    obj = modemm2llgrid(objIn);
    objOut = xygrid(obj);
    return
end

% compute zAir from objIn.Dz(1:objIn.Nza)
zAir(1:objIn.Nza+1,1) = NaN;
zAir(1) = 0;
zAir(2:end) = cumsum(objIn.Dz(objIn.Nza:-1:1));

% call xygrid constructor
objOut = xygrid(objIn.Dx,objIn.Dy,objIn.Dz(objIn.Nza+1:end),...
    objIn.origin,objIn.rotation,objIn.units,objIn.Nza,zAir);

% convert from meters to km
if strcmp(objOut.units,'m')
    objOut.origin = objOut.origin * 1e-3; % km
    objOut.dx = objOut.dx * 1e-3; % km
    objOut.dy = objOut.dy * 1e-3; % km
    objOut.dz = objOut.dz * 1e-3; % km
    objOut.zAir = objOut.zAir * 1e-3; % km
    objOut.units = 'km';
end

% quick final check
if isa(objOut,'xygrid')
    disp('Output grid object is in the format ModelPlot xygrid.');
else
    error('Error converting ModEMM TGrid3D to ModelPlot xygrid');
end
