function objOut = modemm2emfield(objIn)
% objOut = modemm2emfield(objIn)
%
% converter from ModEMM TVector3D_SG or TVector3D_Global to ModelPlot emfield

if ~isa(objIn,'TVector3D_SG') && ~isa(objIn,'TVector3D_Global')
    error('Input to modemm2emfield should be of ModEMM TVector3D derivative class');
end

if strcmp(objIn.grid.type,'xy')
    disp('Input to modemm2emfield is in cartesian coordinates. Using xygrid');
    grid = modemm2xygrid(objIn.grid);
elseif strcmp(objIn.grid.type,'latlon')
    disp('Input to modemm2emfield is in spherical coordinates. Using llgrid');
    grid = modemm2llgrid(objIn.grid);
end

% call emfieldplot constructor
edgeORface = upper(objIn.type);
objOut = emfieldplot(grid,edgeORface);

% update values
objOut.x = objIn.x;
objOut.y = objIn.y;
objOut.z = objIn.z;

% quick final check
if isa(objOut,'emfieldplot')
    disp('Output EM field is in the format ModelPlot emfieldplot.');
else
    error('Error converting ModEMM TVector3D to ModelPlot emfieldplot');
end
