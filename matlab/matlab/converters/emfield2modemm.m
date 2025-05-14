function objOut = emfield2modemm(objIn)
% objOut = emfield2modemm(objIn)
%
% converter from ModelPlot emfieldplot to ModEMM TVector3D_SG or TVector3D_Global

if ~isa(objIn,'emfieldplot')
    error('Input to emfield2modemm should be of ModelPlot emfieldplot derivative class');
end

if isa(objIn.grid,'xygrid')
    disp('Input to emfield2modemm is in cartesian coordinates. Using xy type');
    grid = xygrid2modemm(objIn.grid);
elseif isa(objIn.grid,'llgrid')
    disp('Input to emfield2modemm is in spherical coordinates. Using latlon type');
    grid = llgrid2modemm(objIn.grid);
end

% call emfieldplot constructor
edgeORface = lower(objIn.location);
if objIn.isglobal
    objOut = TVector3D_Global(grid,edgeORface);
else
    objOut = TVector3D_SG(grid,edgeORface);
end

% update values
objOut.x = objIn.x;
objOut.y = objIn.y;
objOut.z = objIn.z;

% quick final check
if isa(objOut,'TVector3D_SG') || isa(objOut,'TVector3D_Global')
    disp('Output EM field is in the format ModEMM TVector3D derivative.');
else
    error('Error converting ModelPlot emfieldplot to ModEMM TVector3D');
end
