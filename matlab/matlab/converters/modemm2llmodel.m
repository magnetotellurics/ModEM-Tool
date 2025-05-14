function objOut = modemm2llmodel(objIn,lat0,lon0)
% objOut = modemm2llmodel(objIn,lat0,lon0)
%
% converter from ModEMM TGrid3D to ModelPlot llmodel
% lat0,lon0 are optional for when the input is a cartesian grid

if ~isa(objIn,'TModelParameterCell3D_EC_SG') && ~isa(objIn,'TModelParameterCell3D_ES_SG')
    error('Input to modemm2llmodel should be of ModEMM TModelParameterCell3D derivative class');
end

if strcmp(objIn.ParamGrid.type,'xy')
    disp('Input to modemm2llmodel is in cartesian coordinates. Converting to llmodel...');
    obj = modemm2xymodel(objIn);
    objOut = llmodel(obj,lat0,lon0);
    return
elseif strcmp(objIn.ParamGrid.type,'latlon')
    grid = modemm2llgrid(objIn.ParamGrid);
end

% call llmodel constructor
value = objIn.v;
paramType = objIn.paramType;
modelType = 'electrical conductivity';
modelUnits = 'S/m';
objOut = llmodel(grid,value,paramType,modelType,modelUnits);

% update other values
objOut.AirCond = objIn.AirCond;

% quick final check
if isa(objOut,'llmodel')
    disp('Output electrical conductivity model is in the format ModelPlot llmodel.');
else
    error('Error converting ModEMM TModelParameterCell3D to ModelPlot llmodel');
end
