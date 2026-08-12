function objOut = modemm2xymodel(objIn)
% objOut = modemm2xymodel(objIn)
%
% converter from ModEMM TGrid3D to ModelPlot xymodel
% lat0,lon0 are optional for when the input is a cartesian grid

if ~isa(objIn,'TModelParameterCell3D_EC_SG') && ~isa(objIn,'TModelParameterCell3D_ES_SG')
    error('Input to modemm2xymodel should be of ModEMM TModelParameterCell3D derivative class');
end

if strcmp(objIn.ParamGrid.type,'latlon')
    disp('Input to modemm2llmodel is in spherical coordinates. Converting to xymodel...');
    obj = modemm2llmodel(objIn);
    objOut = xymodel(obj);
    return
elseif strcmp(objIn.ParamGrid.type,'xy')
    grid = modemm2xygrid(objIn.ParamGrid);
end

% call llmodel constructor
value = objIn.v;
paramType = objIn.paramType;
modelType = 'electrical conductivity';
modelUnits = 'S/m';
objOut = xymodel(grid,value,paramType,modelType,modelUnits);

% update other values
objOut.AirCond = objIn.AirCond;

% quick final check
if isa(objOut,'xymodel')
    disp('Output electrical conductivity model is in the format ModelPlot xymodel.');
else
    error('Error converting ModEMM TModelParameterCell3D to ModelPlot xymodel');
end
