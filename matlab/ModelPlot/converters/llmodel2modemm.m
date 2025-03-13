function objOut = llmodel2modemm(objIn)
% objOut = llmodel2modemm(objIn)
%
% converter from ModelPlot llmodel to ModEMM TModelParameterCell3D (in meters)

if ~isa(objIn,'llmodel')
    error('Input to llmodel2modemm should be of ModelPlot llmodel class');
end

% call the constructor for spherical coordinate model
objOut = TModelParametersCell3D_ES_SG();
objOut.ParamGrid = llgrid2modemm(objIn.grid);

% convert to natural log in the cells for values
obj = objIn.loge;
if strcmp(obj.location,'NODE')
    obj = obj.node2cell;
end
objOut.paramType = obj.paramType;
objOut.AirCond = obj.AirCond;
objOut.v = obj.v;

% quick final check
if isa(objOut,'TModelParameterCell3D')
    disp('Output electrical conductivity model is in the format ModEMM TModelParameterCell3D spherical.');
else
    error('Error converting ModelPlot llmodel to ModEMM TModelParameterCell3D_ES_SG');
end
