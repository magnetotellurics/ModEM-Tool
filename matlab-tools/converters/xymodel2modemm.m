function objOut = xymodel2modemm(objIn)
% objOut = xymodel2modemm(objIn)
%
% converter from ModelPlot xymodel to ModEMM TModelParameterCell3D (in meters)

if ~isa(objIn,'xymodel')
    error('Input to xymodel2modemm should be of ModelPlot xymodel class');
end

% call the constructor for cartesian coordinate model
objOut = TModelParametersCell3D_EC_SG();
objOut.ParamGrid = xygrid2modemm(objIn.grid);

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
    disp('Output electrical conductivity model is in the format ModEMM TModelParameterCell3D cartesian.');
else
    error('Error converting ModelPlot xymodel to ModEMM TModelParameterCell3D_EC_SG');
end
