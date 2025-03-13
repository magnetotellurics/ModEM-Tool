classdef T1DModel < TModel
% Smirnov 2013   
%

properties
end

methods
%*******************************************************************
function obj = T1DModel(Filename)
 %   class constructor ... simple
  obj.LoadFromFile(Filename);
end

function LoadFromFile(obj, Filename)

end;


%*******************************************************************

end     % methods
end    % classdef