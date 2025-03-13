classdef ListData < handle
% Railo 2014
%
%

properties
  Lat, Long, Z   % number of grid cells in 
  sigma          % 3D array of sigma coductivities (log_e scale)
  xUTM, yUTM, zoneUTM % length of the list
  N              % length of the list
end

methods
%*******************************************************************
function obj = ListData()
  %   class constructor ... simple
end
%*******************************************************************
%*******************************************************************

function makeUTMList(obj)
    [obj.xUTM, obj.yUTM, obj.zoneUTM] = deg2utm(obj.Lat,obj.Long);
end

function restrictZone(obj,zone)
    counter = 0;
    for i = 1:obj.N
        temp = obj.zoneUTM(i,:);
        if temp(1) == zone(1) && temp(2) == zone(2) ...
                && double(temp(4)) == double(zone(4))
            counter = counter + 1;
            xTemp(counter) = obj.xUTM(i);
            yTemp(counter) = obj.yUTM(i);
            zTemp(counter) = obj.Z(i);
            sTemp(counter) = obj.sigma(i);
            zoneTemp(counter,:) = temp;
        end
    end
    obj.xUTM = xTemp'; obj.yUTM = yTemp'; obj.Z = zTemp';
    obj.sigma = sTemp'; obj.zoneUTM = zoneTemp; obj.N = length(xTemp);
    [obj.Lat,obj.Long] = utm2deg(xTemp,yTemp,zoneTemp);
end
    
end     % methods
end    % classdef