classdef T2DModel < ListData
% Smirnov 2013, Railo 2014
% class of 2-D data < ListData
% class variables: strikeAngle
% class (user) fuctions: giveLine, findSide, findXY, findZ

properties
    % strike angle of the measurement
    strikeAngle
end

methods
%*******************************************************************
function obj = T2DModel(Filename,varargin)
 %   class constructor ... simple
 if nargin > 1
    obj.strikeAngle = varargin{1}
 else
    obj.strikeAngle = 90;
 end
    obj.LoadFromFile(Filename);
end

function LoadFromFile(obj, Filename)
    % Leads data from the file
    Filename = 'data/Trondelag-Jamtland2008-model.dat'; % for testing, should be removed
    fid = fopen(Filename);
    a = textscan(fid, '%f64%f64%f64%f64','headerlines',0);
    fclose(fid);
    a = cell2mat(a);
    obj.Long = a(:,1);
    obj.Lat = a(:,2);
    obj.Z = 1000*a(:,3);         
    % depth in m, originally in km
    obj.sigma = log(10.^(-a(:,4)));     
    % conductivities, originally in lg(resitivity) -form
    obj.N = length(obj.sigma); % length of the list
    obj.makeUTMList;
end;

function [A,B] = giveLine(obj)
    % point is structure with .x and .y    
    xmin = min(obj.xUTM);
    xmax = max(obj.xUTM);
    ymin = min(obj.yUTM);
    ymax = max(obj.yUTM);
    XX = find(obj.xUTM==xmin);
    YY = find(obj.yUTM==ymin);
  
    %   assume that the line is not perpendicular or parallel to x-axis
    if XX(1) == YY(1)
        xend = xmin; yend = ymin;
        xbeg = xmax; ybeg = ymax;
    elseif XX(1) ~= YY(1)
        xend = xmin; yend = ymax;
        xbeg = xmax; ybeg = ymin;
    else
        disp('Line is perpendicular or parallel to an axis');
    end
    A = [xbeg,ybeg];
    B = [xend,yend];
end

function [A,B] = giveLineLatLong(obj)
    xmin = min(obj.Lat);
    xmax = max(obj.Lat);
    ymin = min(obj.Long);
    ymax = max(obj.Long);
    XX = find(obj.Lat==xmin);
    YY = find(obj.Long==ymin);
  
    %   assume that the line is not perpendicular or parallel to x-axis
    if XX(1) == YY(1)
        xend = xmin; yend = ymin;
        xbeg = xmax; ybeg = ymax;
    elseif XX(1) ~= YY(1)
        xend = xmin; yend = ymax;
        xbeg = xmax; ybeg = ymin;
    else
        disp('Line is perpendicular or parallel to an axis');
    end
    A = [xbeg,ybeg];
    B = [xend,yend];
end
    

function res = findSide(obj,Q)
    % Finds out which side of the line AB is the point Q.
    %
    % INPUT
    % Q, A, B 2d vectors
    %
    % OUTPUT
    % -1 = "under" the line
    %  0 = on the line
    %  1 = "over" the line

    [A,B] = giveLine(obj);
    area = (B(1)-A(1))*(Q(2)-A(2))-(B(2)-A(2))*(Q(1)-A(1));
    res = sign(area);
end

function res = findXY(obj,Q)
    Angle = obj.strikeAngle;
    [A,B] = giveLine(obj);
    side = findSide(obj,Q);
    
    XQ(:,1) = Q(1)-obj.xUTM;
    XQ(:,2) = Q(2)-obj.yUTM; 
    temp = XQ.^2;
    normQ = sqrt(sum(temp'));
    
     if side == -1
        XA(:,1) = A(1)-obj.xUTM;
        XA(:,2) = A(2)-obj.yUTM;
        temp = XQ.*XA;
        dot = sum(temp');
        temp = XA.^2;
        normi = sqrt(sum(temp'));
     elseif side == 1
         XB(:,1) = B(1)-obj.xUTM;
         XB(:,2) = B(2)-obj.yUTM;
         temp = XQ.*XB;
         dot = sum(temp');
        temp = XB.^2;
         normi = sqrt(sum(temp'));
     end
    testValue = abs(dot./normi-normQ.*cos(deg2rad(Angle)));
    minim = min(testValue);
    res = find(testValue == minim);
    % not unique; gives a line with fixed (x,y) coordinates
end

function res = findZ(obj,indmin,z)
    ind = findNearest(z,obj.Z(indmin));
    if length(ind) > 1
        disp('The index was not unique!')
    end
    res = indmin(ind(1));
end

function res = findLatLong(obj,Q)
    Angle = obj.strikeAngle;
    [A,B] = obj.giveLineLatLong;
    
    % See "help SLoC"
    testValue = abs(SLoC([obj.Lat obj.Long],B,Q)-cos(deg2rad(Angle)));
    minim = min(testValue);
    res = find(testValue == minim);
    % not unique; gives a line with fixed (Lat,Long) coordinates
    
end

% Plot functions

function plotLayer(obj, depth)
    h = pcolor(obj.grid.x,obj.grid.y,obj.sig(:,:,depth));
    set(h, 'EdgeColor', 'none');
    colorbar
end

function plotLinedensities(obj, depth)
    h = pcolor(obj.grid.x,obj.grid.y,obj.sig(:,:,depth));
    set(h, 'EdgeColor', 'none');
    colorbar
end



%*******************************************************************

end     % methods
end    % classdef