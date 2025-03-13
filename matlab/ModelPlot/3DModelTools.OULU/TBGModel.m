classdef TBGModel < TModel
% Smirnov 2013   
%

properties
  geodata;  
  Size   %4x2 array   with rows as follow
         % 1 - x range of grid in Lat
         % 2 - y range in Long
         % 3 - z range in m
         % 4    cell size X,Y in meters
end

methods
%*******************************************************************
function obj = TBGModel(Filename)
 %   class constructor ... simple
%  obj.LoadFromFile(Filename);
end

function LoadFromNETCDFFile(obj, Filename)
  % read netcdf file
  % otherwise it can be anything like shape fiels too
  % but shape fiel format should be found
  d = ncread(Filename,'dimension');
  obj.Nx = d(1);  obj.Ny = d(2); obj.Nz = 1;l
  
  size(1,:) = ncread(Filename,'x_range');
  size(2,:) = ncread(Filename,'y_range');
  size(3,:) = ncread(Filename,'z_range');
  size(4,:) = ncread(Filename,'spacing'); 
  
  z = ncread(Filename,'z');
  obj.sigma = reshape(z,obj.Nx,obj.Ny);
  
  % compute Lat, Long
  
  
end;

function LoadFromXYZFile(obj, Filename)
  % read netcdf file
  % otherwise it can be anything like shape fiels too
  % but shape fiel format should be found
  data = importdata(Filename);
  %geodata = geoshape(data(:,1), data(:,2), 'Conductivity' , data(:,3));
  
  % UTM projection
  utmstruct = defaultm('utm'); 
  utmstruct.zone = '34V'; 
  %utmstruct.geoid = ellipsoid; 
  
  utmstruct = defaultm(utmstruct);
  [x,y,sigma] = mfwdtran(utmstruct,data(:,1),data(:,2),data(:,3));
   
  obj.Nx = size(x);
  obj.Ny = size(y);
  obj.Nz = 1;
  
  
  
  [xi,yi] = meshgrid(min(x):5000:max(x), min(y):5000:max(y));
  obj.sigma = griddata(x,y,sigma,xi,yi,'natural');
  obj.Nx = size(xi,1);
  obj.Ny = size(xi,2);
  obj.Nz = 1;
  
 % geodata = geoshape(data(:,1), data(:,2));
 % geoshow(geodata);
 % geodata.Geometry = 'line'  
  %geoshow(geodata);
end;


%*******************************************************************

end     % methods
end    % classdef