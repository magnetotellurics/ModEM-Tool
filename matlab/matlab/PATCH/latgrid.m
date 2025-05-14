classdef latgrid
    %   class to define latitude and longitude of grid
    properties
        grid
        LAT
        LON
        lat0
        lon0
        x
        y
        kmPerDeg  = 1.1119e+02;
    end
    
    methods
        function [obj] = latgrid(grid,lat0,lon0)
            %   class constructor
            obj.grid = grid;
            obj.lat0 = lat0;
            obj.lon0  = lon0;
            obj.x = [0 ; cumsum(grid.dx)]+grid.origin(1);
            obj.x = (obj.x(1:end-1)+obj.x(2:end))/2;
            obj.y = [0 ; cumsum(grid.dy)]+grid.origin(2);
            obj.y = (obj.y(1:end-1)+obj.y(2:end))/2;
            obj.LAT = (obj.x/obj.kmPerDeg)*ones(1,grid.Ny) + lat0;
            
            obj.LON = (ones(grid.Nx,1)*obj.y'/obj.kmPerDeg)./cos(obj.LAT*pi/180) + lon0;
        end
        function result = limits(obj)
            %   find lat, lon, depth limits
            result = struct('latmin',min(min(obj.lat)),'latmax',max(max(obj.lat)),...
                'lonmin',min(min(obj.lon)),'lonmax',max(max(obj.lon)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [xy] = ll2xy(obj,ll)
            xy(1,:) = obj.kmPerDeg*(ll(1,:)-obj.lat0);
            xy(2,:) = obj.kmPerDeg*(ll(2,:)-obj.lon0).*cos(ll(1,:)*pi/180);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ll] = xy2ll(obj,xy)
            ll(1,:) = xy(1,:)/obj.kmPerDeg+obj.lat0;
            ll(2,:) = (xy(2,:)/obj.kmPerDeg)./cos(ll(1,:)*pi/180)+obj.lon0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static)
        function [lat0,lon0] = origin(lat,lon,x,y)
            %   given latitude, longitude of a list of sites, and model
            %   grid coordinates X and Y (in m) find best fitting origin in
            %   latitude and longitude space
            X = [ones(size(x)) x ];
            b = X\lat';
            lat0 = b(1);
            X = [ones(size(y)) y ];
            b = X\lon';
            lon0 = b(1);
            if lon0 < 0
                lon0 = 360+lon0;
            end           
        end      
    end
end