classdef latlontools
    %   class to define latitude and longitude of grid
    %   and to select a geographic region for plotting
    %   (c) Gary Egbert & Anna Kelbert, 2013-2017
    
    properties(Constant)
        EarthRad = 6378.137; % km
        EarthEcc2 = 0.00669437999014; % eccentricity squared
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = latlontools(~)
            % empty constructor            
        end
    end
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = kmPerDeg(lat0,depth)
            % res = kmPerDeg(lat0,depth)
            % depth in km is optional otherwise use Earth's surface
            a = latlontools.EarthRad;
            e2 = latlontools.EarthEcc2;
            if nargin > 0
                if nargin == 1
                    depth = 0;
                end
                res = pi*(a-depth)*(1-e2)/(180*(1-e2*sin(lat0*pi/180)^2)^(3/2));
            else
                res = 1.1119e+02; % for latitude ~ 48 degrees
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [xy] = ll2xy(ll,lat0,lon0,depth)
            % [xy] = ll2xy(ll,lat0,lon0,depth)
            if nargin < 4
                xy(1,:) = latlontools.kmPerDeg(lat0)*(ll(1,:)-lat0);
                xy(2,:) = latlontools.kmPerDeg(lat0)*(ll(2,:)-lon0).*cos(ll(1,:)*pi/180);
            elseif nargin == 4
                xy(1,:) = latlontools.kmPerDeg(lat0,depth)*(ll(1,:)-lat0);
                xy(2,:) = latlontools.kmPerDeg(lat0,depth)*(ll(2,:)-lon0).*cos(ll(1,:)*pi/180);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ll] = xy2ll(xy,lat0,lon0,depth)
            % [ll] = xy2ll(xy,lat0,lon0,depth)
            if nargin < 4
                ll(1,:) = xy(1,:)/latlontools.kmPerDeg(lat0)+lat0;
                ll(2,:) = (xy(2,:)/latlontools.kmPerDeg(lat0))./cos(ll(1,:)*pi/180)+lon0;
            elseif nargin == 4
                ll(1,:) = xy(1,:)/latlontools.kmPerDeg(lat0,depth)+lat0;
                ll(2,:) = (xy(2,:)/latlontools.kmPerDeg(lat0,depth))./cos(ll(1,:)*pi/180)+lon0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [lat0,lon0] = latlon0(lat,lon,x,y)
            %   given latitude, longitude of a list of sites, and model
            %   grid coordinates X and Y (in m) find best fitting origin in
            %   latitude and longitude space
            if size(lat) ~= size(x); lat = lat'; end
            if size(lon) ~= size(y); lon = lon'; end
            X = [ones(size(x)) x ];
            b = X\lat;
            lat0 = b(1);
            X = [ones(size(y)) y ];
            b = X\lon;
            lon0 = b(1);
            if lon0 > 180
                lon0 = lon0-360;
            end           
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function lon = lon360(lon)
            %   given list of longitudes, convert them to [0,360)
            lon(lon<0) = lon(lon<0)+360;
        end       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function lon = lon180(lon)
            %   given list of longitudes, convert them to (-180,180]
            lon(lon>180) = lon(lon>180)-360;
         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function lims = limits(lat,lon,depth)
            %   given list of latitudes, longitudes, depths make limits
            lims.latmin = min(lat);
            lims.latmax = max(lat);
            lims.lonmin = latlontools.lon180(min(lon));
            lims.lonmax = latlontools.lon180(max(lon));
            lims.depthmin = min(depth);
            lims.depthmax = max(depth);
         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function [deltax,nodex] = delta(x)
             % if x are locations on the grid (lat, lon, depth, x, y, z)
             % that define cell centers...
             % then deltax are the corresponding grid cell sizes - this
             % makes an educated guess at the edges assuming equal spacing
             % in the center
             % optionally also outputs the coordinates of the nodes.
             n = length(x);
             if n == 0; deltax = []; nodex = []; return; end
             deltax = zeros(n,1);
             ctri = floor(n/2);
             d = x(ctri+1)-x(ctri);
             deltax(ctri) = d;
             for i=ctri:-1:2
                 deltax(i-1) = (x(i)-x(i-1))*2-deltax(i);
             end
             for i=ctri:1:n-1
                 deltax(i+1) = (x(i+1)-x(i))*2-deltax(i);
             end
             if nargout > 1
                 if size(x,1)==1 
                     % x is horizontal
                     nodex = [x(1)-deltax(1)/2; x'+deltax/2];
                 else
                     % x is vertical
                     nodex = [x(1)-deltax(1)/2; x+deltax/2];
                 end
             end
         end     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function [deltaz,nodez] = deltaz(z)
             % if z are vertical locations on the grid (depth, z)
             % that define cell centers...
             % similar to delta() but assumes that the top is at zero.
             % optionally also outputs the coordinates of the nodes.
             n = length(z);
             if n == 0; deltaz = []; nodez = []; return; end
             deltaz = zeros(n,1);
             deltaz(1) = z(1)*2;
             for i=2:n
                 deltaz(i) = (z(i) - sum(deltaz(1:i-1)))*2;
             end
             if nargout > 1
                 if size(z,1)==1 
                     % z is horizontal
                     nodez = [0  cumsum(deltaz)];
                 else
                     % z is vertical
                     nodez = [0; cumsum(deltaz)];
                 end
             end
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function [lims,posn] = getLimits(region)
            % Define limits for plotting or output for the usual areas;
            % optionally, also output optimal paper position for printing.
            % Regions: global
            %          subduction
            %          NorthAmerica
            %          WesternUS
            %          SouthAmerica
            %          Argentina
            %          Africa
            %          Europe
            %          EasternEurope
            %          NorthernEurope
            %          Eurasia
            %          Japan
            %          India
            %          WesternAsia
            %          Asia
            %          Pacific
            %          Australia
            %          Canada
            %          Antarctica            
            
            if nargin < 1
                region = 'global';
            end
            
            switch lower(region)
                case 'subduction'
                    lims = struct('lonmin',-280,'lonmax',-20,'latmin',-66,'latmax',66);
                    posn = [1 0.2 12 7];
                case 'northamerica'
                    lims = struct('lonmin',-175,'lonmax',-45,'latmin',10,'latmax',80);
                    posn = [1 0.2 11 9];
                case 'usa&canada'
                    lims = struct('lonmin',-175,'lonmax',-45,'latmin',23,'latmax',80);
                    posn = [1 0.2 11 9];
                case 'usa'
                    lims = struct('lonmin',-150,'lonmax',-40,'latmin',5,'latmax',65);
                    posn = [1 0.2 11 9];
                case 'continentalus'
                    lims = struct('lonmin',-130,'lonmax',-65,'latmin',23,'latmax',51);
                    posn = [1 0.2 16 8];
                case 'continentalus+'
                    lims = struct('lonmin',-135,'lonmax',-60,'latmin',23,'latmax',53);
                    posn = [1 0.2 16 9];
                case 'westernus'
                    lims = struct('lonmin',-130,'lonmax',-100,'latmin',24,'latmax',51);
                    posn = [1 0.2 10 12];
                case 'easternus'
                    lims = struct('lonmin',-104,'lonmax',-60,'latmin',24,'latmax',51);
                    posn = [1 0.2 12 10];
                case 'northwesternus'
                    lims = struct('lonmin',-126,'lonmax',-106,'latmin',36,'latmax',52);
                    posn = [1 0.2 12 10];
                case 'alaska'
                    lims = struct('lonmin',-184,'lonmax',-126,'latmin',46,'latmax',76);
                    posn = [1 0.2 12 10];
                case 'cascadia'
                    lims = struct('lonmin',-126,'lonmax',-114,'latmin',36,'latmax',52);
                    posn = [1 0.2 6 10];
                case 'iowa'
                    lims = struct('lonmin',-97.5,'lonmax',-89.5,'latmin',39,'latmax',45);
                    posn = [1 0.2 10 10];
                case 'florida'
                    lims = struct('lonmin',-91,'lonmax',-74,'latmin',18,'latmax',38);
                    posn = [1 0.2 5 12];
                case 'southamerica'
                    lims = struct('lonmin',-90,'lonmax',-20,'latmin',-70,'latmax',30);
                    posn = [1 0.2 8 12];
                case 'argentina'
                    lims = struct('lonmin',-80,'lonmax',-30,'latmin',-50,'latmax',-10);
                    posn = [1 0.2 10 8];
                case 'africa'
                    lims = struct('lonmin',-30,'lonmax',90,'latmin',-50,'latmax',50);
                    posn = [1 0.2 10 8];
                case 'europe'
                    lims = struct('lonmin',-15,'lonmax',30,'latmin',30,'latmax',75);
                    posn = [1 0.2 10 13];
                case 'easterneurope'
                    lims = struct('lonmin',-30,'lonmax',60,'latmin',30,'latmax',65);
                    posn = [1 0.2 10 13];
                case 'northerneurope'
                    lims = struct('lonmin',-30,'lonmax',60,'latmin',45,'latmax',80);
                    posn = [1 0.2 10 13];
                case 'eurasia'
                    lims = struct('lonmin',-30,'lonmax',180,'latmin',-40,'latmax',65);
                    posn = [1 0.2 12 15];
                case 'japan'
                    lims = struct('lonmin',90,'lonmax',180,'latmin',-20,'latmax',65);
                    posn = [1 0.2 9 8];
                case 'philippine'
                    lims = struct('lonmin',110,'lonmax',170,'latmin',4,'latmax',64);
                    posn = [1 0.2 7.2 7];
                case 'india'
                    lims = struct('lonmin',60,'lonmax',110,'latmin',0,'latmax',40);
                    posn = [1 0.2 10 8];
                case 'westernasia'
                    lims = struct('lonmin',60,'lonmax',150,'latmin',-20,'latmax',60);
                    posn = [1 0.2 4.3 10];
                case 'asia'
                    lims = struct('lonmin',30,'lonmax',180,'latmin',-20,'latmax',60);
                    posn = [1 0.2 14 8];
                case 'pacific'
                    lims = struct('lonmin',86,'lonmax',306,'latmin',-40,'latmax',70);
                    posn = [1 0.2 12 7];
                case 'australia&islands'
                    lims = struct('lonmin',110,'lonmax',180,'latmin',-56,'latmax',6);
                    posn = [1 0.2 4.5 10];
                case 'australia'
                    lims = struct('lonmin',111,'lonmax',156,'latmin',-45,'latmax',-10);
                    posn = [1 0.2 8 8];
                case 'canada'
                    lims = struct('lonmin',-180,'lonmax',-15,'latmin',45,'latmax',80);
                    posn = [1 0.2 12 15];
                case 'antarctica'
                    lims = struct('lonmin',-30,'lonmax',200,'latmin',-80,'latmax',0);
                    posn = [1 0.2 12 15];
                otherwise
                    % by default, use global limits
                    lims = struct('lonmin',0,'lonmax',360,'latmin',-90,'latmax',90);
                    posn = [1 0.2 10 5];
            end
            
        end
         
    end
end