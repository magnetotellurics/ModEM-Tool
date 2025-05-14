%% testRoutine V0.1
%% VERSION 0.1 README (testRoutine, testRoutine2): ROUTINES FOR TESTING DATA+SMAP with 3D GRID
% FEATURES
% - Freely defined "flat" 3-D grid i.e. (x,y,z)-grid
% - Possibility to read the SMAP data for the grid (only till 60 km depth)
% - Possibility to read a 2-D data profile for the grid; strike angle can
%   be taken in count when data is combined
% - Possibility to combine these datasets 
%   (with so called Gaussian distribution approach)
% - Can be use inside only one UTM-zone (which can be freely chosen)
% - Handles the data and results with UTM-coordinates
%
% USER NOTES
% - testRoutine code:
%       PART 1. defines the origin (x,y,UTM-zone) and grid
%       PART 2. reads the SMAP data
%       PART 3. reads a 2-D data
%       PART 4. combines the datasets (once)
% - testRoutine2 code:
%       PART 4. combines the datasets (if objects in PARTS 1-3 are created)
%       PART 5. plots the data if the datasets are combined
% - If you test different parameters, run testRoutine2! (saves time)
% - More details on the comments
%
% NEXT VERSION(S) FEATURES (LITTLE-BY-LITTLE)
% 0.    comment all the object class and tip-function codes & make
%       the codes clearer
% 1.    possibility to read/combine any number of 1/2/3D datasets
% 2.    uses weight-functions rather than "Gaussian distributions" in data
%       combinations; easier to define and both lead similar results (it's
%       more and less only a way of thinking)
%       For 1-2 check out "testRoutine3" and class "Grid3DdataTEST"
% 3.    data combination in geographical coordinates; which removes
%       the current UTM-zone restriction (or an other approach to kill the
%       problem)
% 4.    grid possible to define both in geographical and UTM-coords
% 5.    Include more depths to the background model (new layers)
% 6.    Include topography data to the background model (new layers)
%
% By Jesse Railo, jesse.railo@oulu.fi, 9.9.2014

% clear old works
clear all

%% PARTS 1: define the grid and zone
% define the UTM-zone we work in
zone = '33 V';

% define the origin of grid
origin.x = 4.05e5; origin.y = 6.96e6; origin.zone = zone;
% define the grid dimensions
Dx = 3500*ones(1,60); Dy = 2200*ones(1,60); Dz = 1.21.^(1:1:15)*1000;
% number of "air" layers?
Nza = 0;
% create the object GRID3D of TGrid3D
GRID3D = TGrid3D(Dx,Dy,Dz,Nza,origin);
% NOTE: This is the final grid user wants to use

%% PART 2: create the object SMAP of TBGSModel based on the SMAP file
% variable '1' should be path for the SMAP file
SMAP = TBGSModel(1);
% creates the object SMAPList of ListData; conductance -> conductivity etc.
SMAPList = SMAP.makeListData;
% makes the SMAPList shorter (we are only interested one zone)
SMAPList.restrictZone(zone);

%% PART 3: create the object DATA2D of T2DModel based on the 2-D file
% variable '1' should be path for the 2-D file (similar as above)
DATA2D = T2DModel(1);
% see above
DATA2D.restrictZone(zone);

%% PART 4
% define parameters of data combination and combine data

% "background model variance/uncertainty"
c = 0.1;
% "measurement variance/uncertainty"
a = 0.02;
% combination is made differently in different depths
Z = cumsum(Dz)-Dz/2;
% "reaching factor" for each depth 
% (smaller -> 2-D data reaches further) (1e-3...2e-5)
% now kind of crazy function; could be proportional to Z^(-1)
b = log(Z).^2.5/max(log(Z).^2.4).*Z.^(-1);

% creates the object RES of Grid3Ddata
% does the trick...
RES = Grid3Ddata(SMAPList,DATA2D,GRID3D,a,b,c);
% Has e.g.
%     postSigma: [60x60x15 double]
%     likeSigma: [60x60x15 double]
%      priSigma: [60x60x15 double]
%       postVar: [60x60x15 double]
%       likeVar: [60x60x15 double]
%        priVar: [60x60x15 double]
%             D: [60x60x15 double]
%            Nx: 60
%            Ny: 60
%            Nz: 15
%           Nza: 0
%            Dx: [1x60 double]
%            Dy: [1x60 double]
%            Dz: [1x15 double]
%        origin: [1x1 struct]