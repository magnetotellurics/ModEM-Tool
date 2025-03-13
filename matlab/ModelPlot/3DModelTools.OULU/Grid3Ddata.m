classdef Grid3Ddata < TGrid3D
% Railo 2014
%
%

properties
    postSigma,likeSigma,priSigma
    postVar,likeVar,priVar
    D
end

methods
    
%*******************************************************************
% Class constructor
function obj = Grid3Ddata(SMAP,Data2d,Grid3d,varargin)
    obj.D = ones(1,1);
% Combines SMAP "prior" < ListData and 2D data < T2DModel
% ListData from SMAP can be created by the function TBGSModel.makeListData
    if nargin > 3
        a = varargin{1}; % "uncertainty" or variance at the measurement
        b = varargin{2}; % "smoothing" factor
        c = varargin{3}; % prior variance
    else
        a = 1;
        b = 1;
        c = 10;
    end
    xShift = Grid3d.origin.x; % position of the 1st x coord in the UTM
    yShift = Grid3d.origin.y; % position of the 1st y coord in the UTM
    flag1 = 0; flag2 = 0;
    h = waitbar(0,'Initializing waitbar...');
    for i = 1:Grid3d.Nx
        waitbar(i/Grid3d.Nx,h,'In progress...')
        xUTM = sum(Grid3d.Dx(1:i-1)) + Grid3d.Dx(i)/2 + xShift;
        for j = 1:Grid3d.Ny
            yUTM = sum(Grid3d.Dy(1:j-1)) + Grid3d.Dy(j)/2 + yShift;
            % now [xUTM yUTM] gives the centre of cell
            % finds possible XY coordinates from 2D data
             indXY = Data2d.findXY([xUTM yUTM]);
            for k = 1:Grid3d.Nz
                zDepth = sum(Grid3d.Dz(1:k-1)) + Grid3d.Dz(k)/2;
                
                % find corresponding coordinate from 2D data
                indZ = Data2d.findZ(indXY,zDepth);
                d = norm([Data2d.xUTM(indZ)-xUTM,...
                    Data2d.yUTM(indZ)-yUTM, Data2d.Z(indZ)-zDepth],2);
                obj.D(i,j,k) = d;
                % find corresponding coordinate from SMAP data (prior)
                d2 = [SMAP.xUTM-xUTM,SMAP.yUTM-yUTM,SMAP.Z-zDepth];
                temp = sum(d2'.^2);
                indminSMAP = find(temp==min(temp));
                if length(indminSMAP) > 1 % checks possible errors
                	flag2 = 1;
                end
                indSMAP = indminSMAP(1);
                
                % combine the data and prior
                dataVar = a*exp(b(k)*d); % variance based on distance
                priorVar = c;
                if dataVar/priorVar > 1e3
                    psMean = SMAP.sigma(indSMAP);
                    psVar = priorVar;
                elseif priorVar/dataVar > 1e3
                    psMean = Data2d.sigma(indZ);
                    psVar = dataVar;
                else
                    [psMean,psVar] = GausPDFProd(Data2d.sigma(indZ),...
                        dataVar, SMAP.sigma(indSMAP), priorVar);
                end
                % saves values
                obj.postSigma(i,j,k) = psMean;
                obj.likeSigma(i,j,k) = Data2d.sigma(indZ);
                obj.priSigma(i,j,k) = SMAP.sigma(indSMAP);
             
                obj.postVar(i,j,k) = psVar;
                obj.likeVar(i,j,k) = dataVar;
                obj.priVar(i,j,k) = priorVar;

% %               OLD VERSION
%                 dis = [Data2d.xUTM-xUTM,Data2d.yUTM-yUTM,Data2d.Z-zDepth];
%                 temp = sum(dis'.^2);
%                 indmin2d = find(temp==min(temp));
%                 d = sqrt(min(temp));
%                 if length(indmin2d) > 1
%                     flag1 = 1;
%                 end
%                 ind1 = indmin2d(1);
%                 obj.likeSigma(i,j,k) = Data2d.sigma(ind1);
            end
        end
    end
    close(h);
    obj.Dx = Grid3d.Dx; obj.Dy = Grid3d.Dy; obj.Dz = Grid3d.Dz;
    obj.Nza = Grid3d.Nza; obj.origin = Grid3d.origin;
    obj.Nx = Grid3d.Nx; obj.Ny = Grid3d.Ny; obj.Nz = Grid3d.Nz;
    if flag1 == 1
        disp('Warning: the nearest 2d data point(s) was not unique!')
    end
    if flag2 == 1
        disp('Warning: the nearest SMAP point(s) was not unique!')
    end
end
%*******************************************************************


function PlotPlainView(obj, Nz)
  %  plot plain view of sigma for Nz depth
  if Nz < obj.Nz+1
    figure;
    contourf(1:obj.Nx, 1:obj.Ny,  squeeze(obj.sigma(:,:,Nz))');
  else
    disp('Specify another depth layer to plot. Nz is out of range');  
  end
  
end

% Plot functions



end     % methods
end    % classdef