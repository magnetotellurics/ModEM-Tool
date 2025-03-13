%% Class constructor
function obj = Grid3DdataTEST(GRID3D, BGMODEL, DATA1D, DATA2D, DATA3D, varargin)
% INPUTS
%   GRID3D      - the 3d grid of form TGrid3D
%   BGMODEL     - the background model of form ListData
%   DATA1D      - 1d datasets: vector of T1DModels
%   DATA2D      - 2d datasets: vector of T2DModels
%   DATA3D      - 3d datasets: vector of T3DModels
%   varargin    - parameters for data combination
%   
%
% NO OUTPUTS / class constructor
%
% Combines BGMODEL "prior" < ListData and 2D data < T2DModel
% ListData from BGMODEL can be created by the function TBGSModel.makeListData
    if nargin > 3
        a = varargin{1}; % "uncertainty" or variance at the measurement
        b = varargin{2}; % "smoothing" factor
        c = varargin{3}; % prior variance
    else
        a = 1;
        b = 1;
        c = 1/10;
    end
    xShift = GRID3D.origin.x; % position of the 1st x coord in the UTM
    yShift = GRID3D.origin.y; % position of the 1st y coord in the UTM
    flag1 = 0; flag2 = 0;
    h = waitbar(0,'Initializing waitbar...');
    
    L1D = length(DATA1D);
    L2D = length(DATA2D);
    L3D = length(DATA3D);
    
    for i = 1:GRID3D.Nx
        waitbar(i/GRID3D.Nx,h,'In progress...')
        xUTM = sum(GRID3D.Dx(1:i-1)) + GRID3D.Dx(i)/2 + xShift;
        for j = 1:GRID3D.Ny
            yUTM = sum(GRID3D.Dy(1:j-1)) + GRID3D.Dy(j)/2 + yShift;
            % now [xUTM yUTM] gives the centre of cell

            % finds possible XY coordinates from 2D data
            for l = 1:L2D
                temp = DATA2D(l);
                indXY{l} = temp.findXY([xUTM yUTM]);
            end
            
            for k = 1:GRID3D.Nz
                zDepth = sum(GRID3D.Dz(1:k-1)) + GRID3D.Dz(k)/2;
                
                % find corresponding distance from 1D datasets
                d1 = ones(1,L1D); sigma1 = ones(1,L1D); % initialize values
                for l = 1:L1D
                    temp = DATA1D(l);
                    temp2 = [temp.xUTM-xUTM,temp.yUTM-yUTM,temp.Z-zDepth];
                    ss = sum(temp2.^2);
                    indmin = find(ss==min(ss));
                    if length(indminBG) > 1 % checks possible errors
                        disp('The nearest 1D point was not unique!');
                    end
                    d1(l) = norm(d(indmin));
                    sigma1(l) = temp.sigma(indmin);
                end
                
                % find corresponding coords and distances from 2D datasets             
                d2 = ones(1,L2D); sigma2 = ones(1,L2D); % initialize values
                for l = 1:L2D
                    temp = DATA2D(l);
                    indZ = temp.findZ(indXY{l},zDepth);
                    d2(l) = norm([temp.xUTM(indZ)-xUTM,...
                    temp.yUTM(indZ)-yUTM, temp.Z(indZ)-zDepth],2);
                    sigma2(l) = temp.sigma(indZ);
                end
                
                % find corresponding coords and distances from 3D datasets
%                 L = 0; an another way to combine data
                d3 = ones(1,L3D); sigma3 = ones(1,L3D); % initialize values
                for l = 1:L3D
                    temp = DATA3D(l);
                    temp2 = [temp.xUTM-xUTM,temp.yUTM-yUTM,temp.Z-zDepth];
                    ss = sum(temp2.^2);
                    indmin = find(ss==min(ss));
                    if length(indminBG) > 1 % checks possible errors
                        disp('The nearest 3D point was not unique!');
                    end
                    d3(l) = norm(d(indmin));
                    sigma3(l) = temp.sigma(indmin);                   
%                     N = temp.N;
%                     d3(L+1:L+N) = norm([temp.xUTM-xUTM,temp.yUTM-yUTM,...
%                         temp.Z-zDepth],2);
%                     sigma3(L+1:L+N) = temp.sigma;
%                     L = L + N;
                end
                
                % find corresponding coordinate from BGMODEL data (prior)
                dBG = [BGMODEL.xUTM-xUTM,BGMODEL.yUTM-yUTM,BGMODEL.Z-zDepth];
                temp = sum(dBG'.^2);
                indminBG = find(temp==min(temp));
                if length(indminBG) > 1 % checks possible errors
                	flag2 = 1;
                end
                indBG = indminBG(1);
                sigmaBG = BGMODEL.sigma(indBG);
                
                % combine the datasets and bgmodel in a vector
                % Order: [1Ds...2Ds...3Ds...BG]
                L1 = length(d1); L2 = length(d2); L3 = length(d3);
                L = L1 + L2 + L3;
                values = ones(1,L); weights = ones(1,L); % initialize
                for l = 1:L1
                    values(1:L1) = sigma1;
                    weights(1:L1) = a*exp(-d1*b(k))/5;
                end
                for l = 1:L2
                    values(L1+1:L1+L2) = sigma2;
                    weights(L1+1:L1+L2) = a*exp(-d2*b(k));
                end
                for l = 1:L3
                    values(L1+L2+1:L) = sigma3;
                    weights(L1+L2+1:L) = 5*a*exp(-d2*b(k));
                end
                values(L+1) = sigmaBG;
                weights(L+1) = c; 
                
                % calculate the value for a grid point
                mean = wMean(values,weights);
                
                % save the results
                obj.postMean(i,j,k) = mean;
            end
        end
    end
    close(h);
    obj.Dx = GRID3D.Dx; obj.Dy = GRID3D.Dy; obj.Dz = GRID3D.Dz;
    obj.Nza = GRID3D.Nza; obj.origin = GRID3D.origin;
    obj.Nx = GRID3D.Nx; obj.Ny = GRID3D.Ny; obj.Nz = GRID3D.Nz;
    if flag1 == 1
        disp('Warning: the nearest 2d data point(s) was not unique!')
    end
    if flag2 == 1
        disp('Warning: the nearest BGMODEL point(s) was not unique!')
    end
end