classdef TBGSModel < handle
% Railo 2014
%

properties
 S                  % 3D conductance map
 Lat, Long      
 water_S, water_h % water conductance, water thickness
 sedim_S, sedim_h % sediment conductance, sediment thickness
 crust0, crust1, crust2, crust3, crust4, crust5
 % crust from 0..10km, 10..20km,...,50..60km
 grLat, grLong, grZ, Nx, Ny, Nz
 N              % length of the list
end

methods
%*******************************************************************
function obj = TBGSModel(Filename)
 %   class constructor ... simple
    obj.LoadFromFile(Filename);
end

function LoadFromFile(obj, Filename)
    Filename = 'data/Smap_all.xyz'; % for testing, should be removed
    fid = fopen(Filename);
    a = textscan(fid, '%f64%f64%f64%f64%f64%f64%f64%f64%f64%f64%f64%f64','headerlines',20);
    fclose(fid);
    a = cell2mat(a);
    
    % save all the data
    obj.Lat = a(:,1);
    obj.Long = a(:,2);
    obj.water_S = a(:,3);
    obj.water_h = a(:,4);
    obj.sedim_S = a(:,5);
    obj.sedim_h = a(:,6);
    obj.crust0 = a(:,7);
    obj.crust1 = a(:,8);
    obj.crust2 = a(:,9);
    obj.crust3 = a(:,10);
    obj.crust4 = a(:,11);
    obj.crust5 = a(:,12);
    obj.N = length(obj.Lat);
end;

function TBGSlist = makeListData(obj)
    TBGSlist = ListData(); % no file needed
    temp = [obj.crust0 obj.crust1 obj.crust2 ...
        obj.crust3 obj.crust4 obj.crust5];
    counter = 0;
    for i = 1:obj.N
        for j = 1:8
            index = j + 8*(i-1) - counter;
            flag = 0;
            if j == 1
                % depth of 1st layer is water_h
                if obj.water_h(i) ~= 0
                    TBGSlist.Z(index) = obj.water_h(i)/2;
                    TBGSlist.sigma(index) = obj.water_S(i)/obj.water_h(i);
                else
                    counter = counter + 1;
                    flag = 1;
                end
            elseif j == 2
                % depth of 2nd layer is water_h + sedim_h
                if obj.sedim_h(i) ~= 0
                    TBGSlist.Z(index) = obj.water_h(i) + obj.sedim_h(i)/2;
                    TBGSlist.sigma(index) = obj.sedim_S(i)/obj.sedim_h(i);
                else
                    counter = counter + 1;
                    flag = 1;
                end
            else  
                % depth of next layers is constant 10 km
                TBGSlist.Z(index) = obj.water_h(i) + obj.sedim_h(i) +...
                (j-3)*10000 + 5000;
                TBGSlist.sigma(index) = temp(i,(j-2))/10000;
            end
            if flag == 0
                TBGSlist.Lat(index) = obj.Lat(i);
                TBGSlist.Long(index) = obj.Long(i);
            end
        end
    end
    TBGSlist.Lat = TBGSlist.Lat';
    TBGSlist.Long = TBGSlist.Long';
    TBGSlist.Z = TBGSlist.Z';
    TBGSlist.sigma = log(TBGSlist.sigma');
    TBGSlist.makeUTMList;
    TBGSlist.N = index;
end

function make3dGrid(obj)
    % make 3D grid for the data
    obj.grLat = unique(a(:,1));
    obj.grLong = unique(a(:,2));
    % Z SHOULD BE THINKING AGAIN!!!
    obj.grZ = ones(length(a(:,1)),8)*10; % 10 km for crust
    obj.grZ(:,1) = a(:,4); % water
    obj.grZ(:,2) = a(:,6); % sediment
    obj.Nx = length(obj.grLat); % dimension of lat
    obj.Ny = length(obj.grLong); % dimension of long
    obj.Nz = 8; % water, sediment, 6 layers of crust
    Sval = [a(:,3) a(:,5) a(:,7) a(:,8) a(:,9) a(:,10) a(:,11) a(:,12)];
    % conductance of water, sediment, crust 0...60 km with 10 km step
    obj.S = NaN(obj.Nx,obj.Ny,obj.Nz); % 3D conductance map
    for i = 1:obj.Nx
        indX = find(a(:,1)==obj.Lat(i));
        for j = 1:obj.Ny
            indY = find(a(indX,2)==obj.Long(j));
            indY = indX(indY);
            for k = 1:obj.Nz
                if isempty(indY)
                    obj.S(i,j,k) = NaN;
                else
                    obj.S(i,j,k) = Sval(indY,k);
                end
            end
        end
    end
end
%*******************************************************************

end     % methods
end    % classdef