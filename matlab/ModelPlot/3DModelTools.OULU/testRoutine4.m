% ROUTINE FOR TESTING DATA+SMAP with 3D GRID (LAT,LONG,Z)
%% PART 1: Grid, SMAP, Data etc.
clear all
zone = '33 V';

SMAP = TBGSModel(1);
SMAPList = SMAP.makeListData;
%SMAPList.restrictZone(zone);

DATA2D = T2DModel(1);
%DATA2D.restrictZone(zone);

%Dx = 3500*ones(1,60); Dy = 2200*ones(1,60); Dz = 1.21.^(1:1:15)*1000;
Dx = 5/60*ones(1,60);
Dy = 10/60*ones(1,60);
Dz = 1.21.^(1:1:15)*1000;
Nza = 0;

%origin.x = 4.05e5; origin.y = 6.96e6; origin.zone = zone;
origin.x = 65;
origin.y = 15;
origin.zone = zone;
GRID3D = TGrid3D(Dx,Dy,Dz,Nza,origin);

%% PART 2: data combination

a = 1   ; % measurement variance/uncertainty
Z = cumsum(Dz)-Dz/2;
b = log(Z).^2.5/max(log(Z).^2.4).*Z.^(-1); % smoothing factor (smaller -> smoother) (1e-4)
c = 0.1;
RES = Grid3DdataTEST2(GRID3D,SMAPList,ones(0,0),DATA2D,ones(0,0),a,b,c);

%% PART 3: plotting
% axis grid for the plot
X = cumsum(RES.Dx)-Dx/2 + RES.origin.x;
Y = cumsum(RES.Dy)-Dy/2 + RES.origin.y;
Z = cumsum(RES.Dz)-Dz/2;

DATA2Dt = RES.DATA2D{1};

cmin = min(min(min(RES.postMean)));
cmin = min(cmin,min(min(min(RES.BGMODEL))));
cmin = min(cmin,min(min(min(DATA2Dt))));

cmax = max(max(max(RES.postMean)));
cmax = max(cmax,max(max(max(RES.BGMODEL))));
cmax = max(cmax,max(max(max(DATA2Dt))));

% assume there is only one line of sights
[A,B] = DATA2D(1).giveLineLatLong;
cmin = -9;
cmax = 0;
for i = 1:length(RES.Dz)
    figure(1);
    caxis([cmin cmax])
	imagesc(Y,X,RES.postMean(:,:,i)); col = colorbar;
    set(gca,'YDir','normal');
    ylabel(col,'log(conductivity)')
    caxis([cmin cmax])
	axis equal
	M(i) = getframe;
    zone = RES.origin.zone;
    title(['Near UTM-zone: ' zone '   Depth: ' num2str(Z(i)) ' +- cellsize m']);
    xlabel('Longitude');
    ylabel('Latitude');
    hold on
    plot([A(2) B(2)],[A(1) B(1)])
    
    figure(2);
    caxis([cmin cmax])
	imagesc(Y,X,RES.BGMODEL(:,:,i)); col = colorbar;
    set(gca,'YDir','normal');
    ylabel(col,'log(conductivity)')
    caxis([cmin cmax])
	axis equal
	M(i) = getframe;
    zone = RES.origin.zone;
    title(['Near UTM-zone: ' zone '   Depth: ' num2str(Z(i)) ' +- cellsize m']);
    xlabel('Longitude');
    ylabel('Latitude');
    hold on
    plot([A(2) B(2)],[A(1) B(1)])
    
    figure(3);
    caxis([cmin cmax])
	imagesc(Y,X,DATA2Dt(:,:,i)); col = colorbar;
    set(gca,'YDir','normal');
    ylabel(col,'log(conductivity)')
    caxis([cmin cmax])
	axis equal
	M(i) = getframe;
    zone = RES.origin.zone;
    title(['Near UTM-zone: ' zone '   Depth: ' num2str(Z(i)) ' +- cellsize m']);
    xlabel('Longitude');
    ylabel('Latitude');
    hold on
    plot([A(2) B(2)],[A(1) B(1)])
    pause
end