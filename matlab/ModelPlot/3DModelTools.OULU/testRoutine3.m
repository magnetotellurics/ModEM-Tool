% ROUTINE FOR TESTING DATA+SMAP with 3D GRID
%% PART 1
clear all
zone = '33 V';

SMAP = TBGSModel(1);
SMAPList = SMAP.makeListData;
SMAPList.restrictZone(zone);

DATA2D = T2DModel(1);
DATA2D.restrictZone(zone);

Dx = 3500*ones(1,60); Dy = 2200*ones(1,60); Dz = 1.21.^(1:1:15)*1000;
Nza = 0;

origin.x = 4.05e5; origin.y = 6.96e6; origin.zone = zone;
GRID3D = TGrid3D(Dx,Dy,Dz,Nza,origin);

%% PART 2

a = 1; % measurement variance/uncertainty
Z = cumsum(Dz)-Dz/2;
b = log(Z).^2.5/max(log(Z).^2.4).*Z.^(-1); % smoothing factor (smaller -> smoother) (1e-4)
c = 0.1;
RES = Grid3DdataTEST(GRID3D,SMAPList,ones(0,0),DATA2D,ones(0,0),a,b,c);

% axis grid for the plot
X = cumsum(RES.Dx)-Dx/2 + RES.origin.x;
Y = cumsum(RES.Dy)-Dy/2 + RES.origin.y;
Z = cumsum(RES.Dz)-Dz/2;

cmin = min(min(min(RES.postMean)));
cmax = max(max(max(RES.postMean)));

[A,B] = DATA2D.giveLine;

for i = 1:length(RES.Dz)
    figure(1);
	imagesc(X,Y,RES.postMean(:,:,i)'); c = colorbar;
    set(gca,'YDir','normal');
    ylabel(c,'log(conductivity)')
    caxis([cmin cmax])
	axis equal
	M(i) = getframe;
    zone = RES.origin.zone;
    title(['UTM-zone: ' zone '   Depth: ' num2str(Z(i)) ' +- cellsize m']);
    xlabel('m');
    ylabel('m');
    hold on
    plot([A(1) B(1)],[A(2) B(2)])
    pause
end