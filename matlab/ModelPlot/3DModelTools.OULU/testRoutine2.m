%% testRoutine2 V0.1
% If you encounter difficulties, read the README of testRoutine file
%
% By Jesse Railo, jesse.railo@oulu.fi, 9.9.2014

%% PART 4:	define parameters of data combination and combine data

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

%% PART 5:	plotting results

% finds limits for axes
% "post" refers to the results
PostSigma = RES.postSigma;
cmin = min(min(min(PostSigma)));
cmax = max(max(max(PostSigma)));
% "pri" refers to the background model
PriSigma = RES.priSigma;
pmin = min(min(min(PriSigma)));
pmax = max(max(max(PriSigma)));
% "like" refers to the 2-D data
LikeSigma = RES.likeSigma;
lmin = min(min(min(LikeSigma)));
lmax = max(max(max(LikeSigma)));

mini = min(cmin,min(pmin,lmin));
maxi = max(cmax,max(pmax,lmax));

lmin = mini; cmin = mini; pmin = mini;
lmax = maxi; cmax = maxi; pmax = maxi;

% axis grid for the plot
X = cumsum(RES.Dx)-Dx/2 + RES.origin.x;
Y = cumsum(RES.Dy)-Dy/2 + RES.origin.y;
Z = cumsum(RES.Dz)-Dz/2;

% solves the end points of measurement lines
[A,B] = DATA2D.giveLine;

% plots every layer separately
% pause after each loop
for i = 1:length(RES.Dz)
    figure(1);
    % plots the colormap
	imagesc(X,Y,PostSigma(:,:,i)'); c = colorbar;
    % normally x and y axis are reversed...
    set(gca,'YDir','normal');
    % names the colorbar
    ylabel(c,'log(conductivity)')
    % set limits of the axes
    caxis([cmin cmax])
	axis equal
    % saves "frame" in the case you want to play "movie" 
    % (see bottom of the code)
	M(i) = getframe;
    % reads used UTM-zone
    zone = RES.origin.zone;
    % set titel and x, y labes
    title(['UTM-zone: ' zone '   Depth: ' num2str(Z(i)) ' +- cellsize m']);
    xlabel('m');
    ylabel('m');
    hold on
    % plots the line of measurement
    plot([A(1) B(1)],[A(2) B(2)])
    
    % similar as above
    figure(2);
	imagesc(X,Y,PriSigma(:,:,i)'); c = colorbar;
    set(gca,'YDir','normal');
    ylabel(c,'log(conductivity)')
    caxis([pmin pmax])
	axis equal
	M(i) = getframe;
    zone = RES.origin.zone;
    title(['UTM-zone: ' zone '   Depth: ' num2str(Z(i)) ' +- cellsize m']);
    xlabel('m');
    ylabel('m');
    hold on
    plot([A(1) B(1)],[A(2) B(2)])
    
    % similar as above
    figure(3);
	imagesc(X,Y,LikeSigma(:,:,i)'); c = colorbar;
    set(gca,'YDir','normal');
    ylabel(c,'log(conductivity)')
    caxis([lmin lmax])
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
% in the case you want to play the saved frames, see 'help movie'
% movie(M,5,2)