Format = 'ASCII';
Mplot = 1;
%   Driver script for interactive impedance plotting program
%   Root directory to start browsing for impedance files
fdir = '.';
%  prompt user for impedance file name
%  note: if we use suffixes to name impedance files (e.g. *.grd)
%   we can change following line to
filt = [fdir '/*.dat; *.res; *.imp'];
[filename, pathname] = uigetfile(filt, 'Impedance File');
Zfile = [pathname filename];
if findstr('.imp',Zfile) % imp file extension indicates old ascii format
    Format = 'ASCII OLD';
end
switch Format
  case 'OLD KUSH'
      [Imp,loc,T] = rdZall(Zfile);
  case 'NEW KUSH'
      [Imp,loc,T] = rdZ(Zfile);
  case 'ASCII OLD'
      [AllData] = readZ_3D_old(Zfile,'[mV/km]/[nT]');
      [Imp,loc,T,sites,origin] = allData2Imp(AllData,Mplot);
    otherwise
      [AllData] = readZ_3D(Zfile,'[mV/km]/[nT]');
      [Imp,loc,T,sites,origin] = allData2Imp(AllData,Mplot);
end

if isempty(origin)
    Mplot = 0;
end

%[I,J,X,Y] = sortLoc(loc);
SiteX = loc(1,:);
SiteY = loc(2,:);
DiffX = (max(SiteX)-min(SiteX))/100;
DiffY = (max(SiteY)-min(SiteY))/100;
X = min(SiteX)-DiffX/2:DiffX:max(SiteX)+DiffX/2;
Y = min(SiteY)-DiffY/2:DiffY:max(SiteY)+DiffY/2;
nX = length(X);
nY = length(Y);
[nTF,nmode,nSite,nPer] = size(Imp);
Z = zeros(nTF,nmode,nX,nY,nPer);
Z = Z+1i*Z;
[gdY,gdX] = meshgrid(Y,X);
for iPer = 1:nPer
    for iTF = 1:nTF
        for imode = 1:nmode
            gdImp = griddata(SiteY,SiteX,Imp(iTF,imode,:,iPer),gdY,gdX);
            Z(iTF,imode,:,:,iPer) = gdImp;
        end
    end
end
% for k = 1:nX
%    ind = find(I == k);
%    Z(:,:,k,J(ind),:) = Imp(:,:,ind,:);
% end

options.T = T;
options.X = X;
options.Y = Y;
options.Xlim = [min(X)-5*DiffX max(X)+5*DiffX];
options.Ylim = [min(Y)-5*DiffY max(Y)+5*DiffY];
options.nPer = length(T);
options.nX = nX;
options.nY = nY;
options.SiteX = SiteX;
options.SiteY = SiteY;
options.slice = 'T';
options.Np = 1;
options.Mode = 2;
options.Comp = 1;
options.Mplot = Mplot;
options.PlotSites = 0;
options.File = Zfile;

options.iXlim(1) = 1;
options.iXlim(2) = nX;
options.iYlim(1) = 1;
options.iYlim(2) = nY;
options.iZlim(1) = 1;
options.iZlim(2) = nPer;
ImpPlotSet(Z,options);
