function [newlat,newlon,z,cov] = mkGrid3D(minlon,maxlon,minlat,maxlat,top,bottom,nlon,nlat,nz,padding,logdz)

% [newlon,newlat,z,cov] = mkGrid3D(minlon,maxlon,minlat,maxlat,top,bottom,nlon,nlat,nz,padding,logdz)
%
% longitudes & latitudes are in degrees, vertical spacing in km.
% optional padding is a string that contains letters W,E,N,S,U,D
% optional logdz, if specified, makes depth distribution logarithmic;
% (top - bottom) = sum_{k=0}^nz (startdz * logdz^k)
%                = startdz * (logdz^(nz+1) - 1)/(logdz - 1)
% therefore startdz = (top - bottom)*(logdz - 1)/(logdz^(nz+1) - 1)
% and n = 0:nz-1; z = cumsum(startdz*(logdz.^n)).

% [elev lon lat]=m_etopo2([-130 -118 40 50]);
% minlon = -126;
% maxlon = -122;
% minlat = 42;
% maxlat = 48;

avglon = (minlon+maxlon)/2;
avglat = (minlat+maxlat)/2;

% nlon = 60;
% nlat = 140;
% nz = 15;
% top = 3;
% bottom = -3; % km

if nargin <= 9
    padding = '';
end

padWest = 0;
padEast = 0;
padSouth = 0;
padNorth = 0;
padUp = 0;
padDown = 0;
increaseFactor = 1.8;
nPad = 8;
pad = cumsum(exp(log(increaseFactor)*[1:nPad]));

EWfactor = (maxlon-minlon)/nlon;
NSfactor = (maxlat-minlat)/nlat;
UDfactor = (top-bottom)/nz;

if findstr(padding,'W') 
    padWest = EWfactor*pad(end:-1:1); % degrees
end
if findstr(padding,'E') 
    padEast = EWfactor*pad; % degrees
end
if findstr(padding,'S')
    padSouth = NSfactor*pad(end:-1:1); % degrees
end
if findstr(padding,'N')
    padNorth = NSfactor*pad ; % degrees
end
if findstr(padding,'U')
    padUp = [32 16 8 4 2]; % km
end
if findstr(padding,'D')
    padDown = UDfactor*[5 10 20 40 80]; % km
end


[elev lon lat]=m_etopo2([minlon-max(padWest) maxlon+max(padEast) ...
    minlat-max(padSouth) maxlat+max(padNorth)]);

% pcolor(lon,lat,elev); shading interp; colorbar; caxis([-1000 1000])

[dx,temp,newlat] = m_lldist([avglon avglon],[minlat maxlat],nlat); dx = dx/nlat;
[dy,newlon,temp] = m_lldist([minlon maxlon],[avglat avglat],nlon); dy = dy/nlon;

if nargin <= 10
    % logdz = 1;
    dz(1:nz) = (top-bottom)/nz;
    disp(['Spacing: lon ' num2str(dy) ' km, lat ' num2str(dx) ...
        ' km, vertical ' num2str(dz(1)) ' km']);
else
    startdz = (top - bottom)*(logdz - 1)/(logdz^nz - 1);
    n = 0:nz-1; dz = startdz*(logdz.^n);
    disp(['Spacing: lon ' num2str(dy) ' km, lat ' num2str(dx) ...
        ' km, vertical starts with ' num2str(dz(1)) ' km']);    
end


if findstr(padding,'W')
    newlon = [minlon-padWest newlon];
end
if findstr(padding,'E')
    newlon = [newlon maxlon+padEast];
end
if findstr(padding,'S')
    newlat = [minlat-padSouth newlat];
end
if findstr(padding,'N')
    newlat = [newlat maxlat+padNorth];
end

% interpolate for plotting
[y,x] = meshgrid(newlon,newlat);
newelev = interp2(lon,lat,elev,y,x);
pcolor(y,x,newelev); colorbar; caxis(1000*[bottom top]);

% interpolate for elevation comparisons at cell centers
ctrlon = newlon(1:end-1) + diff(newlon)/2;
ctrlat = newlat(1:end-1) + diff(newlat)/2;
[y,x] = meshgrid(ctrlon,ctrlat);
newelev = interp2(lon,lat,elev,y,x);

z = [top top-cumsum(dz)];
if findstr(padding,'U')
    z = [top+padUp z];
end
if findstr(padding,'D')
    z = [z bottom-padDown];
end
ctrz = 1000*(z(1:end-1) + diff(z)/2);

cov = zeros(length(ctrlat),length(ctrlon),length(ctrz));
v   = zeros(length(ctrlat),length(ctrlon));

for k = 1:length(ctrz)
    if ctrz(k)>0
        iAir = find(newelev<ctrz(k)); v(iAir) = 0;
        iEarth = find(newelev>=ctrz(k)); v(iEarth) = 1;
    else       
        iOcean = find(newelev<ctrz(k)); v(iOcean) = 9;
        iEarth = find(newelev>=ctrz(k)); v(iEarth) = 1;
    end
    cov(:,:,k) = v;
end

% for k = 1:nz+1
%     figure; pcolor(y,x,cov(:,:,k)); colorbar
% end