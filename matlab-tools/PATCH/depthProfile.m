function result = depthProfile(CondObj,patchStruct,patchNum)

%    Inputs are Cond_xy_ll object containing conductivity model, and patch
%    struture

%   first need to interpolate log cond onto patch coordinates
%   --> lat/lon of patch cell centers
[ii,jj]  = find(patchStruct.mask==patchStruct.mask_nums(patchNum));
[n,m]  = size(patchStruct.mask);
lonRange = patchStruct.ll_lims(2)-patchStruct.ll_lims(1);
latRange = patchStruct.ll_lims(4)-patchStruct.ll_lims(3);
dLat = latRange/m;
dLon = lonRange/n;
lat = patchStruct.ll_lims(3)+(jj-.5)*dLat;
lon = patchStruct.ll_lims(1)+(ii-.5)*dLon;
nll = length(lat);
ll = [lat lon];
xy = ll2xy(CondObj,ll');
nz = length(CondObj.z);
Z = ones(size(lat))*CondObj.z';
X = xy(1,:)'*ones(1,nz);
Y = xy(2,:)'*ones(1,nz);
sigPatch = interp3(CondObj.Y,CondObj.X,CondObj.Z,CondObj.CONDxyz,Y,X,Z);


%   generate easily rotatble grid (defined on a circle, centered on patch)
%   for evaluation of anisotropic conductivity
%   note that xy is in km
xy(1,:) = xy(1,:)-mean(xy(1,:));
xy(2,:) = xy(2,:)-mean(xy(2,:));
rhoMax = sqrt(max(sum(xy.*xy,1)));

Ny = ceil(2*sqrt(nll/3));
dy = rhoMax/Ny;
yi = -rhoMax:dy:rhoMax;
n = length(yi);
xi = yi;
yi = yi'*ones(1,n);
xi = ones(n,1)*xi;
r = sqrt(yi.*yi+xi.*xi);
phi = atan2(xi,yi);

%  loop over angles  ...
nTheta = 20;
condAvg = zeros(nz,nTheta);
dTheta = pi/nTheta;

%  "anisotropic" avergae as a function of azimuth
for iz = 1:nz
    Fsig = TriScatteredInterp(xy',squeeze(sigPatch(:,iz)));
    k = 0;
    for theta = 0:dTheta:pi-dTheta
        k = k+1;
        yi = r.*cos(theta+phi);
        xi = r.*sin(theta+phi);
        logSig = Fsig(xi,yi);
        rho = 10.^(-logSig);
        use = isnan(rho);
        rho(use)=0;
        use = 1-use; 
        rhoAvg = sum(rho,1)./sum(use,1);
        useAvg = isnan(rhoAvg);
        sigAvg = 1./rhoAvg;
        sigAvg(useAvg) = 0;
        useAvg = 1-useAvg;
        condAvg(iz,k) = sum(sigAvg)./sum(useAvg);
    end
end
condAvg = log10(condAvg); 

%% fit avg, + sine cosine
TH =  0:dTheta:pi-dTheta;
X = [ones(length(TH),1), cos(2*TH)' sin(2*TH)'];
condFit = X\condAvg';
thetaMax = atan2(condFit(3,:)',condFit(2,:)')/2;
thetaMax = thetaMax*180/pi;
sigmaMax = sqrt(condFit(2,:).^2+condFit(3,:).^2);
sigmaMin = condFit(1,:)-sigmaMax;
sigmaMax = condFit(1,:)+sigmaMax;
result = struct('condAvg',condAvg,'sigmaMax',sigmaMax,'sigmaMin',sigmaMin,...
    'thetaMax',thetaMax,'z',CondObj.z,'TH',TH,'xy',xy,'sigPatch',sigPatch);


    