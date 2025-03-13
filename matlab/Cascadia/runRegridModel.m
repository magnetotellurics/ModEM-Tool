%   regrid cascadia inverse solution
%CascInvDir = '/home/server/pi/homes/egbert/EarthScope/INVERSION/';
%Dir = '/home/gauss/egbert/progs/EM3D/EM3Dnew/src/ModEM/examples/3D_MT/';
%ModelNameIn =  [ CascInvDir 'caswithprior_model.02_01'];
%ModelNameOut = [Dir 'CascInv/cascInv2Course.wsm'];
Dir = '../';
ModelNameIn = 'original/caswithprior_model.02_01';
ModelNameOut = 'test.ws';

[dx,dy,dz,rhoIn,nzAir,type,origin,rotation] = read_WS3d_model(ModelNameIn);
gridIn = struct('dx',dx,'dy',dy,'dz',dz,'origin',origin);
maskIn = zeros(size(rhoIn));
maskIn(rhoIn<=0.5) = 1;
dx2 = [dx(1:5);  24300*ones(36,1); dx(end-4:end)];
dy2 = [dy(1:6);  25125*ones(32,1); dy(end-5:end)];
gridOut = struct('dx',dx2,'dy',dy2,'dz',dz,'origin',origin);

[rhoOut,maskOut] = regridModel(gridIn,gridOut,rhoIn,maskIn);

origin = zeros(size(origin));
gridOut.origin = origin;
status = write_WS3d_model(ModelNameOut,dx2,dy2,dz, ...
        rhoOut,nzAir,type,origin,rotation)

%   Load site locations ... grid coordinate system

DataFileIn = [ Dir 'Cascadia/cascad3.imp'];
[allData,info,units,isign] = readZ_3D(DataFileIn);
siteLoc = allData{1}.siteLoc;

x1 = gridOut.origin(1)+[0; cumsum(gridOut.dx)];
y1 = gridOut.origin(2)+[0; cumsum(gridOut.dy)];
z1 = gridOut.origin(3)+[0; cumsum(gridOut.dz)];

x = [0; cumsum(gridIn.dx)];
y = [0; cumsum(gridIn.dy)];
z = [0; cumsum(gridIn.dz)];

figure; pcolor(y,x,zeroPad(squeeze(rhoIn(:,:,1))));
hold on;
plot(siteLoc(:,2),siteLoc(:,1),'y*','markersize',10);

[iOcean,ix,iy] = findSites2D(siteLoc,x,y,rhoIn(:,:,1),0.35)

figure; pcolor(y1,x1,zeroPad(squeeze(rhoOut(:,:,1))));
hold on;
plot(siteLoc(:,2),siteLoc(:,1),'y*','markersize',10);

[iOcean,ix1,iy1] = findSites2D(siteLoc,x1,y1,rhoOut(:,:,1),0.35)
plot(siteLoc(iOcean,2),siteLoc(iOcean,1),'wo','markersize',10);

% 54, 75, 76: move these 3 sites manually inland
siteLoc(54,:) = [siteLoc(54,1), y1(iy1(1)+1)+5000, 0];
siteLoc(75,:) = [x1(ix1(2))-5000, siteLoc(75,2), 0];
siteLoc(76,:) = [x1(ix1(3))-5000, siteLoc(76,2), 0];

for i = 1:length(allData)
    allData{i}.siteLoc = siteLoc;
end

DataFileOut = [ Dir 'Cascadia/cascad3a.imp'];
writeZ_3D(DataFileOut,allData,info,units,isign);


