% regiondir = 'regions/Mid-continent';
% dataType = 'Full_Impedance';
% newUnits = '[mV/km]/[nT]';
% cd(regiondir); 
% obj = mtdata.read('./xml-temp','xml',newUnits,dataType);
% %obj = obj.regroup;
% % takes a long time to read the XMLs, so save the object
% save MidcontinentRiftData obj
% write(obj,'MidcontinentRift_no_error_floors.dat','list',0);

obj = mtdata.read('Full_Impedance.dat','list','[mV/km]/[nT]','Full_Impedance');
dat = mtdataplot(obj);
%dat = dat.regroup;
%apresplt(dat,'ORG10'); 


% [x,y,z,v] = mkGrid3D(-119,-106.5,41,46.5,top,-6,nlon,nlat,12,'NSEW'); 
% NOTE: UNABLE TO FIT THE DATA WITH LOGARITHMIC VERTICAL SPACING;
% THE INTERPOLATED 14 KM MODEL PRODUCES VERY DIFFERENT RESPONSES.
% HERE, USE NON-LOGARITHMIC VERTICAL SPACING WITH 10 KM GRID,
% BUT MAKE IT GO TO 900 KM DEPTH.
clear ll mygrid ans
lims.latmin = min(obj.v.lat);
lims.latmax = max(obj.v.lat);
lims.lonmin = min(obj.v.lon);
lims.lonmax = max(obj.v.lon);
nlat = 62;
nlon = 62;
ll = llgrid.empty(2,0);
% optional topography to 3 km: use 6 layers for vertical spacing 0.5 km
%ll(1) = llgrid(lims,-3,0,nlat,nlon,0,'NSEW');

% vertical spacing 0.5 km
ll(1) = llgrid(lims,0,6,nlat,nlon,12,'NSEW');
% now, make them logarithmic until 900 km depth
ll(2) = llgrid(lims,2,900,nlat,nlon,78,'NSEW',1.05);
% merge the layers into one big grid
mygrid = llgrid.merge(ll);
[COV,ELEV] = mask(mygrid);
plot(mygrid,ELEV);

%%
plot(mygrid,COV,0);
writeCov_3D('MidcontinentRift.cov',COV);

%% set the origin to match that of the data set
lat0 = dat.origin(1);
lon0 = dat.origin(2);
%obj = mygrid;
grid = xygrid(mygrid,lat0,lon0);
% temp = llgrid(grid,lat0,lon0);
% temp.lon - obj.lon

%% try grid center for origin
% lat0 = 43.7500;
% lon0 = -112.7500;
% grid = xygrid(SRPYgrid,lat0,lon0);
% temp = llgrid(grid,lat0,lon0);

%%

% map covariance to prior resistivity values
% two options: uniform halfspace

prior = xymodel.prior(grid,'uniform',200);
write('MidcontinentRift_20km_200ohmm_pr.ws',prior,'WS');

% and radially increasing with depth

prior = xymodel.prior(grid,'radial',[200 1]);
write('MidcontinentRift_20km_radial_pr.ws',prior,'WS');
uiplot(prior);

%%
% dat = mtdata.read('Yellowstone_14freq_paper_errfl5T3.dat','list','[mV/km]/[nT]','Full_Impedance');
% dat = mtdataplot(dat); 
% %dat.write('temp.dat');
% apresplt(dat,'ORG10');
% %newdat=apres(dat);
% res = mtdata.read('Yellowstone_10km_errfl5T3_smooth_NLCG_049.dat','list','[mV/km]/[nT]','Full_Impedance','predicted');
% res = mtdataplot(res); 
% obj = mtdatacompare(dat,res);
% rms = obj.misfit;
% obj.apresplt('ORG10');
% obj.impplt('ORG10');
% obj.plot(4,'ZXY',6);

