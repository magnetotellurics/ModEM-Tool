% [x,y,z,v] = mkGrid3D(-119,-106.5,41,46.5,top,-6,nlon,nlat,12,'NSEW'); 
% NOTE: UNABLE TO FIT THE DATA WITH LOGARITHMIC VERTICAL SPACING;
% THE INTERPOLATED 14 KM MODEL PRODUCES VERY DIFFERENT RESPONSES.
% HERE, USE NON-LOGARITHMIC VERTICAL SPACING WITH 10 KM GRID,
% BUT MAKE IT GO TO 900 KM DEPTH.
clear ll SRPYgrid ans
lims.latmin = 41;
lims.latmax = 46.5;
lims.lonmin = -119;
lims.lonmax = -106.5;
nlat = 60;
nlon = 98;
ll = llgrid.empty(8,0);
% optional topography to 3 km: use 6 layers for vertical spacing 0.5 km
ll(1) = llgrid(lims,-3,0,nlat,nlon,0,'NSEW');
% vertical spacing 0.5 km
ll(2) = llgrid(lims,0,6,nlat,nlon,12,'NSEW');
% increase vertical spacing to 1 km
ll(3) = llgrid(lims,6,10,nlat,nlon,4,'NSEW');
% vertical spacing 2 km
ll(4) = llgrid(lims,10,18,nlat,nlon,4,'NSEW');
% vertical spacing 4 km
ll(5) = llgrid(lims,18,50,nlat,nlon,8,'NSEW');
% the bulk of the model with 8 km vertical spacing
ll(6) = llgrid(lims,50,210,nlat,nlon,20,'NSEW');
% the bottom layers with 16 km vertical spacing
ll(7) = llgrid(lims,210,306,nlat,nlon,6,'NSEW');
% now, make them logarithmic until 900 km depth
ll(8) = llgrid(lims,306,900,nlat,nlon,16,'NSEW',1.05);
% merge the layers into one big grid
SRPYgrid = llgrid.merge(ll);
[COV,ELEV] = mask(SRPYgrid);
plot(SRPYgrid,ELEV);

%%
plot(SRPYgrid,COV,0);
writeCov_3D('Yellowstone_10km.cov',COV);

%% set the origin to match that of the data set
lat0 = 42.016;
lon0 = - 116.477;
obj = SRPYgrid;
grid = xygrid(SRPYgrid,lat0,lon0);
temp = llgrid(grid,lat0,lon0);
temp.lon - obj.lon

%% try grid center for origin
% lat0 = 43.7500;
% lon0 = -112.7500;
% grid = xygrid(SRPYgrid,lat0,lon0);
% temp = llgrid(grid,lat0,lon0);

%%

% map covariance to prior resistivity values
% two options: uniform halfspace

prior = xymodel.prior(grid,'uniform',200);
write('Yellowstone_10km_200ohmm_pr.ws',prior,'WS');

% and radially increasing with depth

prior = xymodel.prior(grid,'radial',[200 1]);
write('Yellowstone_10km_radial_pr.ws',prior,'WS');
uiplot(prior);

%%
dat = mtdata.read('Yellowstone_14freq_paper_errfl5T3.dat','list','[mV/km]/[nT]','Full_Impedance');
dat = mtdataplot(dat); 
%dat.write('temp.dat');
apresplt(dat,'ORG10');
%newdat=apres(dat);
res = mtdata.read('Yellowstone_10km_errfl5T3_smooth_NLCG_049.dat','list','[mV/km]/[nT]','Full_Impedance','predicted');
res = mtdataplot(res); 
obj = mtdatacompare(dat,res);
rms = obj.misfit;
obj.apresplt('ORG10');
obj.impplt('ORG10');
obj.plot(4,'ZXY',6);

