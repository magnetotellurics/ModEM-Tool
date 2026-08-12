function [E1D,model1d] = PrimaryField(dataFile3d,modelFile3d,modelFile1d,sphHarm,falseGrid)
% Function to 1) create an imitation of 2 MT modes, 2) run the 1D global
% forward modeling in Matlab and 3) write it out as primary fields by
% extracting the E-fields that are centered around (0,0) and putting them
% on a spherical grid that is appropriate for the actual model location.
% A. Kelbert, 27 Dec 2022; last mod. 10 Jul 2023 for general sources
%    BUT this is still hard-coded for exactly two sources only

%% create the P10 spherical source with 2 modes - unless source provided
if nargin < 4
    % scaled for a halfspace of 100 Ohmm at 1000 secs (old air layers)
    scalingX = 200;
    scalingY = 200;
    temp = SHCreateVec(1);
    P10 = SHSetValue(temp,1,1,0);
    ModeX = SHRotateVec(P10,0,45,-45)*scalingX;
    ModeY = -P10*scalingY;
    shc(1,:) = ModeX(2:end);
    shc(2,:) = ModeY(2:end);
    %figure; sp=SHPlotProj(ModeX,1);
    
    % overwrite with the complex coeff that seem correct
    shc(1,:) = [0.0000 + 0.0000i 0.0000 - 0.5000i 0.0000 + 0.5000i];
    shc(2,:) = [-1 0 0];

else
    shc = sphHarm;
    
end

% by default, falseGrid is true (to be used for MT)
if nargin < 5
    
    falseGrid = 1;
    
end

%Period = 10.^(0.6:0.6:4); % 6 periods from the unit_testing_setup.m

%% for testing, initially start with 1 period - skip the constant but
% but we no longer rotate to the center of the model, instead we extract
% an area centered around (0,0) from the 1D field for primary
% These files are not currently used but will be used for internal global
% 1D modeling in Fortran.
%shcwrite(ModeX(2:end),Period,sprintf('SphericalPlaneWave_ModeX_%02dper',length(Period)));
%shcwrite(ModeY(2:end),Period,sprintf('SphericalPlaneWave_ModeY_%02dper',length(Period)));

%%
if ischar(dataFile3d)
    dat = mtdata.read(dataFile3d,'list','[mV/km]/[nT]','Full_Impedance');
    dat = mtdataplot(dat);
    Period = dat.v.per;
elseif isnumeric(dataFile3d)
    Period = dataFile3d;
else
    warning('Set of input periods or the data file name not recognized');
end

[modelFilePath,modelFileName,ext] = fileparts(modelFile3d);
if strcmp(ext,'.mat')
    load([modelFile3d]);
    if ~exist('model3d','var')
        warning('Please make sure that the input Matlab file contains model3d');
    end
elseif strcmp(ext,'.nc')
    model3d = llmodel.read(modelFile3d,'netcdf');
elseif strcmp(ext,'.rho')
    model3d = llmodel.read(modelFile3d,'ModEM');
else
    warning('Unknown input model file type in PrimaryField.m');   
end

if nargin >= 3
    [model1dFilePath,model1dFileName,ext] = fileparts(modelFile1d);
    if strcmp(ext,'.mat')
        load(modelFile1d);
        if ~exist('model1d','var')
            warning('Please make sure that the input Matlab file contains model1d');
        end
    elseif strcmp(ext,'.nc')
        model1d = llmodel.read(modelFile1d,'netcdf');
    elseif strcmp(ext,'.rho')
        model1d = llmodel.read(modelFile1d,'ModEM');
    else
        warning('Unknown input model file type for 1D model in PrimaryField.m');
    end
else
    model1d = model3d;
    for k=1:size(model1d.v,3)
        model1d.v(:,:,k) = mean(mean(model3d.v(:,:,k)));
    end
end
layer1d = 1e3*model1d.grid.depth(1:end);
cond1d = 10.^squeeze(squeeze(model1d.v(1,1,:)));

%%
field = {'Es','Er'};
origin = [model3d.grid.lat(1) model3d.grid.lon(1) 0];
dx = model3d.grid.dlat;
dy = model3d.grid.dlon;
dz = model3d.grid.dz*1000;
x = model3d.grid.lat;
y = model3d.grid.lon;
z = model3d.grid.depth*1000;
ctrx = model3d.grid.lat(1:end-1) + model3d.grid.dlat/2;
ctry = model3d.grid.lon(1:end-1) + model3d.grid.dlon/2;
nx = length(dx);ny = length(dy);
grd = TGrid3D(dx,dy,dz,0,origin);
airthck = 10^6;
nzAir = 12;
grd = setAirLayers(grd,'method','fixed height','maxheight',airthck,'nlayers',nzAir);
% use the old default for air layers for now, 'mirror 10 3. 30.' 
%grd = setAirLayers(grd,'method','mirror','nlayers',nzAir);
sz = nzAir + 1;
ez = nzAir + model3d.grid.nzCrust + model3d.grid.nzEarth;

% TRUE lat/lon centered at lat=center(1) lon=center(2)
center(1) = model3d.grid.lat(1)+sum(model3d.grid.dlat)/2;
center(2) = model3d.grid.lon(1)+sum(model3d.grid.dlon)/2;
center(3) = 0;
truelatex = ctrx;
truelonex = y;
truelatey = x;
trueloney = ctry;

if falseGrid
    % FALSE lat/lon centered at lat=0 lon=0 (in [-180,180] range)
    falseorigin(1) = -sum(model3d.grid.dlat)/2;
    falseorigin(2) = -sum(model3d.grid.dlon)/2;
    falseorigin(3) = 0;
    gridshift = falseorigin - origin;
    latex = ctrx+gridshift(1);
    lonex = y+gridshift(2);
    latey = x+gridshift(1);
    loney = ctry+gridshift(2);
else
    % TRUE lat/lon used throughout
    latex = truelatex;
    lonex = truelonex;
    latey = truelatey;
    loney = trueloney;
end

alt_air = flip(cumsum(grd.Dz(grd.Nza:-1:1)));
alt =-1*cumsum([0;grd.Dz(sz:ez)]);
depth_earth = cumsum([0;grd.Dz(grd.Nza+1:end)]);
alt_c = (alt(1:end-1)+alt(2:end))/2.0;

nrAir = length(alt_air);
nr = length(alt);
nt = 90;
lat1D = 0.5:1:180; % global co-latitude
lon1D = 0.5:1:360; % global longitude
sig = [cond1d;cond1d(end)];
ds = [layer1d,sig];
nPer = length(Period);
nMode = 2;
Es = cell(nPer,nMode);
Er = cell(nPer,nMode);
Hs = cell(nPer,nMode);
Hr = cell(nPer,nMode);

%% TS or layered 1D modeling for multiple periods, multiple modes
for iPer = 1:nPer
    mm = TSModel(ds,0,Period(iPer),nt,'uniform',1e-4);
    for iMode = 1:nMode
        scaling = Period(iPer)/5;
        [Es{iPer,iMode}]=mm.ShcInc(shc(iMode,:)*scaling,[alt_air;alt],{'Es'});
        Er{iPer,iMode}(1:nr,1)=radField(mm);
        [~,Hr{iPer,iMode},Hs{iPer,iMode}]=mm.ShcInc(shc(iMode,:)*scaling,0,{'Hr','Hs'});
    end
end

%% save E field from TS model
% % % convert to regional ModEMM EM field file
[Lon1D,Lat1D] = meshgrid(lon1D, lat1D);
[Lonex,Latex] = meshgrid(lonex, latex);
[Loney,Latey] = meshgrid(loney, latey); % latey(2:end-1)
[Lonez,Latez] = meshgrid(lonex, latey); % latey(2:end-1)

[Lonhx,Lathx] = meshgrid(loney,latey); % latey(2:end-1)
[Lonhy,Lathy] = meshgrid(lonex,latex);
[Lonhz,Lathz] = meshgrid(loney,latex);

Ex= cell(nPer,nMode);Ey= cell(nPer,nMode);Ez = cell(nPer,nMode);

for iPer = 1:nPer
    
    nlayers = length(Es{1});
    for iMode = 1:nMode
        for ir= 1:nlayers
            foox = -1i*flip(Es{iPer,iMode}(ir).ft,1);
            fooy = 1i*flip(Es{iPer,iMode}(ir).fp,1);
%             if iPer==1 && ir==nrAir+1
%                 figure; pcolor(Lon1D'-180,90-Lat1D',squeeze(real(foox)')); colorbar
%                 figure; pcolor(Lon1D'-180,90-Lat1D',squeeze(imag(foox)')); colorbar
%                 figure; pcolor(Lon1D'-180,90-Lat1D',squeeze(real(fooy)')); colorbar
%                 figure; pcolor(Lon1D'-180,90-Lat1D',squeeze(imag(fooy)')); colorbar
%             end
            
            Ex{iPer,iMode}(:,:,ir) = interp2(Lon1D-180,90-Lat1D,foox,Lonex,Latex, 'spline');
            Ey{iPer,iMode}(:,:,ir) = interp2(Lon1D-180,90-Lat1D,fooy,Loney,Latey, 'spline');
        end
    end
    
    nlayers = length(Er{1});
    for iMode = 1:nMode
        for ir= 1:nlayers
            fooz = -1i*flip(Er{iPer,iMode}(ir).fr,1);
            Ez{iPer,iMode}(:,:,ir) = interp2(Lon1D-180,90-Lat1D,fooz,Lonez,Latez, 'spline');
        end
    end
end

% conjugate all electric fields - some conventions somewhere don't match!
% for iPer = 1:nPer
%     for iMode = 1:nMode
%         Ex{iPer,iMode}(:,:,:) = conj(Ex{iPer,iMode}(:,:,:));
%         Ey{iPer,iMode}(:,:,:) = conj(Ey{iPer,iMode}(:,:,:));
%         Ez{iPer,iMode}(:,:,:) = conj(Ez{iPer,iMode}(:,:,:));
%     end
% end
for iPer = 1:nPer
    for iMode = 1:nMode
        Ex{iPer,iMode}(:,:,:) = 1i*(Ex{iPer,iMode}(:,:,:));
        Ey{iPer,iMode}(:,:,:) = 1i*(Ey{iPer,iMode}(:,:,:));
        Ez{iPer,iMode}(:,:,:) = 1i*(Ez{iPer,iMode}(:,:,:));
    end
end

%save('NorthAmericaSynthetic_Global_E1D.mat','Ex','Ey','Ez');

%% compute magnetic fields
% Hx= cell(nPer,nMode);Hy= cell(nPer,nMode);Hz = cell(nPer,nMode);
% for iPer = 1:nPer
%     for iMode = 1:nMode
%         foohx = -1*flip(Hs{iPer,iMode}.ft,1);
%         foohy = flip(Hs{iPer,iMode}.fp,1);
%         foohz = -1*flip(Hr{iPer,iMode}.fr,1);
%         
%         Hx{iPer,iMode}(:,:) = interp2(Lon1D,Lat1D,foohx,Lonhx,Lathx, 'spline');
%         Hy{iPer,iMode}(:,:) = interp2(Lon1D,Lat1D,foohy,Lonhy,Lathy, 'spline');
%         Hz{iPer,iMode}(:,:) = interp2(Lon1D,Lat1D,foohz,Lonhz,Lathz, 'spline');
%     end
% end

%% write out E0 file (interpolate to the false grid but write true grid)
E0 = TVector3D_SG(grd);
[Lonz,Latz,Zz] = meshgrid(lonex, latey, alt); % latey(2:end-1)
[Lonz_c,Latz_c,Zz_c] = meshgrid(lonex, latey, alt_c); % latey(2:end-1)
for iPer = 1:nPer
    for iMode = 1:nMode
        E0.x(:,:,:)= Ex{iPer,iMode};
        E0.y(:,:,:)= Ey{iPer,iMode};
        E0.z(:,:,sz:ez) = interp3(Lonz,Latz,Zz,Ez{iPer,iMode},Lonz_c,Latz_c,Zz_c,'spline');
        SOLN0.E{iMode,iPer} = E0;
    end
end
SOLN0.S.CondRead = 1;
SOLN0.grid = E0.grid;
SOLN0.grid.rotation = 0.0;
SOLN0.grid.units = 'm';
SOLN0.Modes = {'X','Y'}; % '1'
SOLN0.Periods = Period;
plotEMsoln(SOLN0);

E1D = SOLN0;

end

%cd(modelFilePath);
%save(['NorthAmericaSynthetic',downsample,'.E1D_matlab.mat'],'SOLN0','-v7.3');
%filename = ['NorthAmericaSynthetic',downsample,'.E1D_matlab.soln'];
%wtExyz(SOLN0,filename)

%I/O testing
%SOLN0_read = readE(filename,1);
%plotEMsoln(SOLN0_read)

%model3d.write(['NorthAmericaSynthetic',downsample,'.3D.rho'],'ModEM');
%model1d.write(['NorthAmericaSynthetic',downsample,'.1D.rho'],'ModEM');

%save(['USA',downsample,'.E1D_matlab.mat'],'SOLN0','-v7.3');
%filename = ['USA',downsample,'.E1D_matlab.soln'];
%wtExyz(SOLN0,filename)
%model3d.write(['USA',downsample,'.3D.rho'],'ModEM');
%model1d.write(['USA',downsample,'.1D.rho'],'ModEM');

%%
%filenameSFF = 'NorthAmericaSynthetic.downsample.halfdegree.output.soln';
%SOLN0_SFF = readE(filenameSFF,nPer);
%plotEMsoln(SOLN0_SFF);

%%
% cd('~/Developer/ModEM-OSU/branches/stable-/unit_testing/3D/');
% halfspaceC = readE('Halfspace_6per_newAir_depth1000km_stable.soln',1);
% plotEMsoln(halfspaceC)
% 
% %% 
% cd('~/Developer/ModEM-OSU/branches/stable-/unit_testing/3D/');
% dat = mtdata.read('Halfspace_1per_newAir_depth1000km_stable.dat','list','[mV/km]/[nT]','Full_Impedance');
% dat = mtdataplot(dat);
% dat.plotap('xy',1);
% rhoxy = 0.2*dat.v.per(1)*abs(mean(dat.v.data(:,1,2)))^2;
% z = 7.908043E+00+1i*7.920041E+00; per = 10^0.6; % 4 secs
% rho = 0.2*per*abs(z)^2
% phase = rad2deg(atan(imag(z)/real(z)))
% 
% cd('~/Developer/ModEM-OSU/branches/stable-/unit_testing/3D_spherical/');
% dats = mtdata.read('Halfspace+Spherical_6per_SFF_depth1000km_output.dat','list','[mV/km]/[nT]','Full_Impedance');
% dats = mtdataplot(dats);
% dats.plotap('xy',1);
% rhoxys = 0.2*dats.v.per(1)*abs(mean(dats.v.data(:,1,2)))^2;
% z = 7.924457E+00+1i*7.937456E+00; per = 10^0.6; % 4 secs
% rhos = 0.2*per*abs(z)^2
% phases = rad2deg(atan(imag(z)/real(z)))

%% read SFF outputs
% SOLN0_read = readE(filename,1);
% plotEMsoln(SOLN0_read);

% USAFWD = readE('USA-2021-08-05-0.5x1-large-MT-FWD.soln',1);
% plotEMsoln(USAFWD);
% 
% USASFF = readE('USA-2021-08-05-0.5x1-large-MT-SFF.soln',1);
% plotEMsoln(USASFF);
% 
% USAGDV = readE('USA-2021-08-05-0.5x1-large-GlobalDV-SFF.soln',1);
% plotEMsoln(USAGDV);
