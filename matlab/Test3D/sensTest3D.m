%  Script sensTest is used to test sensitivity calculation in fortran 90
%% settings and file directory
SENS_MATRIX = 0;
testDir = 'SENSTEST/';
config = [testDir 'config.mat'];
%% create conf structure
if exist(config,'file')
    load(config)
else
    Config3D
end
LOGCOND = conf.LOGCOND;

%% File names
m0_File = [testDir 'm0.ws'];
m1_File = [testDir 'm1.ws'];
dm_File = [testDir 'dm.ws'];
d0_File = [testDir 'd0.dat'];
e0_File = [testDir 'e0.soln'];
d1_File = [testDir 'd1.dat'];
d_File = [testDir 'd.dat'];
Jxdm_File = [testDir 'Jxdm.dat'];
JTxd_File = [testDir 'JTxd.ws'];
J_File = [testDir 'J.sns'];
if exist(testDir,'dir')==0
    mkdir(testDir);
end

%% create the background model and data files and structures
conf.COMPUTE_DATA = 1;
conf.ADD_NOISE = 0;
conf.PLOT_MODEL = 1;
conf.PLOT_DATA = 1;
[m0,d0] = mkCond3D(conf,m0_File,d0_File,e0_File);

%% save the (manually, if required) edited configuration
if conf.SAVE_CONF
    save(config,'conf')
end

%% create a random perturbation in the background model
dm = dCond_3D(m0,0.0002);
writeCond_3D(dm_File,dm,2);

% create the perturbed model and corresponding data
m1 = addCond_3D(m0,dm,LOGCOND);
writeCond_3D(m1_File,m1,2);
Test3D('FORWARD',m1_File,d0_File,d1_File)
d1 = readZ_3D(d1_File);

% Compute perturbation of solution ...
%    (divide by epsilon to compare to Jxdm)
nPer = length(d0);
d = d0;
for iPer = 1:nPer
   d{iPer}.Z =  (d1{iPer}.Z-d0{iPer}.Z);
end
writeZ_3D(d_File,d);

% forward sensitivity calculation
Test3D('MULT_BY_J',m0_File,dm_File,d0_File,Jxdm_File)
[Jxdm] = readZ_3D(Jxdm_File);

% backward sensitivity calculation
Test3D('MULT_BY_J_T',m0_File,d_File,JTxd_File)
[JTxd] = readCond_3D(JTxd_File,2);
%% if required, compute the full sensitivity matrix and products
if (SENS_MATRIX)
    % create the matrix of sensitivities
    Test3D('COMPUTE_J',m0_File,d0_File,J_File)
    % J = cell(2,1);
    % J{1} = d0;
    % J{2} = readCondMTX_3D(J_File);
    J = readCondMTX_3D(J_File);

    % compute Jxdm from sensitivity matrix ...
    JxdmComp = Jxdm;
    if LOGCOND
        temp = dm.v
    else
       %   compute Jxdm from sensitivity matrix, using
       %   --> sensitivity matrix for linear conductivity
       %   --> perturbation dm of log conductivity
       %   --> background conductivity ...
       %       should match Jxdm from log case (to within rounding)
        temp = dm.v .* m0.v
    end
    ii = 0;
    for iPer = 1:nPer
        JxdmComp{iPer}.Z(:,:) = 9999;
        nSites = size(Jxdm{iPer}.Z,2);
        nComp = size(Jxdm{iPer}.Z,1);
        for iSite = 1:nSites
            Zr(1:nComp) = 0; Zi(1:nComp) = 0;
            for iComp = 1:nComp
                ii = ii + 1;
                Zr(iComp) =  sum(sum(sum(J{ii}.v .* temp )));
                ii = ii + 1;
                Zi(iComp) =  sum(sum(sum(J{ii}.v .* temp )));
            end
            JxdmComp{iPer}.Z(:,iSite) = Zr+i*Zi;
        end
    end

    % compute JTxd from sensitivity matrix ...
    JTxdComp = JTxd;
    JTxdComp.v = 9999;
    ii = 0;
    for iPer = 1:nPer
        nSites = size(d{iPer}.Z,2);
        nComp = size(d{iPer}.Z,1);
        for iSite = 1:nSites
            Z = d{iPer}.Z(:,iSite);
            for iComp = 1:nComp
               temp = J{ii+1}.v - i*J{ii+2}.v;
               JTxdComp.v = JTxdComp.v + real(Z(iComp)*temp);
               ii = ii + 2;
            end
        end
    end
end
%%  plot comparison
plotBands = [1 2];
for k = 1:length(plotBands)
   T = d0{plotBands(k)}.T;
   figure('Position',[100,100,400,700])
   axes('Position',[.1,.55,.8,.35])
   plot(real(Jxdm{plotBands(k)}.Z),'ks')
   hold on
   if (SENS_MATRIX)
     plot(real(JxdmComp{plotBands(k)}.Z),'bo')
   end
   plot(real(d{plotBands(k)}.Z),'r+')
   title(['Period = ', num2str(T)])
   axes('Position',[.1,.1,.8,.35])
   plot(imag(Jxdm{plotBands(k)}.Z),'ks')
   hold on
   if (SENS_MATRIX)
     plot(imag(JxdmComp{plotBands(k)}.Z),'bo')
   end
   plot(imag(d{plotBands(k)}.Z),'r+')
end

figure
pcolor(squeeze(JTxd.v(11,:,:)));shading flat; colorbar
set(gca,'FontWeight','demi','FontSize',12)
title('J^T d : computed by call to JmultT')
if (SENS_MATRIX)
    figure
    pcolor(squeeze(JTxdComp.v(11,:,:)));shading flat; colorbar
    set(gca,'FontWeight','demi','FontSize',12)
    title('J^T d : computed in matlab using full sensitiivty')
end