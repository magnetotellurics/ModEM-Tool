%  Script sensTest is used to test sensitivity calculation in fortran 90
%  NOTE: Test2D must be in your path or current directory
%% settings and file directory
SENS_MATRIX = 0;
testDir = 'SENSTEST/';
config = [testDir 'config.mat'];
%% create conf structure
if exist(config,'file')
    load(config)
else
    Config2D
end
MODE = conf.mode;
LOGCOND = conf.LOGCOND;

%% File names
m0_File = [testDir 'm0.rho'];
m1_File = [testDir 'm1.rho'];
dm_File = [testDir 'dm.rho'];
d0_File = [testDir 'd0_' MODE '.dat'];
d1_File = [testDir 'd1_' MODE '.dat'];
d_File = [testDir 'd_' MODE '.dat'];
Jxdm_File = [testDir 'Jxdm_' MODE '.dat'];
JTxd_File = [testDir 'JTxd_' MODE '.rho'];
J_File = [testDir 'J_' MODE '.sns'];
if exist(testDir,'dir')==0
    mkdir(testDir);
end

%% create the background model and data files and structures
conf.COMPUTE_DATA = 1;
conf.ADD_NOISE = 0;
conf.PLOT_MODEL = 1;
conf.PLOT_DATA = 1;
[m0,d0] = mkCond2D(conf,m0_File,d0_File)

%% save the (manually, if required) edited configuration
if conf.SAVE_CONF
    save(config,'conf')
end

%% create a random perturbation in the background model
dm = dCond_2D(m0,0.001);
writeCond_2D(dm_File,dm);

% create the perturbed model and corresponding data
m1 = addCond_2D(m0,dm,LOGCOND);
writeCond_2D(m1_File,m1);
Test2D('FORWARD',m1_File,d0_File,d1_File)
d1 = readZ_2D(d1_File);

% Compute perturbation of solution ...
%    (divide by epsilon to compare to Jxdm)
nPer = length(d0);
d = d0;
for iPer = 1:nPer
   d{iPer}.Z =  (d1{iPer}.Z-d0{iPer}.Z);
end
writeZ_2D(d_File,d);

% forward sensitivity calculation
Test2D('MULT_BY_J',m0_File,dm_File,d0_File,Jxdm_File)
[Jxdm] = readZ_2D(Jxdm_File);

% backward sensitivity calculation
Test2D('MULT_BY_J_T',m0_File,d_File,JTxd_File)
[JTxd] = readCond_2D(JTxd_File);

cax = [-8,8];  %   log10
ctitle = 'JTxd';
OPTIONS = struct('nYskip',7,'nZplot',20,'cax',cax,'title',ctitle);
plotCond(JTxd,JTxd.grid,OPTIONS);
%% if required, compute the full sensitivity matrix and products
if (SENS_MATRIX)
    % create the matrix of sensitivities
    Test2D('COMPUTE_J',m0_File,d0_File,J_File)
    % J = cell(2,1);
    % J{1} = d0;
    % J{2} = readCondMTX_2D(J_File);
    J = readCondMTX_2D(J_File);
    
    % take care or the units conversion
    f = ImpUnits('[V/m]/[T]',d{1}.units);

    % compute Jxdm from sensitivity matrix ...
    JxdmComp = Jxdm;
    ii = 0;
    for iPer = 1:nPer
        nSites = length(Jxdm{iPer}.Z);
        for iSite = 1:nSites
            ii = ii + 1;
            Zr =  sum(sum(J{ii}.v .* dm.v ));
            ii = ii + 1;
            Zi =  sum(sum(J{ii}.v .* dm.v ));
            JxdmComp{iPer}.Z(iSite) = Zr+1i*Zi;
        end
    end

    %   compute Jxdm from sensitivity matrix, using
    %   --> sensitivity matrix for linear conductivity
    %   --> perturbation dm of log conductivity
    %   --> background conductivity ...
    %       should match Jxdm from log case (to within rounding)

    if ~LOGCOND
        %  perturbation to log conductivity
        temp = dm.v.*m0.v;
        JxdmCompLog = Jxdm;
        ii = 0;
        for iPer = 1:nPer
            nSites = length(Jxdm{iPer}.Z);
            for iSite = 1:nSites
                ii = ii + 1;
                Zr =  sum(sum(J{ii}.v .* temp ));
                ii = ii + 1;
                Zi =  sum(sum(J{ii}.v .* temp ));
                JxdmCompLog{iPer}.Z(iSite) = Zr+1i*Zi;
            end
        end
    end

    % compute JTxd from sensitivity matrix ...
    JTxdComp = JTxd;
    JTxdComp.v = 0;
    ii = 0;
    for iPer = 1:nPer
        nSites = length(d{iPer}.Z);
        for iSite = 1:nSites
            Z = d{iPer}.Z(iSite);
            temp = J{ii+1}.v - 1i*J{ii+2}.v;
            JTxdComp.v = JTxdComp.v + real(Z*temp);
            ii = ii + 2;
        end
    end
    cax = [-8,8];  %   log10
    ctitle = 'JTxdComp';
    OPTIONS = struct('nYskip',7,'nZplot',20,'cax',cax,'title',ctitle);
    plotCond(JTxdComp,JTxdComp.grid,OPTIONS);
    
    % plot the difference
    JTdiff = JTxd;
    JTdiff.v = JTxdComp.v - JTxd.v;
    cax = [-1,1];  %   log10
    ctitle = 'JTxdComp - JTxd';
    OPTIONS = struct('nYskip',7,'nZplot',20,'cax',cax,'title',ctitle);
    plotCond(JTdiff,JTdiff.grid,OPTIONS);    
end
%%  plot comparison
plotBands = [1 4 7 10];
for k = 1:length(plotBands)
   T = d0{plotBands(k)}.T;
   figure('Position',[100,100,400,700])
   axes('Position',[.1,.55,.8,.35])
   plot(real(Jxdm{plotBands(k)}.Z),'ks')
   hold on
   if (SENS_MATRIX)
     plot(f*real(JxdmComp{plotBands(k)}.Z),'bo')
   end
   plot(real(d{plotBands(k)}.Z),'r+')
   title(['Period = ', num2str(T)])
   axes('Position',[.1,.1,.8,.35])
   plot(imag(Jxdm{plotBands(k)}.Z),'ks')
   hold on
   if (SENS_MATRIX)
     plot(f*imag(JxdmComp{plotBands(k)}.Z),'bo')
   end
   plot(imag(d{plotBands(k)}.Z),'r+')
end
