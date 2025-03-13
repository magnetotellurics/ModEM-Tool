%%
inFile = 'whiteNoise8.ws';
%Nx = 51;
%Ny = 51;
%Nz = 26;
Nx = 8; Ny = 8; Nz = 8;
NzEarth = Nz;
%dx = ones(Nx,1)*1000;
%dy = ones(Ny,1)*1000;
%dz = ones(Nz,1)*100;
dx = ones(Nx,1);
dy = ones(Ny,1);
dz = ones(Nz,1);
grid = struct('Nx',Nx,'Ny',Ny,'Nz',Nz,'dx',dx,'dy',dy,'dz',dz,...
    'NzAir',0,'NzEarth',NzEarth,'origin',[0, 0 , 0],'rotation',0,'units','m');
m_In = MT3DmodelParam(grid);
m_In.v(4,4,4) = 1;
mCov = MT3DmodelCovLaplace(grid);
%%
mOut = mCov.CovMult(m_In);
plotCond(mOut);
%%
m_In.v = randn(size(m_In.v));
m_In.writeVec(inFile)
%%
plotCond(m_In)
mCov.tol = 1e-2;
mOut = mCov.CovMult(m_In);
plotCond(mOut)
%mOut2 = mCov.CovMult(mOut);
%plotCond(mOut2)
%%   
mCovOld = MT3DmodelCov(grid);
mOutOld = mCovOld.CovMult(m_In);
plotCond(mOutOld)

%%
mCovBo = MT3DmodelCovPoisson(grid);
mOutBo = mCovBo.CovMult(m_In);
plotCond(mOutBo)
