fname = 'cas_model.01';
saveFile = 'CascadiaInv.mat';

fid=fopen(fname,'r');
tline=fgetl(fid);
rms=sscanf(tline,'%*s%*s%*s%*s%*s%s',1);
nx=fscanf(fid,'%d',1);
ny=fscanf(fid,'%d',1);
nz=fscanf(fid,'%d',1);
fscanf(fid,'%d',1);
dx=fscanf(fid,'%f',nx);
dy=fscanf(fid,'%f',ny);
dz=fscanf(fid,'%f',nz);
rho = zeros(nx,ny,nz);
for i=1:nz
  tmp=fscanf(fid,'%f',nx*ny);
  rho(:,:,i)=reshape(tmp,[nx,ny,1]); % The WS model files run N->S,W->E,T->B
end

rho=log10(flipdim(rho,1)); % flipping around so that rho(:,1,1) is ordered south to north

fclose(fid);

fname = 'cas_resp.01';
fid=fopen(fname,'r');
ns=fscanf(fid,'%d',1);
nf=fscanf(fid,'%d',1);
nresp=fscanf(fid,'%d\n',1);
tline=fgetl(fid);
sx=fscanf(fid,'%f\n',ns); % read in N/S site locations
tline=fgetl(fid);
sy=fscanf(fid,'%f',ns); % read in E/W site locations
f=[]; Zp=[];
for i=1:nf
  f=[f fscanf(fid,'%*s%f',1)];
  temp=fscanf(fid,'%f',nresp*ns);
  Zp(i,:,:)=reshape(temp,nresp,ns);
end
fclose(fid);
% convert Zp to more standard phase convention (e^iwt = -1 * e^-iwt, ImZ --> -Im(Z))
if (nresp==8)
  Zp(:,[2 4 6 8],:)=-1.*Zp(:,[2 4 6 8],:);
else
  Zp(:,[2 4],:)=-1.*Zp(:,[2 4],:);
end

% Station locations are relative to mesh center

dx=dx./1000; dy=dy./1000; dz=dz./1000; sx=sx./1000; sy=sy./1000;
xoff=sum(dx(1:nx/2));
yoff=sum(dy(1:ny/2));

xb(1)=0-xoff;
for i=1:nx-1
  xb(i+1)=xb(i)+dx(i);
end
xb(nx+1)= xb(nx) + dx(nx);

yb(1)=0-yoff;
for i=1:ny-1
  yb(i+1)=yb(i)+dy(i);
end
yb(ny+1)= yb(ny) + dy(ny);
zb(1)=0;
for i=1:nz-1
  zb(i+1)=zb(i)+dz(i);
end
zb(nz+1)= zb(nz) + dz(nz);

%   load lat-lons of sites
load sitelocations
latAv = mean(sitelocations(:,1));
lonAv = mean(sitelocations(:,2));
lonAv = lonAv+360;
x0 = mean(sy);
y0 = mean(sx);
Origins = struct('x0',x0,'y0',y0,'lat0',latAv,'lon0',lonAv);
[s_lat,s_lon] = xformToLL(sy,sx,Origins);

X = xb'*ones(size(yb));
Y = ones(size(xb'))*yb;
[Yll,Xll]=xformToLL(Y,X,Origins);
[NX,NY] = size(Yll);
NZ = length(zb);
X = zeros(NX,NY,NZ);
Y = zeros(NX,NY,NZ);
Z = zeros(NX,NY,NZ);
for k = 1:NZ
  X(:,:,k) = Xll;
  Y(:,:,k) = Yll;
  Z(:,:,k) = zb(k);
end

rhoPad = zeros(nx+2,ny+2,nz+2);
rhoPad(2:end-1,2:end-1,2:end-1) = rho;

rhoPad(1,:,:) = rhoPad(2,:,:);
rhoPad(end,:,:) = rhoPad(end-1,:,:);

rhoPad(:,1,:) = rhoPad(:,2,:);
rhoPad(:,end,:) = rhoPad(:,end-1,:);

rhoPad(:,:,1) = rhoPad(:,:,2);
rhoPad(:,:,end) = rhoPad(:,:,end-1);

xW = [dx(1); dx ; dx(end)];
yW = [dy(1); dy ; dy(end)];
zW = [dz(1); dz ; dz(end)];

for k = 1:nx+1
   rhoPad(k,:,:) = (rhoPad(k,:,:)*xW(k)+ ...
	rhoPad(k+1,:,:)*xW(k+1))/(xW(k)+xW(k+1));
end

for k = 1:ny+1
   rhoPad(:,k,:) = (rhoPad(:,k,:)*yW(k)+ ...
	rhoPad(:,k+1,:)*yW(k+1))/(yW(k)+yW(k+1));
end

for k = 1:nz+1
   rhoPad(:,:,k) = (rhoPad(:,:,k)*zW(k)+ ...
	rhoPad(:,:,k+1)*zW(k+1))/(zW(k)+zW(k+1));
end

Res = struct('X',X,'Y',Y,'Z',Z,'lat',Yll,'lon',Xll,'z',zb,...
   'log10_rho',rhoPad(1:end-1,1:end-1,1:end-1),'sLat',s_lat,'sLon',s_lon)
eval(['save ' saveFile ' Res']);

