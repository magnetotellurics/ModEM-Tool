clear all;
rname='cascad'
fname = 'cas_model.01';
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
for i=1:nz
  tmp=fscanf(fid,'%f',nx*ny);
  rho(:,:,i)=reshape(tmp,nx,ny); % The WS model files run N->S,W->E,T->B
end

rho=log10(flipdim(rho,1)); % flipping around so that rho(:,1,1) is ordered south to north

rho(nx+1,:,:)=rho(nx,:,:);
rho(:,ny+1,:)=rho(:,ny,:);
rho(:,:,nz+1)=rho(:,:,nz);

fclose(fid);

% Read the stations locations and model responses

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


%Plotting

close all;        
h=figure;
set(h,'units','normalized','position',[0.25 0.25 0.6 0.6])
colormap(flipud(jet(32)));
    
while 1  % plot slices until user quits

mdlslice=menu('MODEL PLOT','Slice in Z','Quit');  

switch (mdlslice)
    case 1
    n=nz; z=zb;
    [y,x]=ndgrid(xb,yb); 
   xlims=[-500 480]; ylims=[-440 440]; % for cascadia
 %xlims=[-800 800]; ylims=[-740 740]; % for cascadia new mesh
%    xlims=[-220 220]; ylims=[-280 280]; % for pilot
    xlab='Distance (km)';
    ylab='Distance (km)'; 

    modslice=rho;
  case 2
    break;
  end

  for i=1:n
    hold off
    pcolor(x,y,squeeze(modslice(:,:,i)))
    shading flat
    hold on
    switch (mdlslice)  

      case 1
         plot(sy,sx,'k^')
    end  
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-0.5 3.5]) % SP profile
    xlabel(xlab)
    ylabel(ylab)

    text(xlims(1)+0.02*(xlims(2)-xlims(1)),ylims(1)+0.05*abs(ylims(2)-ylims(1)),'')       
    text(xlims(2)-0.04*(xlims(2)-xlims(1)),ylims(1)+0.05*abs(ylims(2)-ylims(1)),'')
    if (abs(z(i)) <= 1)
      title(['MODEL SLICE: ',num2str(z(i)*1000),' m to ',num2str(z(i+1)*1000),' m','  ',rname])
  else
      title(['MODEL SLICE: ',num2str(z(i)),' km to ',num2str(z(i+1)),' km','  ',rname])
  end
    colorbar
%    print('-djpeg',strcat(strrep(rname,'model','slice_'),num2str(i)))
    print('-depsc2',strcat(strrep(rname,'model','slice_'),num2str(i)))
    pause
  end

end 