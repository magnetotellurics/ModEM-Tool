fdir = '.';
filt = [fdir '/*.imp'];
[filename, pathname] = uigetfile(filt, 'Impedance File');
Zfile = [pathname filename];
[AllData] = readZ_3D(Zfile,'[mV/km]/[nT]');
[Imp,loc,T,sites,origin] = allData2Imp(AllData,0);
avgx = mean(loc(1,:));
avgy = mean(loc(2,:));
for k=1,length(AllData)
    x = AllData{k}.siteLoc(:,1);
    y = AllData{k}.siteLoc(:,2);
    origin = AllData{k}.origin;
    [lat,lon] = xy2latlon(x,y,origin(1),origin(2));
    Z = AllData{k}.Z;
    T = AllData{k}.T;
    figure;
    subplot(2,1,1);
    InterpPlot(lon,lat,real(Z(:,1)),'[mV/km]/[nT]',['Re(Zxx): ' num2str(T) ' secs']);
    subplot(2,1,2);
    InterpPlot(lon,lat,imag(Z(:,1)),'[mV/km]/[nT]',['Im(Zxx): ' num2str(T) ' secs']);
end
