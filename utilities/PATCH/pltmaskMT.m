function pltbathMT(hz,lim,hfig,mask,hLims)
%   Usage: pltbath(hz,lim,hfig,mask);
%   FIRST NEED TO CREATE FIGURE WITH HANDLE hfig !
nbath = 32;
figure(hfig);
lon1=lim(1);lon2=lim(2); 
lat1=lim(3);lat2=lim(4); 
[n,m]=size(hz) ;
hRange = hLims(2)-hLims(1);
hz=ceil((hz-hLims(1))*nbath/hRange);
hz = min(hz,32);
hz = hz.*(mask == 0) + mask;
hz = hz';
stlon=(lon2-lon1)/(n-1); 
stlat=(lat2-lat1)/(m-1); 
lon=[lon1:stlon:lon2]; 
lat=[lat1:stlat:lat2];
hold off;
image(lon,lat,hz);
%set(gco,'Tag','GridImage')
axis('xy');
set(gca,'Fontsize',12,'FontWeight','bold','Tag','GridPlotAxis');
