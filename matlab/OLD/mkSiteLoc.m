siteLoc = [y(50) (y(50)+y(51)/2 y(51) (2*y(51)+y(52))/3 y(52)]
siteLoc = [siteLoc' z(11)*ones(5,1)] 
siteLocCenter = siteLoc - ones(5,1)*[347500      444440]
plot(siteLocCenter(:,1), siteLocCenter(:,2),'w*')
cfile = 'TestSites1.loc';
fid = fopen(cfile,'w');
[nSites,dum] = size(siteLoc)
fprintf(fid,'%i\n',nSites)
fprintf(fid,'%12.4f %12.4f\n',siteLoc')
fclose(fid)
