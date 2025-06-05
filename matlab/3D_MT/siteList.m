fdir = '.';
filt = [fdir '/*.dat'];
[filename, pathname] = uigetfile(filt, 'Impedance File');
fdata = [pathname filename];
fresp = [pathname filename];
filename = filename(1:end-4);
[rms,info] = DataFit(fresp,fdata);

% save site list in *.csv
fid = fopen('sites.txt','w');
for i = 1:length(info.code)
    if findstr(info.code(i,:),'SR')
        survey = 1;
    else
        survey = 0;
    end
    fprintf(fid,'%s\t%f\t%f\t%f\n',info.code(i,:),info.lat(i),info.lon(i),survey);
end
fclose(fid);
