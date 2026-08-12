fname = findfiles('xml','SPUD_bundle_2020-03-04T21.07.49');
n = length(fname);
for i = 1:n
    tflist{i} = mttf.read(fname{i},'xml','Full_Impedance');
end
save('USA_impedance_list_q0+_20200304.mat','tflist');

i=1;
for j = 1:length(tflist)
    if str2double(tflist{j}.metadata.QUALITYRATING) >= 3
        tflist3{i} = tflist{j}; i = i+1;
    end
end
save('USA_impedance_list_q3+_20200304.mat','tflist3');
        
i=1;
for j = 1:length(tflist)
    if str2double(tflist{j}.metadata.QUALITYRATING) >= 4
        tflist4{i} = tflist{j}; i = i+1;
    end
end
save('USA_impedance_list_q4+_20200304.mat','tflist4');

%%
fname = findfiles('xml','SPUD_bundle_2020-03-04T21.07.49');
n = length(fname);
for i = 1:n
    tflist{i} = mttf.read(fname{i},'xml','Full_Vertical_Components');
end
save('USA_tipper_list_q0+_20200304.mat','tflist');

i=1;
for j = 1:length(tflist)
    if str2double(tflist{j}.metadata.QUALITYRATING) >= 3
        tflist3{i} = tflist{j}; i = i+1;
    end
end
save('USA_tipper_list_q3+_20200304.mat','tflist3');
        
i=1;
for j = 1:length(tflist)
    if str2double(tflist{j}.metadata.QUALITYRATING) >= 4
        tflist4{i} = tflist{j}; i = i+1;
    end
end
save('USA_tipper_list_q4+_20200304.mat','tflist4');

%% set siteChar to tfname - we don't always want this except for SPUD loaded diverse data sets
for j = 1:length(tflist)
    tmp = strsplit([tflist{j}.tfname tflist{j}.tfext],'/'); 
    tfname = tmp{end}; tfname(1:end-4)
    tflist{j}.siteChar = tfname(1:end-4);
end 

%% 
usarraysample = mttf.read('/Users/akelbert/Developer/EMTF-FCU/original/USGS/Florida/FL001_merged.zmm');
T = usarraysample.T;
T_longtail = logspace(log10(T(end)),log10(300000),7);
T = [T T_longtail(2:end)];
load('USA_impedance_list_q0+_20200304.mat');
imp = mtdata.tflist2mtdata(tflist,'Full_Impedance',T,'latlon');
save('USA_impedance_20200304_45per.mat','imp');
imp = imp.fillNaNs(0.0);
imp.write('USA_impedance_20200304_45per.dat','list');
load('USA_tipper_list_q0+_20200304.mat');
tip = mtdata.tflist2mtdata(tflist,'Full_Vertical_Components',T,'latlon');
save('USA_tipper_20200304_45per.mat','tip');
tip = tip.fillNaNs(0.0);
tip.write('USA_tipper_20200304_45per.dat','list');

%% 
imp = mtdata.read('data/USA_impedance_20200304_45per.dat','list','[mv/km]/[nT]','Full_Impedance',0);
imp25 = mtdata.read('USA_impedance_20200304_45per_0.25x0.25.dat','list','[mv/km]/[nT]','Full_Impedance',1);
impcomp = mtdatacompare(imp,imp25);

%%
sitelatlon = [dat.d{1}.lat dat.d{1}.lon zeros(dat.d{1}.nSite,1)]';
sitexy = dat.d{1}.siteLoc';
newdat = newdat.initSites(dat.d{1}.lat,dat.d{1}.lon,dat.d{1}.siteChar,dat.origin(1),dat.origin(2),0);
newdat.header = 'USA gridded impedances at NOAA locations for e.g., USA2x2_cartesian_regional_lat49.rho';
newdat.predicted = 1;
newdat.write('USA_NOAA_lat49_45per_template.dat','list');


%%
cd /Users/akelbert/Developer/EMTF-FCU/database/xml/; 
clear tflist tflist3 tflist4 fname;

fname = findfiles('xml','USArray_TA'); n1 = length(fname);
for i = 1:n1
    j = i;
    tflist{j} = mttf.read(fname{i},'xml','Full_Vertical_Components');
end
fname = findfiles('xml','USArray_BB'); n2 = length(fname);
for i = 1:n2
    j = n1 + i;
    tflist{j} = mttf.read(fname{i},'xml','Full_Vertical_Components');
end
fname = findfiles('xml','USGS/Florida'); n3 = length(fname);
for i = 1:n3
    j = n1 + n2 + i;
    tflist{j} = mttf.read(fname{i},'xml','Full_Vertical_Components');
end
save('USArray_tipper_list_all_20170919.mat','tflist');

i=1;
for j = 1:length(tflist)
    if str2double(tflist{j}.metadata.QUALITYRATING) >= 3
        tflist3{i} = tflist{j}; i = i+1;
    end
end
save('USArray_tipper_list_q3+_20170919.mat','tflist3');
        
i=1;
for j = 1:length(tflist)
    if str2double(tflist{j}.metadata.QUALITYRATING) >= 4
        tflist4{i} = tflist{j}; i = i+1;
    end
end
save('USArray_tipper_list_q4+_20170919.mat','tflist4');

%%
load('USArray_impedance_list_q3+_20170919.mat');
impedance = mtdata('Full_Impedance');
for j = 1:length(tflist3)
    tf = tflist{j};
    impedance = impedance.set(tf);
end
impedance = impedance.regroup;
save('USArray_impedance_q3+_20170919.mat','impedance');

%%
% NEED TO MAKE ADDITIONAL CHANGES TO SAVE THE MTTF QUALITY INFORMATION 

% FIX BUGS IN INITIAL READ
% impedance = impedance.changeUnits('[V/m]/[T]');
% impedance.units = '[mV/km]/[nT]';
% for i=1:obj.nPeriods
%     impedance.d{i}.units = '[mV/km]/[nT]';
% end
% 
% impedance.units = '[V/m]/[T]';
% for i=1:impedance.nPeriods
%     impedance.d{i}.units = '[V/m]/[T]';
% end
% impedance = impedance.changeUnits('[mV/km]/[nT]');
% 
% for i=1:impedance.nPeriods
%     impedance.d{i}.lat = impedance.d{i}.lat';
%     impedance.d{i}.lon = impedance.d{i}.lon';
% end
% impedance = impedance.regroup;
% 
% for i=1:impedance.nPeriods
%     impedance.d{i}.nSite = length(impedance.d{i}.lat);
% end
%     
% END FIX BUGS IN INITIAL READ

fname = findfiles('xml','USArray_BB');
for i = 1:length(fname)
    tf = mttf.read(fname{i},'xml','Full_Impedance');
    impedance = impedance.set(tf);
end
impedance = impedance.regroup;

save('USArray_impedance_20170228.mat','impedance');

lims = struct('lonmin',-94,'lonmax',-92,'latmin',43,'latmax',45);
bin44_93 = impedance.subset(lims);  obj = bin44_93;
obj = obj.changeUnits('[V/m]/[T]');
obj.units = '[mV/km]/[nT]';
for i=1:obj.nPeriods
    obj.d{i}.units = '[mV/km]/[nT]';
end
obj = obj.regroup;

tfpred = mttf.read();
%tfpred = obj3.get([44 -93]); 
f = tf1d.apresplt([1 2e5]); f2 = tf1d.impplt;
[allsites,allsitesloc,lat,lon] = obj.getSites;
nSites = length(lat);
tf = obj.get(allsites(2,:)); 
%f = tf.apresplt([1 2e5]); f2 = tf.impplt;
tf.apresplt(f); tf.impplt(1,f2);
for i = 2:nSites
    tf = obj.get(allsites(i,:)); tf.apresplt(f); tf.impplt(1,f2);
end

%%
cd /Users/akelbert/Developer/EMTF-FCU/original/;
fname = findfiles('edi','1D_impedances_Fernberg_models');

fernberg1d = mtdata('Full_Impedance');
for i = 1:length(fname)
    tf = mttf.read(fname{i},'edi','Full_Impedance');
    fernberg1d = fernberg1d.set(tf);
end
fernberg1d = fernberg1d.regroup;
fernberg1d.predicted = 1;
tf1d = fernberg1d.get('CP1');
tf1d.apresplt([1 2e5])

%%
PBfile = '1D_impedances_Fernberg_models/Benchmark_region.txt'; clear tmp tmp2;
fid = fopen(PBfile,'r');
tmp = textscan(fid,'%s','HeaderLines',1,'Delimiter','\n');
fclose(fid);
tmp = tmp{1};

for i = 1:length(tmp)
    line{i} = strsplit(char(tmp(i)));
    site.name{i} = char(line{i}(1));
    site.lat(i) = str2num(char(line{i}(2)));
    site.lon(i) = str2num(char(line{i}(3)));
    site.region{i} = char(line{i}(4));
    site.regionid(i) = str2num(char(line{i}(5)));
    %sitename{i} = [sprintf('%03d',i) '-' char(tmp2(i,5))];
end

%% 
region = 'CP1'; % only 5 sites in Paul's file
region = 'PT1'; % missing from Paul's file
region = 'IP1'; 
ireg = find(strncmp(region,site.region,3));
sitelist = char(site.name{:});
regsites = impedance.subset(sitelist(ireg,:));

selectsites = sitelist(ireg,:); nSites = size(selectsites,1);
tf = regsites.get(selectsites(1,:));
[f,p] = tf.apresplt([1 2e5]); [f2,p2] = tf.impplt;
set(p(1),'Color',[0    0.4000    1.0000]); hold on;
set(p(2),'Color',[1.0000    0.6000    0.6000]); hold on;
set(p(3),'Color',[0    0.4000    1.0000]); hold on;
set(p(4),'Color',[1.0000    0.6000    0.6000]); hold on;
for i = 2:nSites
    tf = regsites.get(selectsites(i,:)); [~,p] = tf.apresplt(f); [~,p2] = tf.impplt(1,f2);
    set(p(1),'Color',[0    0.4000    1.0000]); hold on;
    set(p(2),'Color',[1.0000    0.6000    0.6000]); hold on;
    set(p(3),'Color',[0    0.4000    1.0000]); hold on;
    set(p(4),'Color',[1.0000    0.6000    0.6000]); hold on;
end

tf1d = fernberg1d.get(region);
[~,p1d] = tf1d.apresplt(f); [~,p1d2] = tf1d.impplt(1,f2);
set(p1d,'Color','k','Linewidth',2); hold on;
[~,p] = tf.apresplt(f); [~,p2] = tf.impplt(1,f2);
set(p(1),'Color',[0    0.4000    1.0000]); hold on;
set(p(2),'Color',[1.0000    0.6000    0.6000]); hold on;
set(p(3),'Color',[0    0.4000    1.0000]); hold on;
set(p(4),'Color',[1.0000    0.6000    0.6000]); hold on;
set(f.Children(3).Title,'String',region);
set(f.Children(1).Title,'String','');

%%
% 80.4672 km = 50 miles
% 128.748 km = 80 miles
% 160.934 km = 100 miles
WashingtonDC = [38.904722 -77.016389]; region1 = 'CP1'; region2 = 'PT1';
%Minneapolis = [44.983333 -93.266667]; region1 = 'SU1'; region2 = 'IP1'; % skip i=2
%Chicago = [41.836944 -87.684722]; region1 = 'IP3'; region2 = 'IP1';
[tflist,ind] = impedance.get(WashingtonDC,160.934); 
% for  now, a roundabout way to populate the metadata quality info
clear tflist; i = 1; qualityrating = [];
for j = ind
    tflist{i} = mttf.read(fname{j},'xml','Full_Impedance'); i = i+1;
end
for i = 1:length(tflist)
    qualityrating = [qualityrating str2num(tflist{i}.metadata.QUALITYRATING)]; 
end
% end the roundabout way to populate the metadata quality info!
[f,p] = tflist{1}.apresplt([1 2e5],[1 1e4]);  [f2,p2] = tflist{1}.impplt(1);
set(p(1),'Color',[0    0.4000    1.0000]); hold on;
set(p(2),'Color',[1.0000    0.6000    0.6000]); hold on;
set(p(3),'Color',[0    0.4000    1.0000]); hold on;
set(p(4),'Color',[1.0000    0.6000    0.6000]); hold on;
%set(p2(2*(1:4)),'Color','w') % only real parts plotted
for i = 2:length(tflist)
    if qualityrating(i) > 3
    [~,p] = tflist{i}.apresplt(f);  [~,p2] = tflist{i}.impplt(1,f2);
    set(p(1),'Color',[0    0.4000    1.0000]); hold on;
    set(p(2),'Color',[1.0000    0.6000    0.6000]); hold on;
    set(p(3),'Color',[0    0.4000    1.0000]); hold on;
    set(p(4),'Color',[1.0000    0.6000    0.6000]); hold on;
    %set(p2(2*(1:4)),'Color','w') % only real parts plotted
    end
end
tf1d = fernberg1d.get(region1);
[~,p1d] = tf1d.apresplt(f);  [~,p1d2] = tf1d.impplt(1,f2);
set(p1d,'Color',[0 0.5 0],'Linewidth',2); hold on;
set(p1d2,'Color',[0 0.5 0],'Linewidth',2); hold on;
tf1d = fernberg1d.get(region2);
[~,p1d] = tf1d.apresplt(f);  [~,p1d2] = tf1d.impplt(1,f2);
set(p1d,'Color',[0.4 0.2 0],'Linewidth',2); hold on;
set(p1d2,'Color',[0.4 0.2 0],'Linewidth',2); hold on;
[~,p] = tflist{1}.apresplt(f); [~,p2] = tflist{1}.impplt(1,f2);
set(p(1),'Color',[0    0.4000    1.0000]); hold on;
set(p(2),'Color',[1.0000    0.6000    0.6000]); hold on;
set(p(3),'Color',[0    0.4000    1.0000]); hold on;
set(p(4),'Color',[1.0000    0.6000    0.6000]); hold on;
%set(f.Children(3).Title,'String','100 miles around Chicago');
set(f.Children(3).Title,'String','100 miles around Washington DC');
set(f2.Children(2).Title,'String','');
set(f2.Children(4).Title,'String','');
set(f2.Children(6).Title,'String','');
set(f2.Children(8).Title,'String','');
str = impedance.v.code(ind,:); str(qualityrating>3,:)

%%
tipper = mtdata('Full_Vertical_Components');
for i = 1:4%length(fname)
    tf = mttf.read(fname{i},'xml','Full_Vertical_Components');
    tipper = tipper.set(tf);
end
tipper = tipper.regroup;

% testing: write_xml(tf,'test.xml','TEST'); impplt(tf,2);
