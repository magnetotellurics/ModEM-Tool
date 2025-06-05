function [Z,loc,T,sites,origin,lat,lon] = allData2Imp(allData,Mplot)

%  translates allData structure (read in from ModEM ascii or
%   binary file with, e.g., readZ_3D into arrays that were
%  returned by rdZ (used for reading binary files output by
%  old version of forward code, written by Kush Tandon
%  THIS WILL ONLY WORK IF ALL PERIODS HAVE THE SAME NUMBER OF COMPONENTS

nPer = length(allData);
T = zeros(nPer,1);
Loc = allData{1}.siteLoc;
Sites = allData{1}.siteChar;
lat = allData{1}.lat;
lon = allData{1}.lon;
for j = 2:nPer
    %Loc = union(Loc,allData{j}.siteLoc,'rows');
    [temp,ind] = setdiff(allData{j}.siteLoc,Loc,'rows');
    Loc = [Loc; allData{j}.siteLoc(ind,:)];
    Sites = [Sites; allData{j}.siteChar(ind,:)];
    lat = [lat; allData{j}.lat(ind)];
    lon = [lon; allData{j}.lon(ind)];
end
nSites = size(Sites,1);
nResp  = size(allData{1}.Z,2);
nRow = nResp/2;
Z = NaN * zeros(nRow,2,nSites,nPer) + NaN*1i * zeros(nRow,2,nSites,nPer);
for j = 1:nPer
   T(j) = allData{j}.T;
   for k = 1:nSites
       [exists,ind] = ismember(Loc(k,1:2),allData{j}.siteLoc(:,1:2),'rows');
       if exists
           Z(1,1,k,j) = allData{j}.Z(ind,1);
           Z(1,2,k,j) = allData{j}.Z(ind,2);
           Z(2,1,k,j) = allData{j}.Z(ind,3);
           Z(2,2,k,j) = allData{j}.Z(ind,4);
           if nRow ==   3
               Z(3,1,k,j) = allData{j}.Z(ind,5);
               Z(3,2,k,j) = allData{j}.Z(ind,6);
           end
       end
   end
end
loc = Loc';
sites = Sites;
if isfield(allData{1},'origin') % use degrees, not km
    origin = allData{1}.origin;
    if Mplot
        disp(['Using degrees with origin ' num2str(origin)])
        [loc(1,:),loc(2,:)] = xy2latlon(loc(1,:),loc(2,:),origin(1),origin(2));
    end
else
    origin = []; % cannot use degrees: no origin 
end
    
