function [z] = interpSlab(slabs,depths,lat,lon)
% usage: [z] = interpSlab(slabs,depths,lat,lon);
nDepth = length(depths);
for k = nDepth:-1:1
    i1 = find(slabs{k}(:,2)-lat>0, 1, 'last' );
    w2 = (slabs{k}(i1,2)-lat)/(slabs{k}(i1,2)-slabs{k}(i1+1,2)); 
    w1 = 1-w2;
    lonI = w2*slabs{k}(i1+1,1)+w1*slabs{k}(i1,1);
    if lon>lonI
        break
    else
        lon0 = lonI;
    end
end

if k<nDepth
    w1 = (lon0-lon)/(lon0-lonI);
    z = depths(k)*w1+depths(k+1)*(1-w1);
else
    %  need to extrapolate when k = nDepth
    lon0 = lonI;
    k = nDepth-1;
    i1 = find(slabs{k}(:,2)-lat>0, 1, 'last' );
    w2 = (slabs{k}(i1,2)-lat)/(slabs{k}(i1,2)-slabs{k}(i1+1,2)); 
    w1 = 1-w2;
    lonI = w2*slabs{k}(i1+1,1)+w1*slabs{k}(i1,1);
    w1 = (lon0-lon)/(lon0-lonI);
    z = depths(k)*w1+depths(k+1)*(1-w1);
    z = NaN;
end
