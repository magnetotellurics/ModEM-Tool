function [k,name] = whichPatch_ll(lat,lon,patchStruct)
%   returns number and name of patch which point with coordinates (lon,lat)
%       is in.   k = 0 if  no patch contains this point.
%   Usage:  [k,name] = whichPatch_ll(lat,lon,patchStruct)

MASK = zeros(size(patchStruct.mask));
for k = 1:patchStruct.nPatches
    MASK(patchStruct.mask==patchStruct.mask_nums(k)) = k;
end
[n,m] = size(MASK);

dLat = (patchStruct.ll_lims(4)-patchStruct.ll_lims(3))/m;
dLon = (patchStruct.ll_lims(2)-patchStruct.ll_lims(1))/n;

ii = ceil((lon-patchStruct.ll_lims(1))/dLon);
jj = ceil((lat-patchStruct.ll_lims(3))/dLat);
ii = max(ii,1);
jj = max(jj,1);
k = MASK(ii,jj);
if k >0
    name = patchStruct.names{k};
else
    name = '';
end