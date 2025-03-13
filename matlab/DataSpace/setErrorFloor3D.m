function [dOut] = setErrorFloor3D(dIn,relErr,site)
%  Usage : [dOut] = setErrorFloor3D(dIn,relErr,site);
%
%  If relErr has two components, the second is used for Tx & Ty

if length(relErr)==1
    relErr(2) = relErr(1)
end

dOut = dIn;
nTx = length(dIn);
for k = 1:nTx
    data = dIn{k}.Z;
    err = dIn{k}.Zerr;
    sites = dIn{k}.siteChar;
    [nsite,ncomp]  = size(data);
    if nargin > 2
        [str,i] = intersect(sites,site,'rows');
    else
        i = 1:nsite;
    end
    if ~isempty(i)
        if(ncomp == 2)
            errFloor = relErr(1) * sqrt(abs(data(i,1).*data(i,2)));
            errZ = max(err(i,1:2),repmat(errFloor,1,2));
            dOut{k}.Zerr(i,1:2) = errZ;
        elseif(ncomp == 4)
            errFloor = relErr(1) * sqrt(abs(data(i,2).*data(i,3)));
            errZ = max(err(i,1:4),repmat(errFloor,1,4));
            dOut{k}.Zerr(i,1:4) = errZ;
        elseif(ncomp == 6)
            errFloor = relErr(1) * sqrt(abs(data(i,2).*data(i,3)));
            errZ = max(err(i,1:4),repmat(errFloor,1,4));
            dOut{k}.Zerr(i,1:4) = errZ;
            errFloor = relErr(2);
            errTx = max(err(i,5),errFloor);
            dOut{k}.Zerr(i,5) = errTx;
            errFloor = relErr(2);
            errTy = max(err(i,6),errFloor);
            dOut{k}.Zerr(i,6) = errTy;
        end
    end
end
