function [d] = addNoise_3D(d1,relErr)
%  Usage : [d] = addNoise_3D(d1,relErr);
%
% adds "realistic" percentage errors, unique for each data point;
% assumes off diagonal impedance, full impedance or impedance plus Hz.

d = d1;
nTx = length(d);
for k = 1:nTx
    temp = d1{k}.Z;
    [nsite,ncomp]  = size(temp);
    if(ncomp == 2) 
        stderr = sqrt(abs(temp(:,1).*conj(temp(:,2))));
    else
        stderr = sqrt(abs(temp(:,2).*conj(temp(:,3))));
    end
    errScales = repmat(relErr*stderr,1,ncomp);
    if(ncomp == 6) % assume that 5&6 are vertical transfer functions...
       errScales(:,5:6) = relErr*ones(nsite,2);
    end
    d{k}.Z = d1{k}.Z + ...
	errScales.*(randn(size(d1{k}.Z))+1i*randn(size(d1{k}.Z)));
    d{k}.Zerr = errScales;
end

% d = d1;
% nTx = length(d);
% for k = 1:nTx
%     temp = d1{k}.Z;
%     [nsite,ncomp]  = size(temp);
%     if(ncomp == 2) 
%        errScale = nansum(nansum(temp.*conj(temp)));
%        errScale = relErr*sqrt(errScale/(4*nsite));
%        errScales = errScale*ones(nsite,ncomp);
%     elseif(ncomp == 4)
%        errScale = nansum(nansum(temp(:,2:3).*conj(temp(:,2:3))));
%        errScale = relErr*sqrt(errScale/(4*nsite));
%        errScales = errScale*ones(nsite,ncomp);
%     elseif(ncomp == 6)
%        errScale = nansum(nansum(temp(:,2:3).*conj(temp(:,2:3))));
%        errScale = relErr*sqrt(errScale/(4*nsite));
%        errScales = errScale*ones(nsite,4);
%        errScales = [errScales relErr*ones(nsite,2)];
%     end
%     d{k}.Z = d1{k}.Z + ...
% 	errScales.*(randn(size(d1{k}.Z))+1i*randn(size(d1{k}.Z)));
%     d{k}.Zerr = errScales;
% end
