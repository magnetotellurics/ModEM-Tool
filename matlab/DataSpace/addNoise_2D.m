function [d] = addNoise_2D(d1,relErr)
%  Usage : [d] = addNoise_2D(d1,relErr);

d = d1;
nTx = length(d);
for k = 1:nTx
    temp = d1{k}.Z;
    %   make error scale different for each site ... more realistic!
    %errScale = relErr*sqrt(real(temp'*temp)/length(temp));
    %    d{k}.Z = d1{k}.Z + errScale*(randn(size(d1{k}.Z))+...
	%                                 1i*randn(size(d1{k}.Z)));
    % d{k}.Zerr = errScale*ones(size(d1{k}.Z));
    errScale = relErr*abs(temp);
    d{k}.Z = d1{k}.Z + errScale.*(randn(size(d1{k}.Z))+...
	                                 1i*randn(size(d1{k}.Z)));
    d{k}.Zerr = errScale;
end