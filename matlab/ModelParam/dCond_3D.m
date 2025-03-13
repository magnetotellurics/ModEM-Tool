function dm = dCond_3D(m,delta)

% Usage: dm = dCond_3D(m,delta)
% 
% For an input model parameter, outputs a perturbed model.
% Input parameter delta is optional and assumed logarithmic.

if nargin < 2
    delta = .02;
end

% extract info from the input model parameter
if strcmp(m.paramType,'LOGE')
    LOGCOND = 1;
else
    LOGCOND = 0;
end

dm = m;
%  generate random perturbation to log conductivity
randLogCond = delta*randn(size(m.v));
%   zero out top few rows, to eliminate "Q" stuff
%   just zero out first half
% randLogCond(:,:,1:2) = 0;
if LOGCOND
    dm.v = randLogCond + m.v;
else
    dm.v = randLogCond + log(m.v);
    dm.v = exp(dm.v);
end
% subtract the original model to leave the perturbation
dm.v = dm.v - m.v;
%[status] = writeCond_3D(wModelFile,dm);

%% Plot perturbation using the following code:
% nZplot = 18;
% nYskip = 7;
% cax = [-3.5,-.5]  %   +log10(2);  Model 2
% ctitle = ['Relative Error = ' num2str(100*FracError) '%'];
% OPTIONS = struct('nYskip',nYskip,'nZplot',nZplot,...
%                 'cax',cax,'title',ctitle);
% plotCond(dm,dm.grid,OPTIONS);
% eval(['print -djpeg90 ' fullfile(modelPath,modelBaseName)]);
