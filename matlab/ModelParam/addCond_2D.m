function [m] = addCond_2D(m1,m2,LOGCOND)

% Usage: [m] = addCond_2D(m1,m2,LOGCOND)
% 
% For two input model parameter, outputs their sum.

if nargin < 2
    error('Two input model parameters required')
elseif nargin < 3
    LOGCOND = 1;
end

% convert both models to LOGCOND
if (LOGCOND)
    if strcmp(m1.paramType,'LINEAR')
        v1 = log(m1.v);
    else
        v1 = m1.v;
    end
    if strcmp(m2.paramType,'LINEAR')
        v2 = log(m2.v);
    else
        v2 = m2.v;
    end
else
    if strcmp(m1.paramType,'LOGE')
        v1 = exp(m1.v);
    else
        v1 = m1.v;
    end
    if strcmp(m2.paramType,'LOGE')
        v2 = exp(m2.v);
    else
        v2 = m2.v;
    end
end

% if grid exists, check for consistency and save
if isfield(m1,'grid') && isfield(m2,'grid')
    if (m1.grid.Nz ~= m2.grid.Nz) || (m1.grid.Ny ~= m2.grid.Ny)
        error('Grids not consistent in addCond')
    end
end
        
% create the output
if isfield(m1,'grid')
    m.grid = m1.grid;
elseif isfield(m2,'grid')
    m.grid = m2.grid;
end
m.AirCond = m1.AirCond;
m.v = v1 + v2;
if (LOGCOND)
    m.paramType = 'LOGE';
else
    m.paramType = 'LINEAR';
end