% plotStates
% Plot US state boundaries on an existing map with m_map
% To use, install m_map and add a setenv line to startup, e.g.
% setenv('MAPPATH','/home/mt/anya/ModEM/matlab/Toolbox/m_map/data/');
cbndry  = [0.9 0.9 0.9];
%m_gshhs_i('color',cbndry);              % Coastline... higher resolution
%m_gshhs_i('speckle','color',cbndry);    % with speckle added
m_coast('color',cbndry);              % Coastline...
m_coast('speckle','color',cbndry);    % with speckle added
states = m_shaperead('state_bounds');
for i = 1:length(states.ncst)
    stateslon = states.ncst{i}(:,1);
    stateslat = states.ncst{i}(:,2);
    m_line(stateslon,stateslat,'linewi',2,'color',cbndry);
end

% datadir = getenv('MAPPATH');            
% m_plotbndry([datadir 'washington'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'montana'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'oregon'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'idaho'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'wyoming'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'california'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'nevada'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'utah'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'colorado'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'nebraska'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'south_dakota'],'color',cbndry,'LineWidth',2)
% m_plotbndry([datadir 'north_dakota'],'color',cbndry,'LineWidth',2)