    bndry_lon=[-118 -108 -108 -118 -118];
    bndry_lat=[46.5  46.5  41   41 46.5];
    
    clf;
    lonlim = [-126 -70];
    latlim = [24 50];
    %m_proj('miller','long',lonlim,'lat',latlim,'rectbox','on');
    m_proj('lambert','long',lonlim,'lat',latlim);

    m_gshhs_i('color','k','linewi',4);   % Coastline...
    m_gshhs_i('speckle','color','k');    % with speckle added

    m_line(bndry_lon,bndry_lat,'linewi',6,'color','r');     % Area outline ...
    m_hatch(bndry_lon,bndry_lat,'single',30,5,'color','r'); % ...with hatching added.


    m_grid('box','fancy','tickdir','in','fontsize',16);

    %m_grid('linewi',2,'linest','none','tickdir','out','fontsize',12);
    %title('USA','fontsize',14);
    %m_text(-120,48,5,{'SRPY'},'fontsize',18);
    
%     states = m_shaperead('usastatehi',...
%         'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
%     lat = [states.LabelLat];
%     lon = [states.LabelLon];
%     tf = ingeoquad(lat, lon, latlim, lonlim);
    states = m_shaperead('state_bounds');
    for i = 1:length(states.ncst)
        stateslon = states.ncst{i}(:,1);
        stateslat = states.ncst{i}(:,2);
        m_line(stateslon,stateslat,'linewi',2,'color','k');
    end
    %%
    
    for i = 1:length(states)
        switch char(states(i).Name)
            case 'Washington'
                states(i).Abbrev = 'WA';
            case 'Oregon'
                states(i).Abbrev = 'OR';
            case 'California'
                states(i).Abbrev = 'CA';
            case 'Idaho'
                states(i).Abbrev = 'ID';
            case 'Nevada'
                states(i).Abbrev = 'NV';
            case 'Utah'
                states(i).Abbrev = 'UT';
            case 'Arizona'
                states(i).Abbrev = 'AZ';
            case 'Montana'
                states(i).Abbrev = 'MT';
            case 'Wyoming'
                states(i).Abbrev = 'WY';
            case 'Colorado'
                states(i).Abbrev = 'CO';
            case 'New Mexico'
                states(i).Abbrev = 'NM';
            case 'North Dakota'
                states(i).Abbrev = 'ND';
            case 'South Dakota'
                states(i).Abbrev = 'SD';
            case 'Nebraska'
                states(i).Abbrev = 'NE';
            case 'Kansas'
                states(i).Abbrev = 'KS';
            case 'Oklahoma'
                states(i).Abbrev = 'OK';
            case 'Texas'
                states(i).Abbrev = 'TX';
            case 'Minnesota'
                states(i).Abbrev = 'MN';
            case 'Iowa'
                states(i).Abbrev = 'IA';
            case 'Missouri'
                states(i).Abbrev = 'MO';
            case 'Arkansas'
                states(i).Abbrev = 'AR';
            case 'Louisiana'
                states(i).Abbrev = 'LA';
            case 'Wisconsin'
                states(i).Abbrev = 'WI';
            case 'Illinois'
                states(i).Abbrev = 'IL';
            case 'Indiana'
                states(i).Abbrev = 'IN';
            case 'Kentucky'
                states(i).Abbrev = 'KY';
            case 'Tennessee'
                states(i).Abbrev = 'TN';
            case 'Mississippi'
                states(i).Abbrev = 'MS';
            case 'Alabama'
                states(i).Abbrev = 'AL';
            case 'Ohio'
                states(i).Abbrev = 'OH';
            case 'West Virginia'
                states(i).Abbrev = 'WV';
            case 'New York'
                states(i).Abbrev = 'NY';
            case 'Pennsylvania'
                states(i).Abbrev = 'PA';
            case 'Virginia'
                states(i).Abbrev = 'VA';
            case 'North Carolina'
                states(i).Abbrev = 'NC';
            case 'South Carolina'
                states(i).Abbrev = 'SC';
            case 'Georgia'
                states(i).Abbrev = 'GA';
            case 'Florida'
                states(i).Abbrev = 'FL';
           otherwise
                states(i).Abbrev = '  ';
        end
    end
    
    txt=m_text(lon(tf), lat(tf), {states(tf).Abbrev},'fontsize',14,...
    'HorizontalAlignment', 'center');
    
%     figure; ax = usamap({'CA','MT'});
% set(ax, 'Visible', 'off')
% latlim = getm(ax, 'MapLatLimit');
% lonlim = getm(ax, 'MapLonLimit');
% states = shaperead('usastatehi',...
%         'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
% geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
% 
% lat = [states.LabelLat];
% lon = [states.LabelLon];
% tf = ingeoquad(lat, lon, latlim, lonlim);
% textm(lat(tf), lon(tf), {states(tf).Name}, ...
%    'HorizontalAlignment', 'center')