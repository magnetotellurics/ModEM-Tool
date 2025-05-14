function sigDepth(obj.latLims,lonLims,istep,zLims)

ii = find(obj.llgrid.lat>latLims(1) & obj.llgrid.lat<latLims(2));
jj = find(obj.llgrid.lon>lonLims(1) & obj.llgrid.lon<lonLims(2));
i1 = ii(1);i2 = ii(end);
j1 = jj(1);j2 = jj(end);
for i =i1:1.5*istep:i2
    figure
    %for j = j1:istep:j2
        plot(squeeze(obj.CONDxyz(i,j1:istep:j2,:))',squeeze(obj.Z(i,j1:istep:j2,:))');
        hold on
    %end
    set(gca,'ydir','reverse','fontweight','demi','fontsize',14)
    fatlines(gca,2);
    title([ 'Latitude = ' num2str(obj.llgrid.lat(i))])
    ylim([0,300]);
    eval(['print -dpdf sigDepth' num2str(i)])
end

