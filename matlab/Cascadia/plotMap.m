load map
inds = find(abs(map(:,1))>90000);
map(inds,:) = NaN;
map(:,1) = 360-map(:,1);
load idaho.bln;
plot(map(:,1),map(:,2),'k','lineWidth',2)
hold on
plot(idaho(:,1),idaho(:,2),'k','linewidth',2)
