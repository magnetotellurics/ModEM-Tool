%  box_ind  using rubber-band box, find enclosed range of grid indices
%   assumes n1_c, n2_c lims_c are defined
k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');
finalRect = rbbox;
point2 = get(gca,'CurrentPoint');
point1 = point1(1,1:2);
point2 = point2(1,1:2);
p1 = min(point1,point2);
offset = abs(point2-point1);
n_c = n2_c-n1_c+1;
m_c = m2_c-m1_c+1;
dx = (lims_c(2)-lims_c(1))/(n_c-1);
dy = (lims_c(4)-lims_c(3))/(m_c-1);
n2 = n1_c + round((p1(1) + offset(1) - lims_c(1))/dx);
n1 = n1_c + round((p1(1) - lims_c(1))/dx);
m2 = m1_c + round((p1(2) + offset(2) - lims_c(3))/dy);
m1 = m1_c + round((p1(2) - lims_c(3))/dy);

n2 = min(n2,n);
m2 = min(m2,m);
n1 = max(n1,1);
m1 = max(m1,1);
