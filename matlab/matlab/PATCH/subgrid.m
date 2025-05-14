%  current grid is [n1_c:n2_c] x [m1_c:m2_c] with limits lims_c
%   use rubber band box to define new subgrid, limits
box_ind;
% find range for new subgrid to plot, to keep aspect ratio 1;
n_new = n2-n1;m_new = m2-m1;
n0 = (n1+n2)/2; m0 = (m1+m2)/2;
if( n_new/m_new > n/m)
%  increase m_new ...
   m_new = m*n_new/n;
   if( m0-m_new/2 < 1)
      m1 = 1;
      m2 = m1+fix(m_new)-1;
   elseif (m0+m_new/2 > m)
      m2 = m;
      m1 = m2-fix(m_new)+1;
   else
      m1 = fix(m0-m_new/2); 
      m2 = m1+fix(m_new)-1;
   end
else
   n_new = n*m_new/m;
   if( n0-n_new/2 < 1)
      n1 = 1;
      n2 = n1+fix(n_new)-1;
   elseif (n0+n_new/2 > n)
      n2 = n;
      n1 = n2-fix(n_new)+1;
   else
      n1 = fix(n0-n_new/2); 
      n2 = n1+fix(n_new)-1;
   end
end

lims_c = [ lims_c(1) + (n1-n1_c)*dx lims_c(1) + (n2-n1_c+1)*dx ...
           lims_c(3) + (m1-m1_c)*dy lims_c(3) + (m2-m1_c+1)*dy ] ;

%x = [ lims(1) lims(2) lims(2) lims(1) lims(1) ];
%y = [ lims(3) lims(3) lims(4) lims(4) lims(3) ];
%hold on
%plot(x,y)
%hold off

m1_old = m1_c; n1_old = n1_c; m2_old = m2_c; n2_old = n2_c;

m1_c = m1; n1_c = n1; m2_c = m2; n2_c = n2;

pltmaskMT(hz(n1_c:n2_c,m1_c:m2_c),lims_c,gcf,mask(n1_c:n2_c,m1_c:m2_c),hLims)
set(gcbo,'Value',0);
hold on
plot(info{1}.lon,info{1}.lat','k*','Markersize',10,'linewidth',2)