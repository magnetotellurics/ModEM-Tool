[n,m] = size(hz);
n1_c = 1; n2_c = n; m1_c = 1; m2_c = m; lims_c = ll_lims;
m1_old = m1_c; n1_old = n1_c; m2_old = m2_c; n2_old = n2_c;
hfig = gcf;
pltmaskMT(hz,ll_lims,hfig,mask,hLims)
hold on
plot(info{1}.lon,info{1}.lat','k*','Markersize',10,'linewidth',2)