prCond = readCond_3D('model/Yellowstone_10km_nested_200ohmm_pr.ws',2);
oldCond = readCond_3D('Yellowstone2_14km_ErrFl5_14freq_restart3_NLCG_064.rho',2);

newgrid = prCond.grid;

% replace NaN's, if any, with these values... use these for conductivity;
nz = length(newgrid.dz);
bg(1:nz) = 10^(-2.3);
bg(nz-1) = 10^(-2);
bg(nz) = 10^(-1.5);
newCond = interpCond_3D(newgrid,oldCond,bg);

% OR use zero for conductivity perturbations
% newCond = interpCond_3D(newgrid,oldCond,0);

writeCond_3D('model/Yellowstone_14km_ErrFl5_14freq_inv_10km_nested.ws',newCond,2);
