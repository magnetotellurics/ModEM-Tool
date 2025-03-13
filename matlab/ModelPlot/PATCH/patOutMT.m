%  writes a patch file in matlab
%  mimicing a fortran binary sequential file
% USAGE:   pat_out(cfile,npatch,patch_names,patch_mask,mask);

function   pat_outMT(cfile,npatch,patch_names,patch_nums,mask,ll_lims);

eval(['save ' cfile 'npaatch,patch_names,patch_nums,mask,ll_lims'])
end
