%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run lisa_setprefix_ver2 first. ;
% run lisa_setnames_ver2 first. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% extract A_p. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fname_tmp = sprintf('%s/A_p.mda',dir_out_s0000); fcheck(fname_tmp);
A_p_ = mda_read_r8(fname_tmp);
fname_tmp = sprintf('%s/AZ_rsum.mda',dir_out_s0000); fcheck(fname_tmp);
AZ_rsum_ = mda_read_r8(fname_tmp);
