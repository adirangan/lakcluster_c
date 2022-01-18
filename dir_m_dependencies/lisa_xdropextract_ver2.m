%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run lisa_setprefix_ver2 first. ;
% run lisa_setnames_ver2 first. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% extract rdrop_a and cdrop_a. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fname_tmp = sprintf('%s/out_xdrop_a.txt',dir_out_s0000); fcheck(fname_tmp);
out_xdrop_a_ = textread(fname_tmp); rdrop_a_ = out_xdrop_a_(find(out_xdrop_a_(:,1)>-1),1)+1; cdrop_a_ = out_xdrop_a_(find(out_xdrop_a_(:,2)>-1),2)+1;
fname_tmp = sprintf('%s/out_xdrop_b.txt',dir_out_s0000); fcheck(fname_tmp);
out_xdrop_b_ = textread(fname_tmp); rdrop_b_ = out_xdrop_b_(find(out_xdrop_b_(:,1)>-1),1)+1; cdrop_b_ = out_xdrop_b_(find(out_xdrop_b_(:,2)>-1),2)+1;
