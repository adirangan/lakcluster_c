
dir_weilin = '/home/rangan/dir_bcc/dir_weilin';
dir_mat = sprintf('%s/dir_mat',dir_weilin);

tab_max=41;
for tab=0:tab_max-1;
fname_0in = sprintf('%s/count_asmatrix_1024_full_tab%d.mat',dir_mat,tab);
tmp_ = load(fname_0in);
B_ = sparse(tmp_.A);
fname_out = sprintf('%s/count_asmatrix_1024_sparse_tab%d.mat',dir_mat,tab);
save(fname_out,'B_');
end;%for tab=0:tab_max-1;

B_ = sparse(0,0);
B_logn_ = sparse(0,0);
tab_max=41;
for tab=0:tab_max-1;
fname_0in = sprintf('%s/count_asmatrix_1024_sparse_tab%d.mat',dir_mat,tab);
disp(sprintf(' %% %s',fname_0in));
tmp_ = load(fname_0in);
B_ = [B_ ; tmp_.B_];
B_logn_ = [B_logn_ ; sparse(log(1+tmp_.B_))];
end;%for tab=0:tab_max-1;

fname_mat = sprintf('%s/count_sparse_.mat',dir_mat);
save(fname_mat,'B_');
fname_mat = sprintf('%s/count_logn_sparse_.mat',dir_mat);
save(fname_mat,'B_logn_');

B_sum_1_ = sum(B_,1);
B_sum_2_ = sum(B_,2);
n_smp = size(B_,1);
n_umi = size(B_,2);
