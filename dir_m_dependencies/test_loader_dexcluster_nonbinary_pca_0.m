function [tmp_E1_pos_,tmp_E0_pos_] = test_loader_dexcluster_nonbinary_pca_0(E_array_,gamma,niteration,prefix,n_rank);

tmp_n_u = size(E_array_,1);
tmp_n_g = size(E_array_,2);
dir_code = '/data/rangan/dir_bcc/dir_lakcluster_c_dev'; dir_trunk = pwd;
rev_flag = 0; A_n_rind_ = {[1:tmp_n_u]}; A_n_cind = 1:tmp_n_g; Z_n_rind_ = {[]}; T_n_ = {ones(tmp_n_u,1)}; T_n_cind = 1;
GLOBAL_TEST_sparse = 0; B_MLT = 34; Ireq = 0;  verbose_flag = 0; flag_force_create = 0; flag_force_fig=0;
%%%%%%%%;
shuffle_num=0;
test_string = sprintf('%s_%s',prefix,lakcluster_uADZSZDA_xfix_gen_ver1(rev_flag,A_n_rind_,Z_n_rind_,T_n_cind,GLOBAL_TEST_sparse,gamma,B_MLT,Ireq,shuffle_num));
%disp(sprintf(' test_string: %s',test_string));
dir__in = sprintf('%s/dir_%s',dir_trunk,prefix);
%disp(sprintf(' dir__in: %s',dir__in));
dir_out = sprintf('%s/dir_%s',dir__in,test_string); 
%disp(sprintf(' dir_out: %s',dir_out));
dir_trace = sprintf('%s/dir_trace',dir_out); 
if (~exist(dir_trace,'dir')); disp(sprintf(' %% mkdir %s',dir_trace)); mkdir(dir_trace); end;
dir_jpg = sprintf('%s/dir_jpg',dir_out); 
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
str_trace = sprintf('%s/out_trace_dexnb_s%.4d.txt',dir_trace,shuffle_num); 
if ( exist(str_trace,'file') & ~flag_force_create);
%disp(sprintf(' %% %s found, not creating',str_trace));
end;%if ( exist(str_trace,'file') & ~flag_force_create);
if (~exist(str_trace,'file') |  flag_force_create);
disp(sprintf(' %% %s not found, must create!',str_trace));
end;%if (~exist(str_trace,'file') |  flag_force_create);
%%%%%%%%;

%%%%%%%%;
str_xdrop = sprintf('%s/out_xdrop_a_dexnb.txt',dir_out);
out_xdrop_E_0_ = textread(str_xdrop);
rdrop_E_0_ = 1+out_xdrop_E_0_(find(out_xdrop_E_0_(:,1)>-1),1); rkeep_E_0_ = rdrop_E_0_(end:-1:1);
cdrop_E_0_ = 1+out_xdrop_E_0_(find(out_xdrop_E_0_(:,2)>-1),2); ckeep_E_0_ = cdrop_E_0_(end:-1:1);
%%%%%%%%;
n_shuffle = 0;
trace_ = cell(1+n_shuffle,1);
for nshuffle=0:n_shuffle;
trace_{1+nshuffle} = textread(sprintf('%s/out_trace_s%0.4d.txt',dir_trace,nshuffle));
end;%for nshuffle=1:n_shuffle;
n_iter_ = trace_{1}(:,1);
r_rem_ = trace_{1}(:,2);
c_rem_ = trace_{1}(:,3);
%%%%%%%%;

%%%%%%%%;
tmp_r_rem_cut = r_rem_(niteration);
tmp_c_rem_cut = c_rem_(niteration);
%%%%%%%%;
[tmp_U_,~,tmp_V_] = svds( E_array_(rkeep_E_0_(1:tmp_r_rem_cut),ckeep_E_0_(1:tmp_c_rem_cut)) ,n_rank);
tmp_E1_pos_ = E_array_(:,ckeep_E_0_(1:tmp_c_rem_cut))*tmp_V_; %<-- original. ;
[tmp_U_,~,tmp_V_] = svds( 1.0*(E_array_(rkeep_E_0_(1:tmp_r_rem_cut),ckeep_E_0_(1:tmp_c_rem_cut))>0) ,n_rank);
tmp_E0_pos_ = (1.0*(E_array_(:,ckeep_E_0_(1:tmp_c_rem_cut))>0))*tmp_V_; %<-- binarized. ;
%%%%%%%%;
