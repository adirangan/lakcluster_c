function lisa_pca_collect_ver8(cl_num_arm1,cl_num_arm2,flag_reverse_use,flag_rerun);
% This function differs from the previous lisa_pca_collect_ver?.m ;
% Here we extract the pca_scores from files such as pca_proj_arm2_ni%d_tst%d_k2_B44_AnV_.mda ;
% Then we threshold by the bicluster direction (i.e., either +1 or -1). ;
% We calculate the cauc for each cutoff along this threshod. ;
%{

  flag_rerun=0;
  for cl_num_arm1=4;
  for cl_num_arm2=[1,2,3];
  for flag_reverse_use=[0,1];
  lisa_pca_collect_ver8(cl_num_arm1,cl_num_arm2,flag_reverse_use,flag_rerun); 
  end;%for flag_reverse_use=[0,1];
  end;%for cl_num_arm2=[1,2,3];
  end;%for cl_num_arm1=4;

  %}

if nargin<1; cl_num_arm1 = 4; end;
if nargin<2; cl_num_arm2 = 1; end;
if nargin<3; flag_reverse_use = 0; end;
% clear; flag_rerun=0; cl_num_arm1 = 4; cl_num_arm2 = 1; flag_reverse_use = 0; 
pca_rank = 2; pca_b_mlt = 44;
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));
%addpath('/data/rangan/dir_bcc/dir_code_022316');
addpath('/data/rangan/dir_bcc/dir_PGC_20190328/dir_m');
dir_trunk = '/data/rangan/dir_bcc/dir_PGC_20190328';
%%%%%%%%;
verbose = 1; n_figure = 1;
B_MLT = 34; n_mds = 20; n_scramble = 1; n_shuffle = 0;
flag_dex_vs_lak = {'dex'};
if (strcmp(flag_dex_vs_lak{1},'dex')); gamma = 0.004; mc_string = ''; end; 
if (strcmp(flag_dex_vs_lak{1},'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
flag_reverse = flag_reverse_use;
cl_num = cl_num_arm1;
n_maf = 5; n_cov = 2;
%%%%%%%%;
lisa_setprefix_ver3 ;
lisa_setnames_ver2 ;
lisa_xdropextract_ver2 ;
lisa_loadmdsfambim_ver2 ;
lisa_mxcheck_ver2 ;
lisa_studyindex_ver2 ;
lisa_traceextract_ver2 ;
%%%%%%%%;
bim_arm1_ = bim_;
fam_arm1_ = fam_;
fam_name_arm1_ = fam_name_;
mds_sort_arm1_ = mds_sort_;
string_prefix_arm1 = string_prefix;
%string_name_arm1 = string_name;
string_name_r1_arm1 = string_name;
string_name_r0_arm1 = string_name_r0;
dir__in_arm1 = dir__in;
%dir_out_s0000_arm1 = dir_out_s0000;
dir_out_s0000_r1_arm1 = dir_out_s0000;
dir_out_s0000_r0_arm1 = dir_out_s0000_r0;
%dir_out_trace_arm1 = dir_out_trace;
dir_out_trace_r1_arm1 = dir_out_trace;
dir_out_trace_r0_arm1 = dir_out_trace_r0;
trace_r1_arm1_ = trace_;
%%%%%%%%;
clear bim_ fam_ mds_sort_ mds_;
%%%%%%%%;
verbose = 1;
B_MLT = 34; n_mds = 20; n_scramble = 1; n_shuffle = 0;
flag_dex_vs_lak = {'dex'};
if (strcmp(flag_dex_vs_lak{1},'dex')); gamma = 0.004; mc_string = ''; end; 
if (strcmp(flag_dex_vs_lak{1},'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
flag_reverse = flag_reverse_use;
cl_num = cl_num_arm2;
n_maf = 5; n_cov = 2;
%%%%%%%%;
lisa_setprefix_ver3 ;
lisa_setnames_ver2 ;
lisa_xdropextract_ver2 ;
lisa_loadmdsfambim_ver2 ;
lisa_mxcheck_ver2 ;
lisa_studyindex_ver2 ;
lisa_traceextract_ver2 ;
%%%%%%%%;
study_trunk_arm2_ = study_trunk_;
study_name_arm2_ = study_name_;
n_study_arm2 = n_study;
bim_arm2_ = bim_;
fam_arm2_ = fam_;
fam_name_arm2_ = fam_name_;
mds_sort_arm2_ = mds_sort_;
string_prefix_arm2 = string_prefix;
%string_name_arm2 = string_name;
string_name_r1_arm2 = string_name;
string_name_r0_arm2 = string_name_r0;
dir__in_arm2 = dir__in;
%dir_out_s0000_arm2 = dir_out_s0000;
dir_out_s0000_r1_arm2 = dir_out_s0000;
dir_out_s0000_r0_arm2 = dir_out_s0000_r0;
%dir_out_trace_arm2 = dir_out_trace;
dir_out_trace_r1_arm2 = dir_out_trace;
dir_out_trace_r0_arm2 = dir_out_trace_r0;
trace_r1_arm2_ = trace_;
%%%%%%%%;
clear bim_ fam_ mds_sort_ mds_;
%%%%%%%%;

%%%%%%%%%%%%%%%%;
% copy fam file. ;
%%%%%%%%%%%%%%%%;
ni=1;
fam_fam_arm1_ = fam_arm1_{ni}; ni=ni+1;
fam_iid_arm1_ = fam_arm1_{ni}; ni=ni+1;
fam_pat_arm1_ = fam_arm1_{ni}; ni=ni+1;
fam_mat_arm1_ = fam_arm1_{ni}; ni=ni+1;
fam_sex_arm1_ = fam_arm1_{ni}; ni=ni+1;
fam_phe_arm1_ = fam_arm1_{ni}; ni=ni+1;
fam_dir_arm1_ = fam_arm1_{ni}; ni=ni+1;
%%%%%%%%%%%%%%%%;
% extract mc;
%%%%%%%%%%%%%%%%;
if (flag_reverse_use==0);
tmp_fname = sprintf('%s/%s_mr_A_full.b16',dir__in_arm1,string_prefix_arm1); fcheck(tmp_fname);
mr_D_arm1_ = tutorial_binary_uncompress(tmp_fname)>0;
tmp_fname = sprintf('%s/%s_mr_Z_full.b16',dir__in_arm1,string_prefix_arm1); fcheck(tmp_fname);
mr_X_arm1_ = tutorial_binary_uncompress(tmp_fname)>0;
end;%if (flag_reverse_use==0);
if (flag_reverse_use==1);
tmp_fname = sprintf('%s/%s_mr_Z_full.b16',dir__in_arm1,string_prefix_arm1); fcheck(tmp_fname);
mr_D_arm1_ = tutorial_binary_uncompress(tmp_fname)>0;
tmp_fname = sprintf('%s/%s_mr_A_full.b16',dir__in_arm1,string_prefix_arm1); fcheck(tmp_fname);
mr_X_arm1_ = tutorial_binary_uncompress(tmp_fname)>0;
end;%if (flag_reverse_use==1);
%%%%%%%%%%%%%%%%;
if (flag_reverse_use==0);
disp(sprintf('mr_D_arm1_ error: %0.16f',norm(cast(mr_D_arm1_,'double')-cast(fam_phe_arm1_==2,'double'))));
disp(sprintf('mr_X_arm1_ error: %0.16f',norm(cast(mr_X_arm1_,'double')-cast(fam_phe_arm1_==1,'double'))));
end;%if (flag_reverse_use==0);
if (flag_reverse_use==1);
disp(sprintf('mr_D_arm1_ error: %0.16f',norm(cast(mr_D_arm1_,'double')-cast(fam_phe_arm1_==1,'double'))));
disp(sprintf('mr_X_arm1_ error: %0.16f',norm(cast(mr_X_arm1_,'double')-cast(fam_phe_arm1_==2,'double'))));
end;%if (flag_reverse_use==1);

%%%%%%%%%%%%%%%%;
% copy fam file. ;
%%%%%%%%%%%%%%%%;
ni=1;
fam_fam_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_iid_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_pat_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_mat_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_sex_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_phe_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_dir_arm2_ = fam_arm2_{ni}; ni=ni+1;
%%%%%%%%%%%%%%%%;
% extract mc;
%%%%%%%%%%%%%%%%;
if (flag_reverse_use==0);
tmp_fname = sprintf('%s/%s_mr_A_full.b16',dir__in_arm2,string_prefix_arm2); fcheck(tmp_fname);
mr_D_arm2_ = tutorial_binary_uncompress(tmp_fname)>0;
tmp_fname = sprintf('%s/%s_mr_Z_full.b16',dir__in_arm2,string_prefix_arm2); fcheck(tmp_fname);
mr_X_arm2_ = tutorial_binary_uncompress(tmp_fname)>0;
end;%if (flag_reverse_use==0);
if (flag_reverse_use==1);
tmp_fname = sprintf('%s/%s_mr_Z_full.b16',dir__in_arm2,string_prefix_arm2); fcheck(tmp_fname);
mr_D_arm2_ = tutorial_binary_uncompress(tmp_fname)>0;
tmp_fname = sprintf('%s/%s_mr_A_full.b16',dir__in_arm2,string_prefix_arm2); fcheck(tmp_fname);
mr_X_arm2_ = tutorial_binary_uncompress(tmp_fname)>0;
end;%if (flag_reverse_use==1);
%%%%%%%%%%%%%%%%;
if (flag_reverse_use==0);
disp(sprintf('mr_D_arm2_ error: %0.16f',norm(cast(mr_D_arm2_,'double')-cast(fam_phe_arm2_==2,'double'))));
disp(sprintf('mr_X_arm2_ error: %0.16f',norm(cast(mr_X_arm2_,'double')-cast(fam_phe_arm2_==1,'double'))));
end;%if (flag_reverse_use==0);
if (flag_reverse_use==1);
disp(sprintf('mr_D_arm2_ error: %0.16f',norm(cast(mr_D_arm2_,'double')-cast(fam_phe_arm2_==1,'double'))));
disp(sprintf('mr_X_arm2_ error: %0.16f',norm(cast(mr_X_arm2_,'double')-cast(fam_phe_arm2_==2,'double'))));
end;%if (flag_reverse_use==1);

dir_out_s0000_pca_r0_arm1 = sprintf('%s/dir_pca',dir_out_s0000_r0_arm1);
if (~exist(dir_out_s0000_pca_r0_arm1,'dir')); disp(sprintf(' %% creating %s',dir_out_s0000_pca_r0_arm1)); mkdir(dir_out_s0000_pca_r0_arm1); end;
dir_out_s0000_pca_r1_arm1 = sprintf('%s/dir_pca',dir_out_s0000_r1_arm1);
if (~exist(dir_out_s0000_pca_r1_arm1,'dir')); disp(sprintf(' %% creating %s',dir_out_s0000_pca_r1_arm1)); mkdir(dir_out_s0000_pca_r1_arm1); end;

n_iteration_arm1 = length(trace_r1_arm1_{1});
n_iteration_arm2 = length(trace_r1_arm2_{1});
niteration_ = 0:25:n_iteration_arm1; 
niteration_(1)=1; niteration_(end) = min(niteration_(end),n_iteration_arm1-1);
niteration_ = niteration_(length(niteration_):-1:1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
pca_fname_mat = sprintf('%s/%s_pca_collect_ver8_trn%d_tst%d.mat',dir_out_s0000_pca_r0_arm1,string_name_r0_arm1,cl_num_arm1,cl_num_arm2); 
if (~flag_rerun &  exist(pca_fname_mat,'file')); disp(sprintf(' %% found %s, loading',pca_fname_mat)); load(pca_fname_mat); end;
if ( flag_rerun | ~exist(pca_fname_mat,'file')); 
disp(sprintf(' %% could not find %s, creating',pca_fname_mat));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%;
% build libraries for cauc p-value ;
%%%%%%%%%%%%%%%%;
find_X_arm1_ = find(mr_X_arm1_); find_D_arm1_ = find(mr_D_arm1_);
[cauc02_arm1_avg_,cauc02_arm1_min_] = p_cauc_0(mds_sort_arm1_(find_X_arm1_,1: 2),mds_sort_arm1_(find_D_arm1_,1: 2),zeros(length(find_X_arm1_),1),zeros(length(find_D_arm1_),1),1024*8);
[cauc06_arm1_avg_,cauc06_arm1_min_] = p_cauc_0(mds_sort_arm1_(find_X_arm1_,1: 6),mds_sort_arm1_(find_D_arm1_,1: 6),zeros(length(find_X_arm1_),1),zeros(length(find_D_arm1_),1),1024*8);
[cauc19_arm1_avg_,cauc19_arm1_min_] = p_cauc_0(mds_sort_arm1_(find_X_arm1_,[1:6,19]),mds_sort_arm1_(find_D_arm1_,[1:6,19]),zeros(length(find_X_arm1_),1),zeros(length(find_D_arm1_),1),1024*8);
[cauc20_arm1_avg_,cauc20_arm1_min_] = p_cauc_0(mds_sort_arm1_(find_X_arm1_,1:20),mds_sort_arm1_(find_D_arm1_,1:20),zeros(length(find_X_arm1_),1),zeros(length(find_D_arm1_),1),1024*8);
find_X_arm2_ = find(mr_X_arm2_); find_D_arm2_ = find(mr_D_arm2_);
[cauc02_arm2_avg_,cauc02_arm2_min_] = p_cauc_0(mds_sort_arm2_(find_X_arm2_,1: 2),mds_sort_arm2_(find_D_arm2_,1: 2),zeros(length(find_X_arm2_),1),zeros(length(find_D_arm2_),1),1024*8);
[cauc06_arm2_avg_,cauc06_arm2_min_] = p_cauc_0(mds_sort_arm2_(find_X_arm2_,1: 6),mds_sort_arm2_(find_D_arm2_,1: 6),zeros(length(find_X_arm2_),1),zeros(length(find_D_arm2_),1),1024*8);
[cauc19_arm2_avg_,cauc19_arm2_min_] = p_cauc_0(mds_sort_arm2_(find_X_arm2_,[1:6,19]),mds_sort_arm2_(find_D_arm2_,[1:6,19]),zeros(length(find_X_arm2_),1),zeros(length(find_D_arm2_),1),1024*8);
[cauc20_arm2_avg_,cauc20_arm2_min_] = p_cauc_0(mds_sort_arm2_(find_X_arm2_,1:20),mds_sort_arm2_(find_D_arm2_,1:20),zeros(length(find_X_arm2_),1),zeros(length(find_D_arm2_),1),1024*8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%;
%extract rdrop ;
%%%%%%%%%%%%%%%%;
out_xdrop_a_r0_arm1_ = textread(sprintf('%s/out_xdrop_a.txt',dir_out_s0000_r0_arm1));
tmp_ = out_xdrop_a_r0_arm1_(:,1); tmp_ = tmp_(find(tmp_>-1)); %tmp_ = tmp_(end:-1:1); 
tmp_ = tmp_+1; rdrop_a_r0_arm1_ = tmp_; clear out_xdrop_a_r0_arm1_; clear tmp_;
out_xdrop_a_r1_arm1_ = textread(sprintf('%s/out_xdrop_a.txt',dir_out_s0000_r1_arm1));
tmp_ = out_xdrop_a_r1_arm1_(:,1); tmp_ = tmp_(find(tmp_>-1)); %tmp_ = tmp_(end:-1:1); 
tmp_ = tmp_+1; rdrop_a_r1_arm1_ = tmp_; clear out_xdrop_a_r1_arm1_; clear tmp_;
out_xdrop_a_r0_arm2_ = textread(sprintf('%s/out_xdrop_a.txt',dir_out_s0000_r0_arm2));
tmp_ = out_xdrop_a_r0_arm2_(:,1); tmp_ = tmp_(find(tmp_>-1)); %tmp_ = tmp_(end:-1:1); 
tmp_ = tmp_+1; rdrop_a_r0_arm2_ = tmp_; clear out_xdrop_a_r0_arm2_; clear tmp_;
out_xdrop_a_r1_arm2_ = textread(sprintf('%s/out_xdrop_a.txt',dir_out_s0000_r1_arm2));
tmp_ = out_xdrop_a_r1_arm2_(:,1); tmp_ = tmp_(find(tmp_>-1)); %tmp_ = tmp_(end:-1:1); 
tmp_ = tmp_+1; rdrop_a_r1_arm2_ = tmp_; clear out_xdrop_a_r1_arm2_; clear tmp_;
%%%%%%%%%%%%%%%%;
%extract trace ;
%%%%%%%%%%%%%%%%;
out_trace_r0_arm1_ = textread(sprintf('%s/out_trace.txt',dir_out_s0000_r0_arm1));
out_trace_r1_arm1_ = textread(sprintf('%s/out_trace.txt',dir_out_s0000_r1_arm1));
out_trace_r0_arm2_ = textread(sprintf('%s/out_trace.txt',dir_out_s0000_r0_arm2));
out_trace_r1_arm2_ = textread(sprintf('%s/out_trace.txt',dir_out_s0000_r1_arm2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_threshold = 21; 
auc19_r0_arm1_ = zeros(length(niteration_),1);
auc19_r0_arm2_ = zeros(length(niteration_),1);
auc19_r1_arm1_ = zeros(length(niteration_),1);
auc19_r1_arm2_ = zeros(length(niteration_),1);
logp_auc19_r0_arm1_ = zeros(length(niteration_),1);
logp_auc19_r0_arm2_ = zeros(length(niteration_),1);
logp_auc19_r1_arm1_ = zeros(length(niteration_),1);
logp_auc19_r1_arm2_ = zeros(length(niteration_),1);
%%%%%%%%;
auc19_r0_arm1__ = zeros(length(niteration_),n_threshold);
auc19_r0_arm2__ = zeros(length(niteration_),n_threshold);
auc19_r1_arm1__ = zeros(length(niteration_),n_threshold);
auc19_r1_arm2__ = zeros(length(niteration_),n_threshold);
for ni=1:length(niteration_);
niteration = niteration_(ni);
disp(sprintf(' %% ni %d niteration %d',ni,niteration));
%%%%%%%%;
pca_proj_infix_r0_arm1 = sprintf('pca_proj_arm1_ni%d_tst%d',niteration,cl_num_arm2); %<-- Note that here we use the projection of arm1, with snps limited by overlap with arm two ;
DnV_r0_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_AnV_.mda',dir_out_s0000_pca_r0_arm1,pca_proj_infix_r0_arm1,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca_r0 from arm one, since this is where the bicluster is. ;
XnV_r0_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_ZnV_.mda',dir_out_s0000_pca_r0_arm1,pca_proj_infix_r0_arm1,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca_r0 from arm one, since this is where the bicluster is. ;
%if (strcmp(flag_dex_vs_lak,'lak')); DnV_r0_arm1_ = abs(DnV_r0_arm1_); XnV_r0_arm1_ = abs(XnV_r0_arm1_); end;%if (strcmp(flag_dex_vs_lak,'lak'));
mr_D_rmv_r0_arm1_ = mr_D_arm1_*0; mr_D_rmv_r0_arm1_(rdrop_a_r0_arm1_(1:end-out_trace_r0_arm1_(niteration,2)))=1;
mr_D_ret_r0_arm1_ = mr_D_arm1_*0; mr_D_ret_r0_arm1_(rdrop_a_r0_arm1_(end-out_trace_r0_arm1_(niteration,2)+1:end))=1;
%%%%%%%%;
pca_proj_infix_r0_arm2 = sprintf('pca_proj_arm2_ni%d_tst%d',niteration,cl_num_arm2); %<-- Note that here we use the projection of arm2, with snps limited by overlap with arm two ;
DnV_r0_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_AnV_.mda',dir_out_s0000_pca_r0_arm1,pca_proj_infix_r0_arm2,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca_r0 from arm one, since this is where the bicluster is. ;
XnV_r0_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_ZnV_.mda',dir_out_s0000_pca_r0_arm1,pca_proj_infix_r0_arm2,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca_r0 from arm one, since this is where the bicluster is. ;
%if (strcmp(flag_dex_vs_lak,'lak')); DnV_r0_arm2_ = abs(DnV_r0_arm2_); XnV_r0_arm2_ = abs(XnV_r0_arm2_); end;%if (strcmp(flag_dex_vs_lak,'lak'));
%%%%%%%%;
pca_proj_infix_r1_arm1 = sprintf('pca_proj_arm1_ni%d_tst%d',niteration,cl_num_arm2); %<-- Note that here we use the projection of arm1, with snps limited by overlap with arm two ;
DnV_r1_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_AnV_.mda',dir_out_s0000_pca_r1_arm1,pca_proj_infix_r1_arm1,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca_r1 from arm one, since this is where the bicluster is. ;
XnV_r1_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_ZnV_.mda',dir_out_s0000_pca_r1_arm1,pca_proj_infix_r1_arm1,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca_r1 from arm one, since this is where the bicluster is. ;
%if (strcmp(flag_dex_vs_lak,'lak')); DnV_r1_arm1_ = abs(DnV_r1_arm1_); XnV_r1_arm1_ = abs(XnV_r1_arm1_); end;%if (strcmp(flag_dex_vs_lak,'lak'));
mr_D_rmv_r1_arm1_ = mr_D_arm1_*0; mr_D_rmv_r1_arm1_(rdrop_a_r1_arm1_(1:end-out_trace_r1_arm1_(niteration,2)))=1;
mr_D_ret_r1_arm1_ = mr_D_arm1_*0; mr_D_ret_r1_arm1_(rdrop_a_r1_arm1_(end-out_trace_r1_arm1_(niteration,2)+1:end))=1;
%%%%%%%%;
pca_proj_infix_r1_arm2 = sprintf('pca_proj_arm2_ni%d_tst%d',niteration,cl_num_arm2); %<-- Note that here we use the projection of arm2, with snps limited by overlap with arm two ;
DnV_r1_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_AnV_.mda',dir_out_s0000_pca_r1_arm1,pca_proj_infix_r1_arm2,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca_r1 from arm one, since this is where the bicluster is. ;
XnV_r1_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_ZnV_.mda',dir_out_s0000_pca_r1_arm1,pca_proj_infix_r1_arm2,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca_r1 from arm one, since this is where the bicluster is. ;
%if (strcmp(flag_dex_vs_lak,'lak')); DnV_r1_arm2_ = abs(DnV_r1_arm2_); XnV_r1_arm2_ = abs(XnV_r1_arm2_); end;%if (strcmp(flag_dex_vs_lak,'lak'));
flag_disp=0;
if flag_disp;
%%%%%%%%;
figure(n_figure); clf; n_figure = n_figure+1;
%%%%%%%%;
subplot(2,2,1);hold on; 
plot(XnV_r0_arm1_(find_X_arm1_,1),XnV_r0_arm1_(find_X_arm1_,2),'k.','MarkerSize',15);
find_tmp_ = find(mr_D_rmv_r0_arm1_); plot(DnV_r0_arm1_(find_tmp_,1),DnV_r0_arm1_(find_tmp_,2),'g.','MarkerSize',15);
find_tmp_ = find(mr_D_ret_r0_arm1_); plot(DnV_r0_arm1_(find_tmp_,1),DnV_r0_arm1_(find_tmp_,2),'r.','MarkerSize',25);
hold off;
xlabel('PC1');ylabel('PC2'); title(sprintf('ni %d r0 arm1',niteration));
%%%%%%%%;
subplot(2,2,3); hold on;
plot(XnV_r0_arm2_(find_X_arm2_,1),XnV_r0_arm2_(find_X_arm2_,2),'k.','MarkerSize',15);
plot(DnV_r0_arm2_(find_D_arm2_,1),DnV_r0_arm2_(find_D_arm2_,2),'m.','MarkerSize',15);
hold off;
xlabel('PC1');ylabel('PC2'); title(sprintf('ni %d r0 arm2',niteration));
%%%%%%%%;
subplot(2,2,2);hold on; 
plot(XnV_r1_arm1_(find_X_arm1_,1),XnV_r1_arm1_(find_X_arm1_,2),'k.','MarkerSize',15);
find_tmp_ = find(mr_D_rmv_r1_arm1_); plot(DnV_r1_arm1_(find_tmp_,1),DnV_r1_arm1_(find_tmp_,2),'g.','MarkerSize',15);
find_tmp_ = find(mr_D_ret_r1_arm1_); plot(DnV_r1_arm1_(find_tmp_,1),DnV_r1_arm1_(find_tmp_,2),'r.','MarkerSize',25);
hold off;
xlabel('PC1');ylabel('PC2'); title(sprintf('ni %d r1 arm1',niteration));
%%%%%%%%;
subplot(2,2,4); hold on;
plot(XnV_r1_arm2_(find_X_arm2_,1),XnV_r1_arm2_(find_X_arm2_,2),'k.','MarkerSize',15);
plot(DnV_r1_arm2_(find_D_arm2_,1),DnV_r1_arm2_(find_D_arm2_,2),'m.','MarkerSize',15);
hold off;
xlabel('PC1');ylabel('PC2'); title(sprintf('ni %d r1 arm2',niteration));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
set(gcf,'Position',1+[0,0,1024*2,1024]);
%%%%%%%%;
end;%if flag_disp;
%%%
%%
%
%%
%%%
%%%%%%%%;
tmp_auc19_r0_arm1 = cauc_0(XnV_r0_arm1_(find_X_arm1_,1),DnV_r0_arm1_(find_D_arm1_,1),mds_sort_arm1_(find_X_arm1_,[1:6,19]),mds_sort_arm1_(find_D_arm1_,[1:6,19])); flag_flip_r0_arm1 = (tmp_auc19_r0_arm1<0.5);
disp(sprintf(' %% tmp_auc19_r0_arm1 %0.2f flag_flip_r0_arm1 %d',tmp_auc19_r0_arm1,flag_flip_r0_arm1));
tmp_auc19_r0_arm2 = cauc_0(XnV_r0_arm2_(find_X_arm2_,1),DnV_r0_arm2_(find_D_arm2_,1),mds_sort_arm2_(find_X_arm2_,[1:6,19]),mds_sort_arm2_(find_D_arm2_,[1:6,19])); flag_flip_r0_arm2 = (tmp_auc19_r0_arm2<0.5);
disp(sprintf(' %% tmp_auc19_r0_arm2 %0.2f flag_flip_r0_arm2 %d',tmp_auc19_r0_arm2,flag_flip_r0_arm2));
tmp_auc19_r1_arm1 = cauc_0(XnV_r1_arm1_(find_X_arm1_,1),DnV_r1_arm1_(find_D_arm1_,1),mds_sort_arm1_(find_X_arm1_,[1:6,19]),mds_sort_arm1_(find_D_arm1_,[1:6,19])); flag_flip_r1_arm1 = (tmp_auc19_r1_arm1<0.5);
disp(sprintf(' %% tmp_auc19_r1_arm1 %0.2f flag_flip_r1_arm1 %d',tmp_auc19_r1_arm1,flag_flip_r1_arm1));
tmp_auc19_r1_arm2 = cauc_0(XnV_r1_arm2_(find_X_arm2_,1),DnV_r1_arm2_(find_D_arm2_,1),mds_sort_arm2_(find_X_arm2_,[1:6,19]),mds_sort_arm2_(find_D_arm2_,[1:6,19])); flag_flip_r1_arm2 = (tmp_auc19_r1_arm2<0.5);
disp(sprintf(' %% tmp_auc19_r1_arm2 %0.2f flag_flip_r1_arm2 %d',tmp_auc19_r1_arm2,flag_flip_r1_arm2));
%%%%%%%%;
if ~flag_flip_r0_arm1; tmp_DnV_r0_arm1_ = +DnV_r0_arm1_; tmp_XnV_r0_arm1_ = +XnV_r0_arm1_; tmp_DnV_r0_arm2_ = +DnV_r0_arm2_; tmp_XnV_r0_arm2_ = +XnV_r0_arm2_; end;%if ~flag_flip_r0_arm1;
if  flag_flip_r0_arm1; tmp_DnV_r0_arm1_ = -DnV_r0_arm1_; tmp_XnV_r0_arm1_ = -XnV_r0_arm1_; tmp_DnV_r0_arm2_ = -DnV_r0_arm2_; tmp_XnV_r0_arm2_ = -XnV_r0_arm2_; end;%if  flag_flip_r0_arm1;
if ~flag_flip_r1_arm1; tmp_DnV_r1_arm1_ = +DnV_r1_arm1_; tmp_XnV_r1_arm1_ = +XnV_r1_arm1_; tmp_DnV_r1_arm2_ = +DnV_r1_arm2_; tmp_XnV_r1_arm2_ = +XnV_r1_arm2_; end;%if ~flag_flip_r1_arm1;
if  flag_flip_r1_arm1; tmp_DnV_r1_arm1_ = -DnV_r1_arm1_; tmp_XnV_r1_arm1_ = -XnV_r1_arm1_; tmp_DnV_r1_arm2_ = -DnV_r1_arm2_; tmp_XnV_r1_arm2_ = -XnV_r1_arm2_; end;%if  flag_flip_r1_arm1;
%%%%%%%%;
tmp_mds_X_ = mds_sort_arm1_(find_X_arm1_,[1:6,19]); tmp_mds_D_ = mds_sort_arm1_(find_D_arm1_,[1:6,19]);
tmp_auc19_r0_arm1 = cauc_0(tmp_XnV_r0_arm1_(find_X_arm1_,1),tmp_DnV_r0_arm1_(find_D_arm1_,1),tmp_mds_X_,tmp_mds_D_);
tmp_logp_auc19_r0_arm1 = log(max(1,length(find(cauc19_arm1_avg_>tmp_auc19_r0_arm1)))/length(cauc19_arm1_avg_));
auc19_r0_arm1_(ni) = tmp_auc19_r0_arm1;
logp_auc19_r0_arm1_(ni) = tmp_logp_auc19_r0_arm1;
%%%%%%%%;
tmp_mds_X_ = mds_sort_arm2_(find_X_arm2_,[1:6,19]); tmp_mds_D_ = mds_sort_arm2_(find_D_arm2_,[1:6,19]);
tmp_auc19_r0_arm2 = cauc_0(tmp_XnV_r0_arm2_(find_X_arm2_,1),tmp_DnV_r0_arm2_(find_D_arm2_,1),tmp_mds_X_,tmp_mds_D_);
tmp_logp_auc19_r0_arm2 = log(max(1,length(find(cauc19_arm2_avg_>tmp_auc19_r0_arm2)))/length(cauc19_arm2_avg_));
auc19_r0_arm2_(ni) = tmp_auc19_r0_arm2;
logp_auc19_r0_arm2_(ni) = tmp_logp_auc19_r0_arm2;
%%%%%%%%;
tmp_mds_X_ = mds_sort_arm1_(find_X_arm1_,[1:6,19]); tmp_mds_D_ = mds_sort_arm1_(find_D_arm1_,[1:6,19]);
tmp_auc19_r1_arm1 = cauc_0(tmp_XnV_r1_arm1_(find_X_arm1_,1),tmp_DnV_r1_arm1_(find_D_arm1_,1),tmp_mds_X_,tmp_mds_D_);
tmp_logp_auc19_r1_arm1 = log(max(1,length(find(cauc19_arm1_avg_>tmp_auc19_r1_arm1)))/length(cauc19_arm1_avg_));
auc19_r1_arm1_(ni) = tmp_auc19_r1_arm1;
logp_auc19_r1_arm1_(ni) = tmp_logp_auc19_r1_arm1;
%%%%%%%%;
tmp_mds_X_ = mds_sort_arm2_(find_X_arm2_,[1:6,19]); tmp_mds_D_ = mds_sort_arm2_(find_D_arm2_,[1:6,19]);
tmp_auc19_r1_arm2 = cauc_0(tmp_XnV_r1_arm2_(find_X_arm2_,1),tmp_DnV_r1_arm2_(find_D_arm2_,1),tmp_mds_X_,tmp_mds_D_);
tmp_logp_auc19_r1_arm2 = log(max(1,length(find(cauc19_arm2_avg_>tmp_auc19_r1_arm2)))/length(cauc19_arm2_avg_));
auc19_r1_arm2_(ni) = tmp_auc19_r1_arm2;
logp_auc19_r1_arm2_(ni) = tmp_logp_auc19_r1_arm2;
%%%%%%%%;
threshold_r0_arm1_ = prctile( tmp_DnV_r0_arm1_(find_D_arm1_,1), linspace(0,100,n_threshold) );
threshold_r0_arm2_ = prctile( tmp_DnV_r0_arm2_(find_D_arm2_,1), linspace(0,100,n_threshold) );
threshold_r1_arm1_ = prctile( tmp_DnV_r1_arm1_(find_D_arm1_,1), linspace(0,100,n_threshold) );
threshold_r1_arm2_ = prctile( tmp_DnV_r1_arm2_(find_D_arm2_,1), linspace(0,100,n_threshold) );
%%%%%%%%;
for nthreshold=1:n_threshold;
threshold_r0_arm1 = threshold_r0_arm1_(nthreshold);
tmp_find_X_r0_arm1_ = find(tmp_XnV_r0_arm1_(find_X_arm1_,1)>=threshold_r0_arm1); tmp_find_D_r0_arm1_ = find(tmp_DnV_r0_arm1_(find_D_arm1_,1)>=threshold_r0_arm1);
tmp_auc19 = 0.5;
if (length(tmp_find_X_r0_arm1_)==0 & length(tmp_find_D_r0_arm1_)==0); tmp_auc19 = 0.5; end;
if (length(tmp_find_X_r0_arm1_)==0 & length(tmp_find_D_r0_arm1_)> 0); tmp_auc19 = 1.0; end;
if (length(tmp_find_X_r0_arm1_)> 0 & length(tmp_find_D_r0_arm1_)==0); tmp_auc19 = 0.0; end;
if (length(tmp_find_X_r0_arm1_)> 0 & length(tmp_find_D_r0_arm1_)> 0); 
tmp_X_ = tmp_XnV_r0_arm1_(find_X_arm1_(tmp_find_X_r0_arm1_),1); tmp_D_ = tmp_DnV_r0_arm1_(find_D_arm1_(tmp_find_D_r0_arm1_),1);
tmp_mds_X_ = mds_sort_arm1_(find_X_arm1_(tmp_find_X_r0_arm1_),[1:6,19]); tmp_mds_D_ = mds_sort_arm1_(find_D_arm1_(tmp_find_D_r0_arm1_),[1:6,19]);
tmp_auc19 = cauc_0(tmp_X_,tmp_D_,tmp_mds_X_,tmp_mds_D_);
end;%if (length(tmp_find_X_r0_arm1_)> 0 & length(tmp_find_D_r0_arm1_)> 0); 
auc19_r0_arm1__(ni,nthreshold) = tmp_auc19;
end;%for nthreshold=1:n_threshold;
%%%%%%%%;
for nthreshold=1:n_threshold;
threshold_r0_arm2 = threshold_r0_arm2_(nthreshold);
tmp_find_X_r0_arm2_ = find(tmp_XnV_r0_arm2_(find_X_arm2_,1)>=threshold_r0_arm2); tmp_find_D_r0_arm2_ = find(tmp_DnV_r0_arm2_(find_D_arm2_,1)>=threshold_r0_arm2);
tmp_auc19 = 0.5;
if (length(tmp_find_X_r0_arm2_)==0 & length(tmp_find_D_r0_arm2_)==0); tmp_auc19 = 0.5; end;
if (length(tmp_find_X_r0_arm2_)==0 & length(tmp_find_D_r0_arm2_)> 0); tmp_auc19 = 1.0; end;
if (length(tmp_find_X_r0_arm2_)> 0 & length(tmp_find_D_r0_arm2_)==0); tmp_auc19 = 0.0; end;
if (length(tmp_find_X_r0_arm2_)> 0 & length(tmp_find_D_r0_arm2_)> 0); 
tmp_X_ = tmp_XnV_r0_arm2_(find_X_arm2_(tmp_find_X_r0_arm2_),1); tmp_D_ = tmp_DnV_r0_arm2_(find_D_arm2_(tmp_find_D_r0_arm2_),1);
tmp_mds_X_ = mds_sort_arm2_(find_X_arm2_(tmp_find_X_r0_arm2_),[1:6,19]); tmp_mds_D_ = mds_sort_arm2_(find_D_arm2_(tmp_find_D_r0_arm2_),[1:6,19]);
tmp_auc19 = cauc_0(tmp_X_,tmp_D_,tmp_mds_X_,tmp_mds_D_);
end;%if (length(tmp_find_X_r0_arm2_)> 0 & length(tmp_find_D_r0_arm2_)> 0); 
auc19_r0_arm2__(ni,nthreshold) = tmp_auc19;
end;%for nthreshold=1:n_threshold;
%%%%%%%%;
for nthreshold=1:n_threshold;
threshold_r1_arm1 = threshold_r1_arm1_(nthreshold);
tmp_find_X_r1_arm1_ = find(tmp_XnV_r1_arm1_(find_X_arm1_,1)>=threshold_r1_arm1); tmp_find_D_r1_arm1_ = find(tmp_DnV_r1_arm1_(find_D_arm1_,1)>=threshold_r1_arm1);
tmp_auc19 = 0.5;
if (length(tmp_find_X_r1_arm1_)==0 & length(tmp_find_D_r1_arm1_)==0); tmp_auc19 = 0.5; end;
if (length(tmp_find_X_r1_arm1_)==0 & length(tmp_find_D_r1_arm1_)> 0); tmp_auc19 = 1.0; end;
if (length(tmp_find_X_r1_arm1_)> 0 & length(tmp_find_D_r1_arm1_)==0); tmp_auc19 = 0.0; end;
if (length(tmp_find_X_r1_arm1_)> 0 & length(tmp_find_D_r1_arm1_)> 0); 
tmp_X_ = tmp_XnV_r1_arm1_(find_X_arm1_(tmp_find_X_r1_arm1_),1); tmp_D_ = tmp_DnV_r1_arm1_(find_D_arm1_(tmp_find_D_r1_arm1_),1);
tmp_mds_X_ = mds_sort_arm1_(find_X_arm1_(tmp_find_X_r1_arm1_),[1:6,19]); tmp_mds_D_ = mds_sort_arm1_(find_D_arm1_(tmp_find_D_r1_arm1_),[1:6,19]);
tmp_auc19 = cauc_0(tmp_X_,tmp_D_,tmp_mds_X_,tmp_mds_D_);
end;%if (length(tmp_find_X_r1_arm1_)> 0 & length(tmp_find_D_r1_arm1_)> 0); 
auc19_r1_arm1__(ni,nthreshold) = tmp_auc19;
end;%for nthreshold=1:n_threshold;
%%%%%%%%;
for nthreshold=1:n_threshold;
threshold_r1_arm2 = threshold_r1_arm2_(nthreshold);
tmp_find_X_r1_arm2_ = find(tmp_XnV_r1_arm2_(find_X_arm2_,1)>=threshold_r1_arm2); tmp_find_D_r1_arm2_ = find(tmp_DnV_r1_arm2_(find_D_arm2_,1)>=threshold_r1_arm2);
tmp_auc19 = 0.5;
if (length(tmp_find_X_r1_arm2_)==0 & length(tmp_find_D_r1_arm2_)==0); tmp_auc19 = 0.5; end;
if (length(tmp_find_X_r1_arm2_)==0 & length(tmp_find_D_r1_arm2_)> 0); tmp_auc19 = 1.0; end;
if (length(tmp_find_X_r1_arm2_)> 0 & length(tmp_find_D_r1_arm2_)==0); tmp_auc19 = 0.0; end;
if (length(tmp_find_X_r1_arm2_)> 0 & length(tmp_find_D_r1_arm2_)> 0); 
tmp_X_ = tmp_XnV_r1_arm2_(find_X_arm2_(tmp_find_X_r1_arm2_),1); tmp_D_ = tmp_DnV_r1_arm2_(find_D_arm2_(tmp_find_D_r1_arm2_),1);
tmp_mds_X_ = mds_sort_arm2_(find_X_arm2_(tmp_find_X_r1_arm2_),[1:6,19]); tmp_mds_D_ = mds_sort_arm2_(find_D_arm2_(tmp_find_D_r1_arm2_),[1:6,19]);
tmp_auc19 = cauc_0(tmp_X_,tmp_D_,tmp_mds_X_,tmp_mds_D_);
end;%if (length(tmp_find_X_r1_arm2_)> 0 & length(tmp_find_D_r1_arm2_)> 0); 
auc19_r1_arm2__(ni,nthreshold) = tmp_auc19;
end;%for nthreshold=1:n_threshold;
%%%%%%%%;
clear tmp_X_ tmp_D_ tmp_mds_X_ tmp_mds_D_ tmp_auc19
clear tmp_DnV_r0_arm1_; clear tmp_XnV_r0_arm1_; clear tmp_DnV_r0_arm2_; clear tmp_XnV_r0_arm2_;
clear tmp_DnV_r1_arm1_; clear tmp_XnV_r1_arm1_; clear tmp_DnV_r1_arm2_; clear tmp_XnV_r1_arm2_;
end;%for ni=1:length(niteration_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
save(pca_fname_mat,...
'niteration_',...
'n_threshold',...
'threshold_r0_arm1_',...
'threshold_r0_arm2_',...
'threshold_r1_arm1_',...
'threshold_r1_arm2_',...
'auc19_r0_arm1_',...
'auc19_r0_arm2_',...
'auc19_r1_arm1_',...
'auc19_r1_arm2_',...
'logp_auc19_r0_arm1_',...
'logp_auc19_r0_arm2_',...
'logp_auc19_r1_arm1_',...
'logp_auc19_r1_arm2_',...
'auc19_r0_arm1__',...
'auc19_r0_arm2__',...
'auc19_r1_arm1__',...
'auc19_r1_arm2__');
end;%if ( flag_rerun | ~exist(pca_fname_mat,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

disp(sprintf(' %% cl_num_arm1 %d cl_num_arm2 %d',cl_num_arm1,cl_num_arm2));
disp(sprintf(' %% niteration_')); disp(num2str(niteration_));
disp(sprintf(' %% auc19_r0_arm1_')); disp(num2str(transpose(auc19_r0_arm1_)));
disp(sprintf(' %% logp_auc19_r0_arm1_')); disp(num2str(transpose(logp_auc19_r0_arm1_)));
disp(sprintf(' %% auc19_r0_arm2_')); disp(num2str(transpose(auc19_r0_arm2_)));
disp(sprintf(' %% logp_auc19_r0_arm2_')); disp(num2str(transpose(logp_auc19_r0_arm2_)));
disp(sprintf(' %% auc19_r1_arm1_')); disp(num2str(transpose(auc19_r1_arm1_)));
disp(sprintf(' %% logp_auc19_r1_arm1_')); disp(num2str(transpose(logp_auc19_r1_arm1_)));
disp(sprintf(' %% auc19_r1_arm2_')); disp(num2str(transpose(auc19_r1_arm2_)));
disp(sprintf(' %% logp_auc19_r1_arm2_')); disp(num2str(transpose(logp_auc19_r1_arm2_)));
