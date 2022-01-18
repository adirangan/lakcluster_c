% Note the following asymmetry regarding sparse matrices: ;
% summing over columns is *much* faster than summing over rows. ;
% (i.e., about 25 times faster for the example below). ;
%{
tmp_N = 1024*8;
tmp_ = sparse(randperm(tmp_N),randperm(tmp_N),1,tmp_N,tmp_N);
disp('sum(*,2)'); tic;for nl=1:1024;sum(tmp_(:,max(1,min(tmp_N,floor(tmp_N*rand())))),2); end;toc;
disp('sum(*,1)'); tic;for nl=1:1024;sum(tmp_(max(1,min(tmp_N,floor(tmp_N*rand()))),:),1); end;toc;
 %}

clear; setup; verbose=1; 
cl_num_arm1 = 4; cl_num_arm2 = 1;
flag_dex_vs_lak = 'dex'; if (strcmp(flag_dex_vs_lak,'lak')); gamma = 0.001; mc_string = ''; else; gamma = 0.004; mc_string = ''; end; 
flag_reverse_use = 0;
B_MLT=34;n_mds=20;
pca_rank = 2; pca_b_mlt = 44;
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));
addpath('/data/rangan/dir_bcc/dir_code_022316');
addpath('/data/rangan/dir_bcc/dir_PGC_20180304/dir_m');
dir_trunk = '/data/rangan/dir_bcc/dir_PGC_20180304';
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
lisa_setprefix_ver2 ;
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
lisa_setprefix_ver2 ;
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

niteration = 275;
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
%%%%%%%%%%%%%%%%;
disp(sprintf(' %% niteration %d',niteration));
%%%%%%%%;
pca_infix_r0_base = sprintf('pca_ni%d_tst%d',niteration,cl_num_arm1);
V_r0_base_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_V_.mda',dir_out_s0000_pca_r0_arm1,pca_infix_r0_base,pca_rank,pca_b_mlt)); %<-- This is the V-vector associated with (all) the snps in arm1, ordered to apply to arm1. ;
pca_infix_r0_arm1 = sprintf('pca_ni%d_tst%d',niteration,cl_num_arm2);
V_r0_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_V_.mda',dir_out_s0000_pca_r0_arm1,pca_infix_r0_arm1,pca_rank,pca_b_mlt)); %<-- This is the V-vector associated with the snps in both arm1 and arm2, ordered to apply to arm1. ;
V_r0_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_V_arm2_.mda',dir_out_s0000_pca_r0_arm1,pca_infix_r0_arm1,pca_rank,pca_b_mlt)); %<-- This is the V-vector associated with the snps in both arm1 and arm2, ordered to apply to arm2. ;
% Note that V_r0_arm1_ and V_r0_arm2_ are simply permutations of one another. ;
% While V_r0_base_ is not a permutation of V_r0_arm1_, these two vectors are quite correlated. ; plot(V_r0_base_(:,1),V_r0_arm1_(:,1),'.') ;
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
% establish snp-indices for each allelic combination. ;
%%%%%%%%;
b0_ = bim_arm1_{2}; allele1_ = bim_arm1_{5}; allele2_ = bim_arm1_{6}; alleletype_ = bim_arm1_{7}; atfr_base_ = bim_arm1_{9};
[b1_,b0_to_b1_,b1_to_b0_] = unique(b0_);
b0_by_b1_xref_ = sparse(1:length(b0_),b1_to_b0_,1,length(b0_),length(b1_));
%%%%%%%%;
flag_check=0;
if flag_check;
disp(sprintf(' %% checking b0_by_b1_xref_'));
nb1=max(1,min(length(b1_),floor(length(b1_)*rand())));
disp(sprintf(' %% nb1 %.9d: b1_{nb1} %s',nb1,b1_{nb1}));
tmp_b0_ = find(sum(b0_by_b1_xref_(:,nb1),2));
for nl0=1:length(tmp_b0_);
nb0 = tmp_b0_(nl0);
disp(sprintf(' %% nb0 %.9d: b0_{nb0} %s (%s)',nb0,b0_{nb0},alleletype_{nb0}));
assert(strcmp(b0_{nb0},b1_{nb1}));
end;%for nl0=1:length(tmp_b0_);
end;%if flag_check;
%%%%%%%%;
c0_ = find(strcmp(alleletype_,'and'));
[c1_,c0_to_c1_,c1_to_c0_] = unique(b0_(c0_));
c0_by_c1_xref_ = sparse(1:length(c0_),c1_to_c0_,1,length(c0_),length(c1_));
[d1_,d1_to_b1_,d1_to_c1_] = intersect(b1_,c1_,'stable');
b1_by_d1_xref_ = sparse(d1_to_b1_,1:length(d1_),1,length(b1_),length(d1_));
c1_by_d1_xref_ = sparse(d1_to_c1_,1:length(d1_),1,length(c1_),length(d1_));
V_r0_base_and_ = zeros(length(b1_),2);
atfr_base_and_ = zeros(length(b1_),1);
for nd1=1:length(d1_);
nb1 = find(sum(b1_by_d1_xref_(:,nd1),2));
nc1 = find(sum(c1_by_d1_xref_(:,nd1),2));
nc0 = find(sum(c0_by_c1_xref_(:,nc1),2));
nb0 = c0_(nc0);
V_r0_base_and_(nb1,:) = V_r0_base_(nb0,:);
atfr_base_and_(nb1) = atfr_base_(nb0);
end;%for nd1=1:length(d1_);
flag_check=0;
if flag_check;
disp(sprintf(' %% checking V_r0_base_and_'));
nb1=max(1,min(length(b1_),floor(length(b1_)*rand())));
tmp_snp = b1_{nb1};
nb0_and = find(strcmp(b0_,tmp_snp) & strcmp(alleletype_,'and'));
if ~isempty(nb0_and);
assert(length(intersect(nb0_and,find(sum(b0_by_b1_xref_(:,nb1),2))))>0);
disp(sprintf(' %% tmp_snp %s <-- %s (%s)',tmp_snp,b0_{nb0_and},alleletype_{nb0_and}));
disp(sprintf(' %% V_r0_base_ (%0.6f,%0.6f) , V_r0_base_and_ (%0.6f,%0.6f)',V_r0_base_(nb0_and,:),V_r0_base_and_(nb1,:)));
end;%if ~isempty(nb0_and);
end;%if flag_check;
clear c0_ c1_ c0_to_c1_ c1_to_c0_ d1_ d1_to_b1_ d1_to_c1_ b1_by_d1_xref_ c1_by_d1_xref_ ;
%%%%%%%%;
c0_ = find(strcmp(alleletype_,'xor'));
[c1_,c0_to_c1_,c1_to_c0_] = unique(b0_(c0_));
c0_by_c1_xref_ = sparse(1:length(c0_),c1_to_c0_,1,length(c0_),length(c1_));
[d1_,d1_to_b1_,d1_to_c1_] = intersect(b1_,c1_,'stable');
b1_by_d1_xref_ = sparse(d1_to_b1_,1:length(d1_),1,length(b1_),length(d1_));
c1_by_d1_xref_ = sparse(d1_to_c1_,1:length(d1_),1,length(c1_),length(d1_));
V_r0_base_xor_ = zeros(length(b1_),2);
atfr_base_xor_ = zeros(length(b1_),1);
for nd1=1:length(d1_);
nb1 = find(sum(b1_by_d1_xref_(:,nd1),2));
nc1 = find(sum(c1_by_d1_xref_(:,nd1),2));
nc0 = find(sum(c0_by_c1_xref_(:,nc1),2));
nb0 = c0_(nc0);
V_r0_base_xor_(nb1,:) = V_r0_base_(nb0,:);
atfr_base_xor_(nb1) = atfr_base_(nb0);
end;%for nd1=1:length(d1_);
flag_check=0;
if flag_check;
disp(sprintf(' %% checking V_r0_base_xor_'));
nb1=max(1,min(length(b1_),floor(length(b1_)*rand())));
tmp_snp = b1_{nb1};
nb0_xor = find(strcmp(b0_,tmp_snp) & strcmp(alleletype_,'xor'));
if ~isempty(nb0_xor);
assert(length(intersect(nb0_xor,find(sum(b0_by_b1_xref_(:,nb1),2))))>0);
disp(sprintf(' %% tmp_snp %s <-- %s (%s)',tmp_snp,b0_{nb0_xor},alleletype_{nb0_xor}));
disp(sprintf(' %% V_r0_base_ (%0.6f,%0.6f) , V_r0_base_xor_ (%0.6f,%0.6f)',V_r0_base_(nb0_xor,:),V_r0_base_xor_(nb1,:)));
end;%if ~isempty(nb0_xor);
end;%if flag_check;
clear c0_ c1_ c0_to_c1_ c1_to_c0_ d1_ d1_to_b1_ d1_to_c1_ b1_by_d1_xref_ c1_by_d1_xref_ ;
%%%%%%%%;
c0_ = find(strcmp(alleletype_,'nor'));
[c1_,c0_to_c1_,c1_to_c0_] = unique(b0_(c0_));
c0_by_c1_xref_ = sparse(1:length(c0_),c1_to_c0_,1,length(c0_),length(c1_));
[d1_,d1_to_b1_,d1_to_c1_] = intersect(b1_,c1_,'stable');
b1_by_d1_xref_ = sparse(d1_to_b1_,1:length(d1_),1,length(b1_),length(d1_));
c1_by_d1_xref_ = sparse(d1_to_c1_,1:length(d1_),1,length(c1_),length(d1_));
V_r0_base_nor_ = zeros(length(b1_),2);
atfr_base_nor_ = zeros(length(b1_),1);
for nd1=1:length(d1_);
nb1 = find(sum(b1_by_d1_xref_(:,nd1),2));
nc1 = find(sum(c1_by_d1_xref_(:,nd1),2));
nc0 = find(sum(c0_by_c1_xref_(:,nc1),2));
nb0 = c0_(nc0);
V_r0_base_nor_(nb1,:) = V_r0_base_(nb0,:);
atfr_base_nor_(nb1) = atfr_base_(nb0);
end;%for nd1=1:length(d1_);
flag_check=0;
if flag_check;
disp(sprintf(' %% checking V_r0_base_nor_'));
nb1=max(1,min(length(b1_),floor(length(b1_)*rand())));
tmp_snp = b1_{nb1};
nb0_nor = find(strcmp(b0_,tmp_snp) & strcmp(alleletype_,'nor'));
if ~isempty(nb0_nor);
assert(length(intersect(nb0_nor,find(sum(b0_by_b1_xref_(:,nb1),2))))>0);
disp(sprintf(' %% tmp_snp %s <-- %s (%s)',tmp_snp,b0_{nb0_nor},alleletype_{nb0_nor}));
disp(sprintf(' %% V_r0_base_ (%0.6f,%0.6f) , V_r0_base_nor_ (%0.6f,%0.6f)',V_r0_base_(nb0_nor,:),V_r0_base_nor_(nb1,:)));
end;%if ~isempty(nb0_nor);
end;%if flag_check;
clear c0_ c1_ c0_to_c1_ c1_to_c0_ d1_ d1_to_b1_ d1_to_c1_ b1_by_d1_xref_ c1_by_d1_xref_ ;
%%%%%%%%;
% Plotting allelic-types against one another. ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
tmp_ij_ = find( ~ ( (V_r0_base_and_(:,1)==0) & (V_r0_base_xor_(:,1)==0) & (V_r0_base_nor_(:,1)==0) ) );
plot3(V_r0_base_and_(tmp_ij_,1),V_r0_base_xor_(tmp_ij_,1),V_r0_base_nor_(tmp_ij_,1),'.');
xlabel('and'); ylabel('xor'); zlabel('nor'); axis vis3d;
end;% if flag_plot;
flag_plot=0;
if flag_plot;
tmp_ij_ = find( ~ ( (V_r0_base_and_(:,1)==0) & (V_r0_base_xor_(:,1)==0) & (V_r0_base_nor_(:,1)==0) ) );
subplot(3,3,1); plot(V_r0_base_and_(tmp_ij_,1),V_r0_base_xor_(tmp_ij_,1),'b.'); xlabel('and'); ylabel('xor'); title('and-vs-xor 1');
subplot(3,3,2); plot(V_r0_base_and_(tmp_ij_,1),V_r0_base_nor_(tmp_ij_,1),'b.'); xlabel('and'); ylabel('nor'); title('and-vs-nor 1');
subplot(3,3,3); plot(V_r0_base_xor_(tmp_ij_,1),V_r0_base_nor_(tmp_ij_,1),'b.'); xlabel('xor'); ylabel('nor'); title('xor-vs-nor 1');
subplot(3,3,4); plot(V_r0_base_and_(tmp_ij_,2),V_r0_base_xor_(tmp_ij_,2),'r.'); xlabel('and'); ylabel('xor'); title('and-vs-xor 2');
subplot(3,3,5); plot(V_r0_base_and_(tmp_ij_,2),V_r0_base_nor_(tmp_ij_,2),'r.'); xlabel('and'); ylabel('nor'); title('and-vs-nor 2');
subplot(3,3,6); plot(V_r0_base_xor_(tmp_ij_,2),V_r0_base_nor_(tmp_ij_,2),'r.'); xlabel('xor'); ylabel('nor'); title('xor-vs-nor 2');
subplot(3,3,7); plot(V_r0_base_and_(tmp_ij_,1),V_r0_base_and_(tmp_ij_,2),'g.'); xlabel('and 1'); ylabel('and 2'); title('and 1-vs-2');
subplot(3,3,8); plot(V_r0_base_xor_(tmp_ij_,1),V_r0_base_xor_(tmp_ij_,2),'g.'); xlabel('xor 1'); ylabel('xor 2'); title('xor 1-vs-2');
subplot(3,3,9); plot(V_r0_base_nor_(tmp_ij_,1),V_r0_base_nor_(tmp_ij_,2),'g.'); xlabel('nor 1'); ylabel('nor 2'); title('nor 1-vs-2');
end;% if flag_plot;
flag_plot=0;
if flag_plot;
tmp_ij_ = find( ~ ( (V_r0_base_and_(:,1)==0) & (V_r0_base_xor_(:,1)==0) & (V_r0_base_nor_(:,1)==0) ) );
subplot(3,1,1); plot(atfr_base_and_(tmp_ij_,1),atfr_base_xor_(tmp_ij_,1),'b.'); xlabel('and'); ylabel('xor'); title('and-vs-xor');
subplot(3,1,2); plot(atfr_base_and_(tmp_ij_,1),atfr_base_nor_(tmp_ij_,1),'b.'); xlabel('and'); ylabel('nor'); title('and-vs-nor');
subplot(3,1,3); plot(atfr_base_xor_(tmp_ij_,1),atfr_base_nor_(tmp_ij_,1),'b.'); xlabel('xor'); ylabel('nor'); title('xor-vs-nor');
end;% if flag_plot;

dosage_DnV_arm2_ = zeros(size(fam_arm2_{1},1),2);
dosage_XnV_arm2_ = zeros(size(fam_arm2_{1},1),2);

%for nstudy=1:length(study_name_arm2_);
% test file. ;
nstudy = 2;
study_name_arm2 = study_name_arm2_{nstudy};
tmp_ij = strfind(study_name_arm2,'bip_');
study_name_arm2_prefix = study_name_arm2(tmp_ij+(4:7));
fname_dosage = sprintf('%s/dir_misc/%s_eur_sr-qc.hg19.ch.fl.out.dosage_%s_out_cdrop_ni%d.out.dosage',dir_trunk,study_name_arm2_prefix,string_name_r0_arm1,niteration);
%fname_dosage = sprintf('%s/dir_misc/%s_eur_sr-qc.hg19.ch.fl.out.dosage_%s_out_cdrop_ni%d.out.dosage.cropped',dir_trunk,study_name_arm2_prefix,string_name_r0_arm1,niteration);
fcheck(fname_dosage);
%tic;[tmp_data_,tmp_varnames_,tmp_casenames_] = tblread(fname_dosage);toc;
n_dosage_snp = wc_0(fname_dosage)-1;

fid = fopen(fname_dosage);
tmp_line_ = fgetl(fid);
tmp_dosage_patient_ = strsplit(tmp_line_);
assert(strcmp(tmp_dosage_patient_{1},'SNP'));
assert(strcmp(tmp_dosage_patient_{2},'A1'));
assert(strcmp(tmp_dosage_patient_{3},'A2'));
n_dosage_patient = floor((length(tmp_dosage_patient_)-3)/2);
dosage_fam_name_ = cell(n_dosage_patient,1);
for ndosage_patient=1:n_dosage_patient;
dosage_fam_name_{ndosage_patient} = sprintf('%s&%s',strtrim(tmp_dosage_patient_{2 + 2*ndosage_patient}),strtrim(tmp_dosage_patient_{3 + 2*ndosage_patient}));
end;%for ndosage_patient=1:n_dosage_patient;
dosage_bim_snp_ = cell(n_dosage_snp,1);
dosage_bim_al1_ = cell(n_dosage_snp,1);
dosage_bim_al2_ = cell(n_dosage_snp,1);
for ndosage_snp=1:n_dosage_snp;
if (mod(ndosage_snp,1000)==0); disp(sprintf(' %% ndosage_snp %d/%d',ndosage_snp,n_dosage_snp)); end;
tmp_line_ = fgetl(fid);
tmp_dosage_snp_ = textscan(tmp_line_,'%s',3);%tmp_dosage_snp_ = strsplit(tmp_line_);
tmp_dosage_snp_ = tmp_dosage_snp_{1};
flag_check=0;
if flag_check;
tmp_ij_ = find(strcmp(bim_arm1_{2},tmp_dosage_snp_{1}));
for nl=1:length(tmp_ij_);
disp(sprintf(' %% found snp %s (%s,%s) <-- bim %s (%c,%c): %s',tmp_dosage_snp_{1},tmp_dosage_snp_{2},tmp_dosage_snp_{3},bim_arm1_{2}{tmp_ij_(nl)},bim_arm1_{5}(tmp_ij_(nl)),bim_arm1_{6}(tmp_ij_(nl)),alleletype_{tmp_ij_(nl)}));
end;%for nl=1:length(tmp_ij_);
end;%if flag_check;
dosage_bim_snp_{ndosage_snp} = tmp_dosage_snp_{1};
dosage_bim_al1_{ndosage_snp} = tmp_dosage_snp_{2};
dosage_bim_al2_{ndosage_snp} = tmp_dosage_snp_{3};
end;%for ndosage_snp=1:n_dosage_snp;
fclose(fid);

%%%%%%%%;
d1_ = dosage_bim_snp_ ; %<-- we assume this is already unique. ;
[e1_,e1_to_d1_,e1_to_b1_] = intersect(d1_,b1_,'stable');
d1_by_e1_xref_ = sparse(e1_to_d1_,1:length(e1_),1,length(d1_),length(e1_));
e1_by_d1_xref_ = sparse(1:length(e1_),e1_to_d1_,1,length(e1_),length(d1_));
b1_by_e1_xref_ = sparse(e1_to_b1_,1:length(e1_),1,length(b1_),length(e1_));
%%%%%%%%;
flag_check=0;
if flag_check;
ne1=max(1,min(length(e1_),floor(length(e1_)*rand())));
disp(sprintf(' %% ne1 %.9d: e1_{ne1} %s',ne1,e1_{ne1}));
tmp_d1_ = find(sum(d1_by_e1_xref_(:,ne1),2));
for nl1=1:length(tmp_d1_);
nd1 = tmp_d1_(nl1);
disp(sprintf(' %% nd1 %.9d: d1_{nd1} %s',nd1,d1_{nd1}));
assert(strcmp(d1_{nd1},e1_{ne1}));
end;%for nl1=1:length(tmp_d1_);
tmp_b1_ = find(sum(b1_by_e1_xref_(:,ne1),2));
for nl1=1:length(tmp_b1_);
nb1 = tmp_b1_(nl1);
disp(sprintf(' %% nb1 %.9d: b1_{nb1} %s and (%0.6f %0.6f) , xor (%0.6f %0.6f) , nor (%0.6f %0.6f)',nb1,b1_{nb1},V_r0_base_and_(nb1,:),V_r0_base_xor_(nb1,:),V_r0_base_nor_(nb1,:)));
assert(strcmp(b1_{nb1},e1_{ne1}));
end;%for nl1=1:length(tmp_b1_);
tmp_b0_ = find(sum(b0_by_b1_xref_(:,nb1),2));
for nl1=1:length(tmp_b0_);
nb0 = tmp_b0_(nl1);
disp(sprintf(' %% nb0 %.9d: b0_{nb0} %s (%s) (%0.6f %0.6f)',nb0,b0_{nb0},alleletype_{nb0},V_r0_base_(nb0,:)));
assert(strcmp(b0_{nb0},e1_{ne1}));
end;%for nl0=0:length(tmp_b0_);
end;%if flag_check;

% Ordering V_r0_base_xxx_ to align with dosage file. ;
V_r0_perm_and_ = zeros(length(d1_),2);
V_r0_perm_xor_ = zeros(length(d1_),2);
V_r0_perm_nor_ = zeros(length(d1_),2);
atfr_perm_and_ = zeros(length(d1_),1);
atfr_perm_xor_ = zeros(length(d1_),1);
atfr_perm_nor_ = zeros(length(d1_),1);
for ne1=1:length(e1_);
tmp_d1_ = find(sum(d1_by_e1_xref_(:,ne1),2));
if (length(tmp_d1_)>1); disp(sprintf(' %% Warning! tmp_d1_ length %d in lisa_dosage_pca_ver1',length(tmp_d1_))); end;
nd1 = tmp_d1_(1);
tmp_b1_ = find(sum(b1_by_e1_xref_(:,ne1),2));
if (length(tmp_b1_)>1); disp(sprintf(' %% Warning! tmp_b1_ length %d in lisa_dosage_pca_ver1',length(tmp_b1_))); end;
nb1 = tmp_b1_(1);
V_r0_perm_and_(nd1,:) = V_r0_base_and_(nb1,:);
V_r0_perm_xor_(nd1,:) = V_r0_base_xor_(nb1,:);
V_r0_perm_nor_(nd1,:) = V_r0_base_nor_(nb1,:);
atfr_perm_and_(nd1) = atfr_base_and_(nb1);
atfr_perm_xor_(nd1) = atfr_base_xor_(nb1);
atfr_perm_nor_(nd1) = atfr_base_nor_(nb1);
end;%for ne1=1:length(e1_);
atfr_perm_and_alpha_ = atfr_perm_and_ - (1 - atfr_perm_and_); atfr_perm_and_delta_ = 1./( 4 .* atfr_perm_and_ .* (1 - atfr_perm_and_) );
atfr_perm_xor_alpha_ = atfr_perm_xor_ - (1 - atfr_perm_xor_); atfr_perm_xor_delta_ = 1./( 4 .* atfr_perm_xor_ .* (1 - atfr_perm_xor_) );
atfr_perm_nor_alpha_ = atfr_perm_nor_ - (1 - atfr_perm_nor_); atfr_perm_nor_delta_ = 1./( 4 .* atfr_perm_nor_ .* (1 - atfr_perm_nor_) );

f1_ = fam_name_arm2_; %<-- We assume this is unique. ;
g1_ = dosage_fam_name_; %<-- We assume this is unique. ;
[h1_,h1_to_g1_,h1_to_f1_] = intersect(g1_,f1_,'stable');
g1_by_h1_xref_ = sparse(h1_to_g1_,1:length(h1_),1,length(g1_),length(h1_));
f1_by_h1_xref_ = sparse(h1_to_f1_,1:length(h1_),1,length(f1_),length(h1_));
g1_to_f1_ = zeros(length(g1_),1);
for nh1=1:length(h1_);
tmp_f1_ = find(sum(f1_by_h1_xref_(:,nh1),2));
if (length(tmp_f1_)>1); disp(sprintf(' %% Warning! tmp_f1_ length %d in lisa_dosage_pca_ver1',length(tmp_f1_))); end;
nf1 = tmp_f1_(1);
tmp_g1_ = find(sum(g1_by_h1_xref_(:,nh1),2));
if (length(tmp_g1_)>1); disp(sprintf(' %% Warning! tmp_g1_ length %d in lisa_dosage_pca_ver1',length(tmp_g1_))); end;
ng1 = tmp_g1_(1);
assert(strcmp(f1_{nf1},g1_{ng1}));
g1_to_f1_(ng1) = nf1;
end;%for nh1=1:length(h1_);

%%%%%%%%%%%%%%%%;
% Now step through the dosage file line by line. ;
% updating the inner-products as we go. ;
%%%%%%%%%%%%%%%%;
tmp_AnV_dosage_ = zeros(n_dosage_patient,2);
n_mismatch_allele1 = 0;
n_mismatch_allele2 = 0;
fid = fopen(fname_dosage);
tmp_line_ = fgetl(fid); %<-- skip header line. ;
for ndosage_snp=1:n_dosage_snp;
if (mod(ndosage_snp,1000)==0); disp(sprintf(' %% ndosage_snp %d/%d',ndosage_snp,n_dosage_snp)); end;
tmp_line_ = fgetl(fid);
tmp_dosage_snp_ = textscan(tmp_line_,'%s',3);
tmp_dosage_al1_ = tmp_dosage_snp_{1}{2};
tmp_dosage_al2_ = tmp_dosage_snp_{1}{3};
nd1 = ndosage_snp;
ne1 = find(sum(e1_by_d1_xref_(:,nd1),2)); assert(length(ne1)<=1);
if (length(ne1)==1);
nb1 = find(sum(b1_by_e1_xref_(:,ne1),2)); assert(length(nb1)==1);
tmp_b0_ = find(sum(b0_by_b1_xref_(:,nb1),2)); assert(length(tmp_b0_)>=1);
nb0 = tmp_b0_(1);
assert(strcmp(b0_{nb0},d1_{nd1}));
assert(strcmp(b0_{nb0},e1_{ne1}));
assert(strcmp(b0_{nb0},b1_{nb1}));
assert(strcmp(b0_{nb0},tmp_dosage_snp_{1}{1}));
flag_allele1 = strcmp(tmp_dosage_al1_,allele1_(nb0));
flag_allele2 = strcmp(tmp_dosage_al2_,allele2_(nb0));
if (flag_allele1~=1); n_mismatch_allele1 = n_mismatch_allele1 + 1; end;
if (flag_allele2~=1); n_mismatch_allele2 = n_mismatch_allele1 + 2; end;
tmp_length = sum(cellfun('length',tmp_dosage_snp_{1})) + length(tmp_dosage_snp_{1});
tmp_val_ = sscanf(tmp_line_(tmp_length:end),'%f');
assert(length(tmp_val_)/2==n_dosage_patient);
tmp_val_and_ = tmp_val_(1:2:end-1);
tmp_val_xor_ = tmp_val_(2:2:end-0);
tmp_val_nor_ = ones(n_dosage_patient,1) - tmp_val_and_ - tmp_val_xor_;
tmp_val_and_normalized_ = (2*tmp_val_and_ - 1 - atfr_perm_and_alpha_(nd1)) / sqrt(atfr_perm_and_delta_(nd1));
tmp_val_xor_normalized_ = (2*tmp_val_xor_ - 1 - atfr_perm_xor_alpha_(nd1)) / sqrt(atfr_perm_xor_delta_(nd1));
tmp_val_nor_normalized_ = (2*tmp_val_nor_ - 1 - atfr_perm_nor_alpha_(nd1)) / sqrt(atfr_perm_nor_delta_(nd1));

%%%%%%%%%%%%%%%%;
% check against b16 file. ;
%%%%%%%%%%%%%%%%;
flag_check=0;
if flag_check;
%%%%%%%%;
tmp_ij_ = find(strcmp(bim_arm2_{2},tmp_dosage_snp_{1}{1}));
if (length(tmp_ij_)>0);
tmp_ij= tmp_ij_(1);
disp(sprintf(' %% snp %s ; dosage [%s %s] , arm1 [%s %s] (%s), arm2 [%s %s] (%s)',tmp_dosage_snp_{1}{1},tmp_dosage_al1_,tmp_dosage_al2_,allele1_(nb0),allele2_(nb0),alleletype_{nb0},bim_arm2_{5}(tmp_ij),bim_arm2_{6}(tmp_ij),bim_arm2_{7}{tmp_ij}));
fname_A_n_arm2 = sprintf('%s/%s_A_full_n.b16',dir__in_arm2,string_prefix_arm2); fcheck(fname_A_n_arm2);
[tmp_nbins,tmp_nrows,tmp_ncols] = tutorial_binary_getsize(fname_A_n_arm2);
%%%%%%%%;
tmp_ij_and = find(strcmp(bim_arm2_{2},tmp_dosage_snp_{1}{1}) & strcmp(bim_arm2_{7},'and')); if (length(tmp_ij_and)>1); disp(sprintf(' %% Warning! tmp_ij_and length %d in lisa_dosage_pca_ver2.m',length(tmp_ij_and))); end;
if (length(tmp_ij_and)>0);
tmp_arm2_and_ = tutorial_binary_uncompress(fname_A_n_arm2,1:tmp_nrows,min(tmp_ncols,tmp_ij_and));
tmp_arm2_perm_and_ = zeros(length(g1_),1);
for nh1=1:length(h1_);
tmp_f1_ = find(sum(f1_by_h1_xref_(:,nh1),2)); if (length(tmp_f1_)>1); disp(sprintf(' %% Warning! tmp_f1_ length %d in lisa_dosage_pca_ver1',length(tmp_f1_))); end; nf1 = tmp_f1_(1);
tmp_g1_ = find(sum(g1_by_h1_xref_(:,nh1),2)); if (length(tmp_g1_)>1); disp(sprintf(' %% Warning! tmp_g1_ length %d in lisa_dosage_pca_ver1',length(tmp_g1_))); end; ng1 = tmp_g1_(1);
assert(strcmp(f1_{nf1},g1_{ng1}));
tmp_arm2_perm_and_(ng1) = tmp_arm2_and_(nf1);
end;%for nh1=1:length(h1_);
disp(sprintf(' %% nd1 %d tmp_ij_and %d %s (%s) arm2_and error %0.16f',nd1,tmp_ij_and,tmp_dosage_snp_{1}{1},bim_arm2_{7}{tmp_ij_and},norm((tmp_arm2_perm_and_>0) - tmp_val_and_)));
end;%if (length(tmp_ij_and)>0);
%%%%%%%%;
tmp_ij_xor = find(strcmp(bim_arm2_{2},tmp_dosage_snp_{1}{1}) & strcmp(bim_arm2_{7},'xor')); if (length(tmp_ij_xor)>1); disp(sprintf(' %% Warning! tmp_ij_xor length %d in lisa_dosage_pca_ver2.m',length(tmp_ij_xor))); end;
if (length(tmp_ij_xor)>0);
tmp_arm2_xor_ = tutorial_binary_uncompress(fname_A_n_arm2,1:tmp_nrows,min(tmp_ncols,tmp_ij_xor));
tmp_arm2_perm_xor_ = zeros(length(g1_),1);
for nh1=1:length(h1_);
tmp_f1_ = find(sum(f1_by_h1_xref_(:,nh1),2)); if (length(tmp_f1_)>1); disp(sprintf(' %% Warning! tmp_f1_ length %d in lisa_dosage_pca_ver1',length(tmp_f1_))); end; nf1 = tmp_f1_(1);
tmp_g1_ = find(sum(g1_by_h1_xref_(:,nh1),2)); if (length(tmp_g1_)>1); disp(sprintf(' %% Warning! tmp_g1_ length %d in lisa_dosage_pca_ver1',length(tmp_g1_))); end; ng1 = tmp_g1_(1);
assert(strcmp(f1_{nf1},g1_{ng1}));
tmp_arm2_perm_xor_(ng1) = tmp_arm2_xor_(nf1);
end;%for nh1=1:length(h1_);
disp(sprintf(' %% nd1 %d tmp_ij_xor %d %s (%s) arm2_xor error %0.16f',nd1,tmp_ij_xor,tmp_dosage_snp_{1}{1},bim_arm2_{7}{tmp_ij_xor},norm((tmp_arm2_perm_xor_>0) - tmp_val_xor_)));
end;%if (length(tmp_ij_xor)>0);
%%%%%%%%;
tmp_ij_nor = find(strcmp(bim_arm2_{2},tmp_dosage_snp_{1}{1}) & strcmp(bim_arm2_{7},'nor')); if (length(tmp_ij_nor)>1); disp(sprintf(' %% Warning! tmp_ij_nor length %d in lisa_dosage_pca_ver2.m',length(tmp_ij_nor))); end;
if (length(tmp_ij_nor)>0);
tmp_arm2_nor_ = tutorial_binary_uncompress(fname_A_n_arm2,1:tmp_nrows,min(tmp_ncols,tmp_ij_nor));
tmp_arm2_perm_nor_ = zeros(length(g1_),1);
for nh1=1:length(h1_);
tmp_f1_ = find(sum(f1_by_h1_xref_(:,nh1),2)); if (length(tmp_f1_)>1); disp(sprintf(' %% Warning! tmp_f1_ length %d in lisa_dosage_pca_ver1',length(tmp_f1_))); end; nf1 = tmp_f1_(1);
tmp_g1_ = find(sum(g1_by_h1_xref_(:,nh1),2)); if (length(tmp_g1_)>1); disp(sprintf(' %% Warning! tmp_g1_ length %d in lisa_dosage_pca_ver1',length(tmp_g1_))); end; ng1 = tmp_g1_(1);
assert(strcmp(f1_{nf1},g1_{ng1}));
tmp_arm2_perm_nor_(ng1) = tmp_arm2_nor_(nf1);
end;%for nh1=1:length(h1_);
disp(sprintf(' %% nd1 %d tmp_ij_nor %d %s (%s) arm2_nor error %0.16f',nd1,tmp_ij_nor,tmp_dosage_snp_{1}{1},bim_arm2_{7}{tmp_ij_nor},norm((tmp_arm2_perm_nor_>0) - tmp_val_nor_)));
end;%if (length(tmp_ij_nor)>0);
%%%%%%%%;
end;%if (length(tmp_ij_)>0);
end;%if flag_check;

tmp_AnV_dosage_ = tmp_AnV_dosage_ ...
   + tmp_val_and_normalized_*V_r0_perm_and_(nd1,:) ...
   + tmp_val_xor_normalized_*V_r0_perm_xor_(nd1,:) ...
   + tmp_val_nor_normalized_*V_r0_perm_nor_(nd1,:) ...
  ;
end;%if (length(ne1)==1);
end;%for ndosage_snp=1:n_dosage_snp;
fclose(fid);

tmp_AnV_arm2_ = zeros(length(f1_),2);
for nh1=1:length(h1_);
tmp_f1_ = find(sum(f1_by_h1_xref_(:,nh1),2));
if (length(tmp_f1_)>1); disp(sprintf(' %% Warning! tmp_f1_ length %d in lisa_dosage_pca_ver1',length(tmp_f1_))); end;
nf1 = tmp_f1_(1);
tmp_g1_ = find(sum(g1_by_h1_xref_(:,nh1),2));
if (length(tmp_g1_)>1); disp(sprintf(' %% Warning! tmp_g1_ length %d in lisa_dosage_pca_ver1',length(tmp_g1_))); end;
ng1 = tmp_g1_(1);
assert(strcmp(f1_{nf1},g1_{ng1}));
tmp_AnV_arm2_(nf1,:) = tmp_AnV_dosage_(ng1,:);
end;%for nh1=1:length(h1_);

tmp_DnV_arm2_ = zeros(length(f1_),2); tmp_Dij_ = find(mr_D_arm2_); tmp_DnV_arm2_(tmp_Dij_,:) = tmp_AnV_arm2_(tmp_Dij_,:);
tmp_XnV_arm2_ = zeros(length(f1_),2); tmp_Xij_ = find(mr_X_arm2_); tmp_XnV_arm2_(tmp_Xij_,:) = tmp_AnV_arm2_(tmp_Xij_,:);

flag_plot=1;
if flag_plot;
subplot(2,2,1);
tmp_Dij_ = intersect(g1_to_f1_,find(mr_D_arm2_));
plot(DnV_r0_arm2_(tmp_Dij_,1),tmp_DnV_arm2_(tmp_Dij_,1),'.');
title('D 1');
subplot(2,2,2);
tmp_Xij_ = intersect(g1_to_f1_,find(mr_X_arm2_));
plot(XnV_r0_arm2_(tmp_Xij_,1),tmp_XnV_arm2_(tmp_Xij_,1),'.');
title('X 1');
subplot(2,2,3);
tmp_Dij_ = intersect(g1_to_f1_,find(mr_D_arm2_));
plot(DnV_r0_arm2_(tmp_Dij_,2),tmp_DnV_arm2_(tmp_Dij_,2),'.');
title('D 2');
subplot(2,2,4);
tmp_Xij_ = intersect(g1_to_f1_,find(mr_X_arm2_));
plot(XnV_r0_arm2_(tmp_Xij_,2),tmp_XnV_arm2_(tmp_Xij_,2),'.');
title('X 2');
end;%if flag_plot;

dosage_DnV_arm2_ = dosage_DnV_arm2_ + tmp_DnV_arm2_ ;
dosage_XnV_arm2_ = dosage_XnV_arm2_ + tmp_XnV_arm2_ ;

%end;%for nstudy=1:length(study_name_arm2_);