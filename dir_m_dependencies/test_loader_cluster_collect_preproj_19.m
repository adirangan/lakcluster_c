%function test_loader_cluster_collect_preproj_19();
% intended for use with test_loader_38.m ;
% loading some of the preproj data. ;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

verbose=1;
flag_recalc = 0;
flag_replot = 0;
tolerance_master = 1e-2;
nf=0;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root);
str_code = sprintf('%s/halfloop_dev',dir_code);
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jamison',string_root);
dir_data = sprintf('%s/data_summary_20190730',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_cluster = sprintf('%s/dir_loader_cluster',dir_trunk);

flag_compute = strcmp(platform,'access1');
flag_plot_make = strcmp(platform,'OptiPlex');

if flag_compute;

%%%%%%%%;
% checking clustering. ;
%%%%%%%%;
date_diff_threshold = 1.0;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;
fp_label_A_ = fopen(sprintf('%s/str_CLabel_sub_.nsv',dir_mat),'r');
str_label_A_ = textscan(fp_label_A_,'%s'); fclose(fp_label_A_);
n_u = numel(str_label_A_{1});
label_A_ = label_str_to_enum_0(str_label_A_{1});
u_label_A_ = unique(str_label_A_{1});
n_label_A = length(u_label_A_);
n_u_label_A_ = zeros(n_label_A,1);
label_A_each__ = zeros(n_u,n_label_A);
for nlabel_A=0:n_label_A-1;
label_A_each__(:,1+nlabel_A) = zeros(n_u,1);
tmp_index_ = efind(strcmp(str_label_A_{1},u_label_A_{1+nlabel_A}));
label_A_each__(1+tmp_index_,1+nlabel_A) = 1;
n_u_label_A_(1+nlabel_A) = numel(tmp_index_);
end;%for nlabel_A=0:n_label_A-1;
u_label_A_enum_ = zeros(n_label_A,1); 
u_label_A_enum_(1:n_label_A-1) = cellfun(@str2num,u_label_A_(1:n_label_A-1));
u_label_A_enum_(end) = n_label_A;
[~,index_label_A_] = sort(u_label_A_enum_,'ascend'); index_label_A_ = index_label_A_ - 1;
%%%%%%%%;
gamma = 0.01; n_shuffle = 64; p_set = 0.05; n_member_lob = 3;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d); str_g0 = sprintf('_g%.3d',0*gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_dex = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
str_preproj_dex = sprintf('preproj_dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g0,str_p,str_nml);
str_hnb = sprintf('hnbtZRgumb%s%s%s',str_g,str_p,str_nml);
str_preproj_hnb = sprintf('preproj_hnbtZRgumb%s%s%s',str_g0,str_p,str_nml);
str_preproj_hnbr0 = sprintf('preproj_hnbr0tZRgumb%s%s%s',str_g0,str_p,str_nml);
%%%%%%%%;
str_X_ = {'E','I'}; n_X = numel(str_X_);
prefix_normalization_ = {'logn','liXf','liXf_clogr','fill'}; n_normalization = numel(prefix_normalization_);
%prefix_normalization_ = {'logn'}; n_normalization = numel(prefix_normalization_);
prefix_covariate_ = {'C','G'}; n_covariate = numel(prefix_covariate_);
prefix_method_ = ...
{ ...
 'spectral_isosplit5' ...
,'tsne00_isosplit5' ...
,'tsne50_isosplit5' ...
,'umap00_default' ...
,'umap00_isosplit5' ...
,'umap00_hdbscan' ...
,'louvain00_default' ...
,str_dex ...
,str_hnb ...
,'preproj_spectral_isosplit5' ...
,'preproj_tsne00_isosplit5' ...
,'preproj_tsne50_isosplit5' ...
,'preproj_umap00_default' ...
,'preproj_umap00_isosplit5' ...
,'preproj_umap00_hdbscan' ...
,'preproj_louvain00_default' ...
,str_preproj_dex ...
,str_preproj_hnb ...
,str_preproj_hnbr0 ...
}; 
n_method = numel(prefix_method_);
%%%%%%%%;
n_correction___ = zeros(n_X,n_normalization,n_covariate);
prefix_correction___ = cell(n_X,n_normalization,n_covariate);
rank_estimate___ = cell(n_X,n_normalization,n_covariate);
for nX=0:n_X-1;
str_X = str_X_{1+nX};
for nnormalization=0:n_normalization-1;
prefix_normalization = prefix_normalization_{1+nnormalization};
for ncovariate=0:n_covariate-1;
prefix_covariate = prefix_covariate_{1+ncovariate};
prefix_correction_ = {''};
rank_estimate_ = [0];
n_correction = 1;
%%%%;
prefix_correction_{n_correction+1} = sprintf('%s_call',prefix_covariate);
tmp_fname_mat = sprintf('%s/S_call_%s_original.mat',dir_mat,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate = tmp_.rank_estimate_sample; clear tmp_; rank_estimate_(n_correction+1) = rank_estimate;
n_correction = n_correction+1;
%%%%;
tmp_fname_tsv = sprintf('%s/n_QCluster_%s_original.tsv',dir_mat,prefix_covariate);
n_QCluster_original = textread(tmp_fname_tsv); n_QCluster_original = n_QCluster_original(1);
if (verbose); disp(sprintf(' %% %s n_QCluster_original %d',prefix_covariate,n_QCluster_original)); end;
for nQCluster=0:n_QCluster_original-1;
prefix_correction_{n_correction+1} = sprintf('%s_c%d',prefix_covariate,1+nQCluster);
tmp_fname_mat = sprintf('%s/S_c%d_%s_original.mat',dir_mat,1+nQCluster,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate = tmp_.rank_estimate_sample; clear tmp_; rank_estimate_(n_correction+1) = rank_estimate;
n_correction = n_correction+1;
end;%for nQCluster=0:n_QCluster_original-1;
%%%%;
prefix_correction_{n_correction+1} = sprintf('UX_%s_call',prefix_covariate);
tmp_fname_mat = sprintf('%s/S_call_U%s_%s_%s_mahalanobis.mat',dir_mat,str_X,prefix_normalization,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate = tmp_.rank_estimate_sample; clear tmp_; rank_estimate_(n_correction+1) = rank_estimate;
n_correction = n_correction+1;
tmp_fname_tsv = sprintf('%s/n_QCluster_U%s_%s_%s_mahalanobis.tsv',dir_mat,str_X,prefix_normalization,prefix_covariate);
n_QCluster_mahalanobis = textread(tmp_fname_tsv); n_QCluster_mahalanobis = n_QCluster_mahalanobis(1);
if (verbose); disp(sprintf(' %% %s_%s_%s n_QCluster_mahalanobis %d',str_X,prefix_normalization,prefix_covariate,n_QCluster_mahalanobis)); end;
for nQCluster=0:n_QCluster_mahalanobis-1;
prefix_correction_{n_correction+1} = sprintf('UX_%s_c%d',prefix_covariate,1+nQCluster);
tmp_fname_mat = sprintf('%s/S_c%d_U%s_%s_%s_mahalanobis.mat',dir_mat,1+nQCluster,str_X,prefix_normalization,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate = tmp_.rank_estimate_sample; clear tmp_; rank_estimate_(n_correction+1) = rank_estimate;
n_correction = n_correction+1;
end;%for nQCluster=0:n_QCluster_mahalanobis-1;
%%%%;
n_correction___(1+nX,1+nnormalization,1+ncovariate) = n_correction;
rank_estimate___{1+nX,1+nnormalization,1+ncovariate} = rank_estimate_;
prefix_correction___{1+nX,1+nnormalization,1+ncovariate} = prefix_correction_;
end;%for ncovariate=0:n_covariate-1;
end;%for nnormalization=0:n_normalization-1;
end;%for nX=0:n_X-1;
%%%%%%%%;
fname_sub_meta_ = sprintf('%s/tlcc_preproj_meta_17.mat',dir_mat);
save(fname_sub_meta_);
%%%%%%%%;
for nX=0:n_X-1;
str_X = str_X_{1+nX};
for nnormalization=0:n_normalization-1;
prefix_normalization = prefix_normalization_{1+nnormalization};
for ncovariate=0:n_covariate-1;
prefix_covariate = prefix_covariate_{1+ncovariate};
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
prefix_correction_ = prefix_correction___{1+nX,1+nnormalization,1+ncovariate};
%%%%;
fname_sub_cm__ = sprintf('%s/tlcc_preproj_%s_%s_%s_17.mat',dir_mat,str_X,prefix_normalization,prefix_covariate);
if (verbose); disp(sprintf(' %% collecting fname_sub_cm__: %s',fname_sub_cm__)); end;
if ( exist(fname_sub_cm__,'file'));
if (verbose); disp(sprintf(' %% %s found, not creating',fname_sub_cm__)); end;
end;%if ( exist(fname_sub_cm__,'file'));
if (flag_recalc | ~exist(fname_sub_cm__,'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_sub_cm__)); end;
%%%%;
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
flag_stage0_sub_cm__ = zeros(n_correction,n_method);
flag_stage1_sub_cm__ = zeros(n_correction,n_method);
n_cluster_stage0_sub_cm__ = cell(n_correction,n_method);
lpv_A_vs_B_sub_cm__ = cell(n_correction,n_method);
lP0_A_vs_B_sub_cm__ = cell(n_correction,n_method);
n_X_GENE_sub_cm__ = zeros(n_correction,1);
markergene_A_sub_cm__ = cell(n_correction,1);
n_Z_index_sub_cm__ = zeros(n_correction,n_method);
lpv_A_each_vs_B_sub_cm__ = cell(n_correction,n_method);
nlpv_A_each_vs_B_each_sub_cm__ = cell(n_correction,n_method);
%%%%;
for ncorrection=0:n_correction-1;
prefix_correction = prefix_correction_{1+ncorrection};
if 0;
elseif (numel(prefix_correction)==0); tmp_dir_cluster = sprintf('%s/dir_%s_%s_cluster',dir_cluster,str_X,prefix_normalization); 
elseif (numel(prefix_correction)> 0); tmp_dir_cluster = sprintf('%s/dir_%s_%s_%s_cluster',dir_cluster,str_X,prefix_normalization,prefix_correction);
end;%if;
if (verbose); disp(sprintf(' %% dir: %s',tmp_dir_cluster)); end;
tmp_fname_mat = sprintf('%s/markergene_A__.mat',tmp_dir_cluster);
if (~exist(tmp_fname_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',tmp_fname_mat)); end;
end;%if (~exist(tmp_fname_mat,'file'));
if ( exist(tmp_fname_mat,'file'));
tmp_ = load(tmp_fname_mat);
n_X_GENE_sub_cm__(1+ncorrection) = tmp_.n_X_GENE;
markergene_A_sub_cm__{1+ncorrection} = tmp_.markergene_auz_A__;
clear tmp_;
end;%if ( exist(tmp_fname_mat,'file'));
for nmethod=0:n_method-1;
prefix_method = prefix_method_{1+nmethod};
tmp_fname_mat = sprintf('%s/%s.mat',tmp_dir_cluster,prefix_method);
if (~exist(tmp_fname_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',tmp_fname_mat)); end;
end;%if (~exist(tmp_fname_mat,'file'));
if ( exist(tmp_fname_mat,'file'));
flag_stage0_sub_cm__(1+ncorrection,1+nmethod) = 1;
tmp_ = load(tmp_fname_mat);
n_Z_index = numel(tmp_.label_B__);
n_cluster_stage0_sub_cm__{1+ncorrection,1+nmethod} = zeros(n_Z_index,1);
lpv_A_vs_B_sub_cm__{1+ncorrection,1+nmethod} = tmp_.lpv_;
lP0_A_vs_B_sub_cm__{1+ncorrection,1+nmethod} = tmp_.lP0_;
for nZ_index=0:n_Z_index-1;
n_cluster_stage0_sub_cm__{1+ncorrection,1+nmethod}(1+nZ_index) = numel(unique(tmp_.label_B__{1+nZ_index}));
end;%for nZ_index=0:n_Z_index-1;
clear tmp_;
end;%if ( exist(tmp_fname_mat,'file'));
tmp_fname_mat = sprintf('%s/%s_markergene_B___.mat',tmp_dir_cluster,prefix_method);
if (~exist(tmp_fname_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',tmp_fname_mat)); end;
end;%if (~exist(tmp_fname_mat,'file'));
if ( exist(tmp_fname_mat,'file'));
flag_stage1_sub_cm__(1+ncorrection,1+nmethod) = 1;
tmp_ = load(tmp_fname_mat);
n_Z_index_sub_cm__(1+ncorrection,1+nmethod) = tmp_.n_Z_index;
lpv_A_each_vs_B_sub_cm__{1+ncorrection,1+nmethod} = tmp_.lpv_A_each_vs_B__;
nlpv_A_each_vs_B_each_sub_cm__{1+ncorrection,1+nmethod} = tmp_.nlpv_each___;
clear tmp_;
end;%if ( exist(tmp_fname_mat,'file'));
end;%for nmethod=0:n_method-1;
end;%for ncorrection=0:n_correction-1;
save(fname_sub_cm__ ...
     ,'nX','str_X' ...
     ,'nnormalization','prefix_normalization' ...
     ,'ncovariate','prefix_covariate','n_correction','prefix_correction_' ...
     ,'flag_stage0_sub_cm__' ...
     ,'flag_stage1_sub_cm__' ...
     ,'n_cluster_stage0_sub_cm__' ...
     ,'lpv_A_vs_B_sub_cm__' ...
     ,'lP0_A_vs_B_sub_cm__' ...
     ,'n_X_GENE_sub_cm__' ...
     ,'markergene_A_sub_cm__' ...
     ,'n_Z_index_sub_cm__' ...
     ,'lpv_A_each_vs_B_sub_cm__' ...
     ,'nlpv_A_each_vs_B_each_sub_cm__' ...
     , '-v7.3' ...
     );
end;%if (~exist(fname_sub_cm__,'file'));
%%%%;
end;%for ncovariate=0:n_covariate-1;
end;%for nnormalization=0:n_normalization-1;
end;%for nX=0:n_X-1;
%%%%%%%%;

end;%if flag_compute;

if flag_plot_make;

%%%%%%%%;
% copy to OptiPlex.; 
%%%%%%%%;
% see test_loader_cluster_collect_17_log.sh;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

verbose=1;
flag_recalc = 0;
flag_replot = 1;
tolerance_master = 1e-2;
nf=0;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root);
str_code = sprintf('%s/halfloop_dev',dir_code);
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jamison',string_root);
dir_data = sprintf('%s/data_summary_20190730',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_cluster = sprintf('%s/dir_loader_cluster',dir_trunk);

%%%%%%%%;
% collect individual collections. ;
%%%%%%%%;
fname_sub_meta_ = sprintf('%s/tlcc_preproj_meta_17.mat',dir_mat);
load(fname_sub_meta_);
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); string_root = 'home'; end;
if (strcmp(platform,'eval1')); string_root = 'home'; end;
if (strcmp(platform,'rusty')); string_root = 'mnt/home'; end;
dir_code = sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root);
str_code = sprintf('%s/halfloop_dev',dir_code);
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jamison',string_root);
dir_data = sprintf('%s/data_summary_20190730',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_cluster = sprintf('%s/dir_loader_cluster',dir_trunk);
%%%%%%%%;
flag_stage0___ = cell(n_X,n_normalization,n_covariate);
flag_stage1___ = cell(n_X,n_normalization,n_covariate);
n_cluster_stage0___ = cell(n_X,n_normalization,n_covariate);
lpv_A_vs_B___ = cell(n_X,n_normalization,n_covariate);
lP0_A_vs_B___ = cell(n_X,n_normalization,n_covariate);
n_X_GENE___ = cell(n_X,n_normalization,n_covariate);
markergene_A___ = cell(n_X,n_normalization,n_covariate);
n_Z_index___ = cell(n_X,n_normalization,n_covariate);
lpv_A_each_vs_B___ = cell(n_X,n_normalization,n_covariate);
nlpv_A_each_vs_B_each___ = cell(n_X,n_normalization,n_covariate);
for nX=0:n_X-1;
for nnormalization=0:n_normalization-1;
for ncovariate=0:n_covariate-1;
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
flag_stage0___{1+nX,1+nnormalization,1+ncovariate} = zeros(n_correction,n_method);
flag_stage1___{1+nX,1+nnormalization,1+ncovariate} = zeros(n_correction,n_method);
n_cluster_stage0___{1+nX,1+nnormalization,1+ncovariate} = cell(n_correction,n_method);
lpv_A_vs_B___{1+nX,1+nnormalization,1+ncovariate} = cell(n_correction,n_method);
lP0_A_vs_B___{1+nX,1+nnormalization,1+ncovariate} = cell(n_correction,n_method);
n_X_GENE___{1+nX,1+nnormalization,1+ncovariate} = zeros(n_correction,1);
markergene_A___{1+nX,1+nnormalization,1+ncovariate} = cell(n_correction,1);
n_Z_index___{1+nX,1+nnormalization,1+ncovariate} = zeros(n_correction,n_method);
lpv_A_each_vs_B___{1+nX,1+nnormalization,1+ncovariate} = cell(n_correction,n_method);
nlpv_A_each_vs_B_each___{1+nX,1+nnormalization,1+ncovariate} = cell(n_correction,n_method);
end;%for ncovariate=0:n_covariate-1;
end;%for nnormalization=0:n_normalization-1;
end;%for nX=0:n_X-1;
%%%%%%%%;
for nX=0:n_X-1;
str_X = str_X_{1+nX};
for nnormalization=0:n_normalization-1;
prefix_normalization = prefix_normalization_{1+nnormalization};
for ncovariate=0:n_covariate-1;
prefix_covariate = prefix_covariate_{1+ncovariate};
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
prefix_correction_ = prefix_correction___{1+nX,1+nnormalization,1+ncovariate};
%%%%;
fname_sub_cm__ = sprintf('%s/tlcc_preproj_%s_%s_%s_17.mat',dir_mat,str_X,prefix_normalization,prefix_covariate);
if (verbose); disp(sprintf(' %% collecting fname_sub_cm__: %s',fname_sub_cm__)); end;
if ( exist(fname_sub_cm__,'file'));
if (verbose); disp(sprintf(' %% %s found, not skipping',fname_sub_cm__)); end;
tmp_ = load(fname_sub_cm__);
assert(nX==tmp_.nX);
assert(nnormalization==tmp_.nnormalization);
assert(ncovariate==tmp_.ncovariate);
flag_stage0___{1+nX,1+nnormalization,1+ncovariate} = tmp_.flag_stage0_sub_cm__;
flag_stage1___{1+nX,1+nnormalization,1+ncovariate} = tmp_.flag_stage1_sub_cm__;
n_cluster_stage0___{1+nX,1+nnormalization,1+ncovariate} = tmp_.n_cluster_stage0_sub_cm__;
lpv_A_vs_B___{1+nX,1+nnormalization,1+ncovariate} = tmp_.lpv_A_vs_B_sub_cm__;
lP0_A_vs_B___{1+nX,1+nnormalization,1+ncovariate} = tmp_.lP0_A_vs_B_sub_cm__;
n_X_GENE___{1+nX,1+nnormalization,1+ncovariate} = tmp_.n_X_GENE_sub_cm__;
markergene_A___{1+nX,1+nnormalization,1+ncovariate} = tmp_.markergene_A_sub_cm__;
n_Z_index___{1+nX,1+nnormalization,1+ncovariate} = tmp_.n_Z_index_sub_cm__;
lpv_A_each_vs_B___{1+nX,1+nnormalization,1+ncovariate} = tmp_.lpv_A_each_vs_B_sub_cm__;
nlpv_A_each_vs_B_each___{1+nX,1+nnormalization,1+ncovariate} = tmp_.nlpv_A_each_vs_B_each_sub_cm__;
clear tmp_;
end;%if ( exist(fname_sub_cm__,'file'));
if (~exist(fname_sub_cm__,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',fname_sub_cm__)); end;
end;%if (~exist(fname_sub_cm__,'file'));
%%%%;
end;%for ncovariate=0:n_covariate-1;
end;%for nnormalization=0:n_normalization-1;
end;%for nX=0:n_X-1;
%%%%%%%%;

%%%%%%%%;
% pull out a few methods. ;
%%%%%%%%;
index_method_ = [ ...
 efind(strcmp(prefix_method_,'umap00_default')) ...
,efind(strcmp(prefix_method_,'louvain00_default')) ...
,efind(strcmp(prefix_method_,'dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3')) ...
,efind(strcmp(prefix_method_,'hnbtZRgumb_g010_p050_nml3')) ...
,efind(strcmp(prefix_method_,'preproj_umap00_default')) ...
,efind(strcmp(prefix_method_,'preproj_louvain00_default')) ...
,efind(strcmp(prefix_method_,'preproj_dexcluster_nonbinary_trace_ZRimax_g000_p050_nml3')) ...
,efind(strcmp(prefix_method_,'preproj_hnbtZRgumb_g000_p050_nml3')) ...
,efind(strcmp(prefix_method_,'preproj_hnbr0tZRgumb_g000_p050_nml3')) ...
];

%%%%%%%%;
% Plot size of AIBS clusters. ;
%%%%%%%%;
n_u_label_A_ = zeros(n_label_A,1);
for nlabel_A=0:n_label_A-1;
tmp_index_ = efind(strcmp(str_label_A_{1},u_label_A_{1+nlabel_A}));
n_u_label_A_(1+nlabel_A) = numel(tmp_index_);
end;%for nlabel_A=0:n_label_A-1;
fname_fig = sprintf('%s/n_u_label_A_',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
bar(n_u_label_A_(1+index_label_A_),'k');
set(gca,'XTick',1:n_label_A,'XtickLabel',u_label_A_(1+index_label_A_));
xlabel('AIBS cluster label');
ylabel('number of nuclei');
%%%%%%%%;
figbig;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));

%{
%%%%%%%%;
% Plot rank_estimate_sample for covariates. ;
%%%%%%%%;
A__ = [ones(n_u,1) , mean_center_0(G_gaus_)];
[rank_estimate_sample,s_B_nlp_,svd_sample__,eig_tw__,eig_B_,kta_opt__,h_x__,h_eig__,h_opt__] = rank_estimate_sample_1(A__,0.01,1024);
figure(1); set(gca,'Ytick',1:rank_estimate_sample,'YTickLabel',1:rank_estimate_sample); set(gcf,'Position',1+[0,0,512*1.25,1024]);
print('-depsc','~/dir_bcc/dir_jamison/dir_presentations/rank_estimate_sample_1_G_gaus_n1024.eps');
print('-djpeg','~/dir_bcc/dir_jamison/dir_presentations/rank_estimate_sample_1_G_gaus_n1024.jpg');
%%%%%%%%;
A__ = [ones(n_u,1) , mean_center_0(C_rank_)];
[rank_estimate_sample,s_B_nlp_,svd_sample__,eig_tw__,eig_B_,kta_opt__,h_x__,h_eig__,h_opt__] = rank_estimate_sample_1(A__,0.01,1024);
figure(1); set(gca,'Ytick',1:rank_estimate_sample,'YTickLabel',1:rank_estimate_sample); set(gcf,'Position',1+[0,0,512*1.25,1024]);
print('-depsc','~/dir_bcc/dir_jamison/dir_presentations/rank_estimate_sample_1_C_rank_n1024.eps');
print('-djpeg','~/dir_bcc/dir_jamison/dir_presentations/rank_estimate_sample_1_C_rank_n1024.jpg');
%}

fname_tsv = sprintf('%s/C_VariableName_.tsv',dir_mat);
fp=fopen(fname_tsv); C_VariableName_ = textscan(fp,'%s'); C_VariableName_ = C_VariableName_{1}; fclose(fp);
n_CCOV = numel(C_VariableName_);
fname_tsv = sprintf('%s/G_VariableName_.tsv',dir_mat);
fp=fopen(fname_tsv); G_VariableName_ = textscan(fp,'%s'); G_VariableName_ = G_VariableName_{1}; fclose(fp);
n_GCOV = numel(G_VariableName_);
%%%%%%%%;
% Plot original covariate clusters. ;
%%%%%%%%;
for ncovariate=0:n_covariate-1;
prefix_covariate = prefix_covariate_{1+ncovariate};
tmp_fname_tsv = sprintf('%s/n_QCluster_%s_original.tsv',dir_mat,prefix_covariate);
n_QCluster_original = textread(tmp_fname_tsv); n_QCluster_original = n_QCluster_original(1);
if (verbose); disp(sprintf('%s --> %d',prefix_covariate,n_QCluster_original)); end;
if (prefix_covariate=='C'); n_HCOV = n_CCOV; H_VariableName_ = C_VariableName_; end;
if (prefix_covariate=='G'); n_HCOV = n_GCOV; H_VariableName_ = G_VariableName_; end;
str_tmp = sprintf('%s',prefix_covariate);
fname_fig = sprintf('%s/test_loader_cluster_collect_17p_%s_covariate_tree',dir_jpg,str_tmp);
if (flag_replot | ~exist(sprintf('%s_FIGA.jpg',fname_fig)));
%%%%%%%%;
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 2;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('%s_%s',prefix_covariate,str_xfix);
dir_out = []; E_array_base_ = transpose(zeros(n_u,n_HCOV)); E_array_r_ij_ = []; E_array_c_ij_ = [];
%try;
%%%%%%%%;
[ ...
 ZRimax_output_label_ ...
,ZRimax_lpFmax_label_ ...
,ZRimax_lpnext_label_ ...
] = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_4( ...
 dir_cluster ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,[] ...
,flag_force_create_mat ...
);
label_B_ = label_str_to_enum_0(ZRimax_output_label_); 
%%%%%%%%;
figure(1);
label_plot_recursive_1(ZRimax_output_label_,ZRimax_lpFmax_label_,ZRimax_lpnext_label_,H_VariableName_,[0,27]);
title(str_tmp,'Interpreter','none');
figbig;
print('-depsc',sprintf('%s_FIGA.eps',fname_fig));
print('-djpeg',sprintf('%s_FIGA.jpg',fname_fig));
close(gcf);
%%%%%%%%;
%catch; disp(sprintf(' %% could not process %s',prefix_base)); ;end;%try;
%%%%%%%%;
end;%if (flag_replot | ~exist(sprintf('%s_FIGA.jpg',fname_fig)));
end;%for ncovariate=0:n_covariate-1;
%%%%%%%%;

fname_tsv = sprintf('%s/C_VariableName_.tsv',dir_mat);
fp=fopen(fname_tsv); C_VariableName_ = textscan(fp,'%s'); C_VariableName_ = C_VariableName_{1}; fclose(fp);
n_CCOV = numel(C_VariableName_);
fname_tsv = sprintf('%s/G_VariableName_.tsv',dir_mat);
fp=fopen(fname_tsv); G_VariableName_ = textscan(fp,'%s'); G_VariableName_ = G_VariableName_{1}; fclose(fp);
n_GCOV = numel(G_VariableName_);
%%%%%%%%;
% Plot mahalanobis covariate clusters. ;
%%%%%%%%;
for nX=0:n_X-1;
str_X = str_X_{1+nX};
for nnormalization=0:n_normalization-1;
prefix_normalization = prefix_normalization_{1+nnormalization};
for ncovariate=0:n_covariate-1;
prefix_covariate = prefix_covariate_{1+ncovariate};
str_tmp = sprintf('U%s_%s_%s',str_X,prefix_normalization,prefix_covariate);
fname_fig = sprintf('%s/test_loader_cluster_collect_17p_%s_covariate_tree',dir_jpg,str_tmp);
if (flag_replot | ~exist(sprintf('%s_FIGA.jpg',fname_fig)));
tmp_fname_tsv = sprintf('%s/n_QCluster_U%s_%s_%s_mahalanobis.tsv',dir_mat,str_X,prefix_normalization,prefix_covariate);
n_QCluster_mahalanobis = textread(tmp_fname_tsv); n_QCluster_mahalanobis = n_QCluster_mahalanobis(1);
if (verbose); disp(sprintf('%s_%s_%s --> %d',str_X,prefix_normalization,prefix_covariate,n_QCluster_mahalanobis)); end;
if (prefix_covariate=='C'); n_HCOV = n_CCOV; H_VariableName_ = C_VariableName_; end;
if (prefix_covariate=='G'); n_HCOV = n_GCOV; H_VariableName_ = G_VariableName_; end;
%%%%%%%%;
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 2;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('U%s_%s_%s_%s',str_X,prefix_normalization,prefix_covariate,str_xfix);
dir_out = []; E_array_base_ = transpose(zeros(n_u,n_HCOV)); E_array_r_ij_ = []; E_array_c_ij_ = [];
%try;
%%%%%%%%;
[ ...
 ZRimax_output_label_ ...
,ZRimax_lpFmax_label_ ...
,ZRimax_lpnext_label_ ...
] = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_4( ...
 dir_cluster ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,[] ...
,flag_force_create_mat ...
);
label_B_ = label_str_to_enum_0(ZRimax_output_label_); 
%%%%%%%%;
figure(1);
label_plot_recursive_1(ZRimax_output_label_,ZRimax_lpFmax_label_,ZRimax_lpnext_label_,H_VariableName_,[0,27]);
title(str_tmp,'Interpreter','none');
figbig;
print('-depsc',sprintf('%s_FIGA.eps',fname_fig));
print('-djpeg',sprintf('%s_FIGA.jpg',fname_fig));
close(gcf);
%%%%%%%%;
%catch; disp(sprintf(' %% could not process %s',prefix_base)); ;end;%try;
%%%%%%%%;
end;%if (flag_replot | ~exist(sprintf('%s_FIGA.jpg',fname_fig)));
end;%for ncovariate=0:n_covariate-1;
end;%for nnormalization=0:n_normalization-1;
end;%for nX=0:n_X-1;

%%%%%%%%;
% Collect markergene_A__. ;
%%%%%%%%;
n_skip_sml = 1;
n_skip_med = 2;
n_skip_big = 4;
n_cgnx = 0;
for nX=0:n_X-1;
for nnormalization=0:n_normalization-1;
for ncovariate=0:n_covariate-1;
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
for ncorrection=0:n_correction-1;
n_cgnx = n_cgnx+1;
end;%for ncorrection=0:n_correction-1;
if (ncovariate< n_covariate-1); n_cgnx = n_cgnx+n_skip_sml; end;
end;%for ncovariate=0:n_covariate-1;
if (nnormalization< n_normalization-1); n_cgnx = n_cgnx+n_skip_med; end;
end;%for nnormalization=0:n_normalization-1;
if (nX< n_X-1); n_cgnx = n_cgnx+n_skip_big; end;
end;%for nX=0:n_X-1;
if (verbose); disp(sprintf(' %% n_cgnx %d',n_cgnx)); end;
%%%%%%%%;

prefix_cgnx_ = cell(n_cgnx,1);
ticklabel_cgnx_ = cell(n_cgnx,1);
ncgnx=0;
for nX=0:n_X-1;
str_X = str_X_{1+nX};
for nnormalization=0:n_normalization-1;
prefix_normalization = prefix_normalization_{1+nnormalization};
for ncovariate=0:n_covariate-1;
prefix_covariate = prefix_covariate_{1+ncovariate};
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
prefix_correction_ = prefix_correction___{1+nX,1+nnormalization,1+ncovariate};
for ncorrection=0:n_correction-1;
prefix_correction = prefix_correction_{1+ncorrection};
if (numel(prefix_correction)==0); prefix_cgnx = sprintf('%s_%s',str_X,prefix_normalization); end;
if (numel(prefix_correction)> 0); prefix_cgnx = sprintf('%s_%s_%s',str_X,prefix_normalization,prefix_correction); end;
prefix_cgnx_{1+ncgnx} = prefix_cgnx; 
ticklabel_cgnx = prefix_cgnx; ticklabel_cgnx(strfind(prefix_cgnx,'_')) = ' ';
rank_estimate = rank_estimate___{1+nX,1+nnormalization,1+ncovariate}(1+ncorrection);
if ((ncorrection> 0) & (rank_estimate<=1)); disp(sprintf(' %% dir_%s_cluster ; <-- %d',prefix_cgnx,rank_estimate)); end;
if ((ncorrection> 0) & (rank_estimate> 0)); ticklabel_cgnx = sprintf('%s (r%d)',ticklabel_cgnx,rank_estimate); end;
if ((ncorrection> 0) & (rank_estimate<=0)); ticklabel_cgnx = sprintf('%s (r??)',ticklabel_cgnx); end;
if ((ncorrection==0) & (rank_estimate> 0)); ticklabel_cgnx = sprintf('%s (r??)',ticklabel_cgnx); end;
if ((ncorrection==0) & (rank_estimate<=0)); ticklabel_cgnx = sprintf('%s (ori)',ticklabel_cgnx); end;
ticklabel_cgnx_{1+ncgnx} = ticklabel_cgnx;
ncgnx = ncgnx+1;
end;%for ncorrection=0:n_correction-1;
if (ncovariate< n_covariate-1);
for nl=0:n_skip_sml-1;
prefix_cgnx_{1+ncgnx} = '';
ticklabel_cgnx_{1+ncgnx} = '';
ncgnx = ncgnx+1;
end;%for nl=0:n_skip_sml-1;
end;%if (ncovariate< n_covariate-1);
end;%for ncovariate=0:n_covariate-1;
if (nnormalization< n_normalization-1);
for nl=0:n_skip_med-1;
prefix_cgnx_{1+ncgnx} = '';
ticklabel_cgnx_{1+ncgnx} = '';
ncgnx = ncgnx+1;
end;%for nl=0:n_skip_med-1;
end;%if (nnormalization< n_normalization-1);
end;%for nnormalization=0:n_normalization-1;
if (nX< n_X-1);
for nl=0:n_skip_big-1;
prefix_cgnx_{1+ncgnx} = '';
ticklabel_cgnx_{1+ncgnx} = '';
ncgnx = ncgnx+1;
end;%for nl=0:n_skip_big-1;
end;%if (nX< n_X-1);
end;%for nX=0:n_X-1;
assert(ncgnx==n_cgnx);

%%%%%%%%;
markergene_correlation_A_cgnx__ = -Inf*ones(n_label_A,n_cgnx);
markergene_top01_A_cgnx__ = -Inf*ones(n_label_A,n_cgnx);
ncgnx=0;
for nX=0:n_X-1;
for nnormalization=0:n_normalization-1;
for ncovariate=0:n_covariate-1;
n_X_GENE = n_X_GENE___{1+nX,1+nnormalization,1+ncovariate}(1+0);
n_top01 = ceil(n_X_GENE*0.01);
tmp_markergene_A_auz_0__ = markergene_A___{1+nX,1+nnormalization,1+ncovariate}{1+0};
[~,tmp_ij_0__] = sort(abs(tmp_markergene_A_auz_0__),1,'descend');
tmp_ij_0_top01__ = tmp_ij_0__(1:n_top01,:);
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
for ncorrection=0:n_correction-1;
tmp_markergene_A_auz_1__ = markergene_A___{1+nX,1+nnormalization,1+ncovariate}{1+1*ncorrection};
if (numel(tmp_markergene_A_auz_0__)==numel(tmp_markergene_A_auz_1__));
tmp_C__ = corr(tmp_markergene_A_auz_0__,tmp_markergene_A_auz_1__);
markergene_correlation_A_cgnx__(:,1+ncgnx) = diag(tmp_C__);
[~,tmp_ij_1__] = sort(abs(tmp_markergene_A_auz_1__),1,'descend');
tmp_ij_1_top01__ = tmp_ij_1__(1:n_top01,:);
for nlabel_A=0:n_label_A-1;
markergene_top01_A_cgnx__(1+nlabel_A,1+ncgnx) = numel(intersect(tmp_ij_0_top01__(:,1+nlabel_A),tmp_ij_1_top01__(:,1+nlabel_A)))/n_top01;
end;%for nlabel_A=0:n_label_A-1;
end;%if (numel(tmp_markergene_A_auz_0__)==numel(tmp_markergene_A_auz_1__));
ncgnx = ncgnx+1;
end;%for ncorrection=0:n_correction-1;
if (ncovariate< n_covariate-1); ncgnx = ncgnx+n_skip_sml; end;
end;%for ncovariate=0:n_covariate-1;
if (nnormalization< n_normalization-1); ncgnx = ncgnx+n_skip_med; end;
end;%for nnormalization=0:n_normalization-1;
if (nX< n_X-1); ncgnx = ncgnx+n_skip_big; end;
end;%for nX=0:n_X-1;
assert(ncgnx==n_cgnx);
%%%%%%%%;

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_markergene_A__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpv(64,0.85,0.15,0.85); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
subplot(2,1,1);
tmp_lim_ = [0,1]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = markergene_correlation_A_cgnx__(1+index_label_A_,:);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_));
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title('correlation');
set(gca,'FontSize',6);
axis image;
subplot(2,1,2);
tmp_lim_ = [0,1]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = markergene_top01_A_cgnx__(1+index_label_A_,:);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_));
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title('top01');
set(gca,'FontSize',6);
axis image;
figbig; set(gcf,'Color',[1,1,1]);
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% collect lpv_A_each_vs_B__. ; 
% ordering: ;
% Aa: [all , each(sorted) ]. ;
% cgnx: ncorrection, nnormalization, nX ;
% zm: % nZ_index, nmethod ;
%%%%%%%%;
n_Z_index_max_ = zeros(n_method,1);
for nX=0:n_X-1;
for nnormalization=0:n_normalization-1;
for ncovariate=0:n_covariate-1;
n_Z_index_max_ = max(n_Z_index_max_,transpose(max(n_Z_index___{1+nX,1+nnormalization,1+ncovariate},[],1)));
end;%for ncovariate=0:n_covariate-1;
end;%for nnormalization=0:n_normalization-1;
end;%for nX=0:n_X-1;
n_zm = sum(n_Z_index_max_);
prefix_zm_ = cell(n_zm,1);
%%%%%%%%;
nzm=0;
for nmethod=0:n_method-1;
n_Z_index = n_Z_index_max_(1+nmethod);
for nZ_index=0:n_Z_index-1;
if ( isempty(strfind(prefix_method_{1+nmethod},'dexcluster')));
prefix_zm = sprintf('%s (r%d)',prefix_method_{1+nmethod},1+nZ_index);
prefix_zm(strfind(prefix_zm,'_')) = ' ';
prefix_zm_{1+nzm} = prefix_zm;
end;%if ( isempty(strfind(prefix_method_{1+nmethod},'dexcluster')));
%%%%;
if (~isempty(strfind(prefix_method_{1+nmethod},'dexcluster'))); 
if (~isempty(strfind(prefix_method_{1+nmethod},'preproj'))); 
if (nZ_index==0); prefix_zm_{1+nzm} = sprintf('preproj loop (pre)'); end;
if (nZ_index==1); prefix_zm_{1+nzm} = sprintf('preproj loop (pos)'); end;
else;
if (nZ_index==0); prefix_zm_{1+nzm} = sprintf('loop (pre)'); end;
if (nZ_index==1); prefix_zm_{1+nzm} = sprintf('loop (pos)'); end;
end;%if (~isempty(strfind(prefix_method_{1+nmethod},'preproj'))); 
end;%if (~isempty(strfind(prefix_method_{1+nmethod},'dexcluster'))); 
%%%%;
if (~isempty(strfind(prefix_method_{1+nmethod},'hnbr0t'))); 
if (~isempty(strfind(prefix_method_{1+nmethod},'preproj'))); 
if (nZ_index==0); prefix_zm_{1+nzm} = sprintf('preproj hnbr0t (pre)'); end;
if (nZ_index==1); prefix_zm_{1+nzm} = sprintf('preproj hnbr0t (pos)'); end;
else;
if (nZ_index==0); prefix_zm_{1+nzm} = sprintf('hnbr0t (pre)'); end;
if (nZ_index==1); prefix_zm_{1+nzm} = sprintf('hnbr0t (pos)'); end;
end;%if (~isempty(strfind(prefix_method_{1+nmethod},'preproj'))); 
end;%if (~isempty(strfind(prefix_method_{1+nmethod},'hnbr0t'))); 
%%%%;
if (~isempty(strfind(prefix_method_{1+nmethod},'hnbtZ'))); 
if (~isempty(strfind(prefix_method_{1+nmethod},'preproj'))); 
if (nZ_index==0); prefix_zm_{1+nzm} = sprintf('preproj hnbt (pre)'); end;
if (nZ_index==1); prefix_zm_{1+nzm} = sprintf('preproj hnbt (pos)'); end;
else;
if (nZ_index==0); prefix_zm_{1+nzm} = sprintf('hnbt (pre)'); end;
if (nZ_index==1); prefix_zm_{1+nzm} = sprintf('hnbt (pos)'); end;
end;%if (~isempty(strfind(prefix_method_{1+nmethod},'preproj'))); 
end;%if (~isempty(strfind(prefix_method_{1+nmethod},'hnbt'))); 
%%%%;
nzm = nzm+1;
end;%for nZ_index=0:n_Z_index-1;
end;%for nmethod=0:n_method-1;
assert(nzm==n_zm);
%%%%%%%%;

%%%%%%%%;
% pull out a few zm. ;
%%%%%%%%;
index_zm_ = [ ...
 efind(strcmp(prefix_zm_,'umap00 default (r1)')) ...
,efind(strcmp(prefix_zm_,'louvain00 default (r1)')) ...
,efind(strcmp(prefix_zm_,'loop (pre)')) ...
,efind(strcmp(prefix_zm_,'loop (pos)')) ...
,efind(strcmp(prefix_zm_,'hnbt (pre)')) ...
,efind(strcmp(prefix_zm_,'hnbt (pos)')) ...
,efind(strcmp(prefix_zm_,'preproj umap00 default (r1)')) ...
,efind(strcmp(prefix_zm_,'preproj louvain00 default (r1)')) ...
,efind(strcmp(prefix_zm_,'preproj loop (pre)')) ...
,efind(strcmp(prefix_zm_,'preproj loop (pos)')) ...
,efind(strcmp(prefix_zm_,'preproj hnbt (pre)')) ...
,efind(strcmp(prefix_zm_,'preproj hnbt (pos)')) ...
,efind(strcmp(prefix_zm_,'preproj hnbr0t (pre)')) ...
,efind(strcmp(prefix_zm_,'preproj hnbr0t (pos)')) ...
];

n_cluster_cgnx_zm__ = -Inf*ones(n_cgnx,n_zm);
nlpv_A_vs_B_Aa_cgnx_zm___ = -Inf*ones(n_label_A+1,n_cgnx,n_zm);
ncgnx=0;
for nX=0:n_X-1;
for nnormalization=0:n_normalization-1;
for ncovariate=0:n_covariate-1;
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
for ncorrection=0:n_correction-1;
nzm=0;
for nmethod=0:n_method-1;
n_Z_index = n_Z_index_max_(1+nmethod);
for nZ_index=0:n_Z_index-1;
n_cluster_stage0_ = n_cluster_stage0___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection,1+nmethod};
if (~isempty(n_cluster_stage0_));
n_cluster_cgnx_zm__(1+ncgnx,1+nzm) = n_cluster_stage0_(1+nZ_index);
if (strcmp(prefix_zm_(1+nzm),'loop (pos)')); disp(sprintf(' %% %s (%s) <-- n_cluster %d',prefix_cgnx_{1+ncgnx},prefix_zm_{1+nzm},n_cluster_cgnx_zm__(1+ncgnx,1+nzm))); end;
end;%if (~isempty(n_cluster_stage0_));
tmp_lpv__ = lpv_A_each_vs_B___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection,1+nmethod};
if (~isempty(tmp_lpv__)); nlpv_A_vs_B_Aa_cgnx_zm___(1:n_label_A,1+ncgnx,1+nzm) = -tmp_lpv__(1+index_label_A_,1+nZ_index); end;%if (~isempty(tmp_lpv__));
tmp_lpv_ = lpv_A_vs_B___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection,1+nmethod};
if (~isempty(tmp_lpv_)); nlpv_A_vs_B_Aa_cgnx_zm___(1+n_label_A,1+ncgnx,1+nzm) = -tmp_lpv_(1+nZ_index); end;%if (~isempty(tmp_lpv_));
nzm = nzm+1;
end;%for nZ_index=0:n_Z_index-1;
end;%for nmethod=0:n_method-1;
ncgnx = ncgnx+1;
end;%for ncorrection=0:n_correction-1;
if (ncovariate< n_covariate-1); ncgnx = ncgnx+n_skip_sml; end;
end;%for ncovariate=0:n_covariate-1;
if (nnormalization< n_normalization-1); ncgnx = ncgnx+n_skip_med; end;
end;%for nnormalization=0:n_normalization-1;
if (nX< n_X-1); ncgnx = ncgnx+n_skip_big; end;
end;%for nX=0:n_X-1;
assert(nzm==n_zm);
assert(ncgnx==n_cgnx);
%%%%%%%%;
ticklabel_Aa_ = u_label_A_(1+index_label_A_);
ticklabel_Aa_{n_label_A+1} = 'all';

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_lpv_A_each_vs_B__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
%c_ = colormap('hot'); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
np=0;n_prows=4;n_pcols=ceil(numel(index_zm_)/n_prows);
for nizm=0:numel(index_zm_)-1;
nzm = index_zm_(1+nizm);
subplot(n_prows,n_pcols,1+np); cla;
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+nzm);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'YTick',[],'YTickLabel',[]); %<-- labels off. ;
set(gca,'XTick',[],'XTickLabel',[]); %<-- labels off. ;
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('nizm%d(nzm%d):%s',nizm,nzm,prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
np = np+1;
end;%for nizm=0:numel(index_zm_)-1;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_lpv_A_each_vs_B_max__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = max(nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+[0:n_zm-1]),[],3) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('max -log(p) (bonf)'),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% plot lpv_A_each_vs_B_uma__ for a few methods. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_lpv_A_each_vs_B_umap__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
nzm = efind(strcmp(prefix_zm_,'umap00 default (r1)'));
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+nzm) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s -log(p) (bonf)',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_lpv_A_each_vs_B_dex__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
nzm = efind(strcmp(prefix_zm_,'loop (pos)'));
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+nzm) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s -log(p) (bonf)',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_lpv_A_each_vs_B_hnb__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
nzm = efind(strcmp(prefix_zm_,'hnbt (pos)'));
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+nzm) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s -log(p) (bonf)',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_lpv_A_each_vs_B_preproj_umap__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
nzm = efind(strcmp(prefix_zm_,'preproj umap00 default (r1)'));
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+nzm) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s -log(p) (bonf)',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_lpv_A_each_vs_B_preproj_dex__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
nzm = efind(strcmp(prefix_zm_,'preproj loop (pos)'));
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+nzm) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s -log(p) (bonf)',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_lpv_A_each_vs_B_preproj_hnb__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
nzm = efind(strcmp(prefix_zm_,'preproj hnbt (pos)'));
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+nzm) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s -log(p) (bonf)',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17p_lpv_A_each_vs_B_preproj_hnbr0__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
nzm = efind(strcmp(prefix_zm_,'preproj hnbr0t (pos)'));
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+nzm) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s -log(p) (bonf)',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

symbol_zm_ = cell(n_zm,1);
msize_zm_ = zeros(n_zm,1);
nzm=0;
symbol_zm_{1+nzm} = 'go'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'gx'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'g^'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'gs'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'g*'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'gh'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'co'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'cx'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'bo'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'bx'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'kx'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'ko'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'k^'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'rx'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'mo'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'ms'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'mo'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'ms'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'go'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'co'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'cx'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'bo'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'bx'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'kx'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'ko'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'k^'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'rx'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'mo'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'ms'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'mo'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'ms'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'mo'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
symbol_zm_{1+nzm} = 'ms'; msize_zm_(1+nzm) = 15; nzm=nzm+1;
assert(nzm==n_zm);

for tmp_infix_ = {'logn','liXf'};
fname_fig = sprintf('%s/test_loader_cluster_collect_17p_n_cluster_%s__',dir_jpg,tmp_infix_{1});
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
subplot(2,3,1); cla;
hold on;
for nzm=0:n_zm-1;
plot(1:n_cgnx,reshape(nlpv_A_vs_B_Aa_cgnx_zm___(n_label_A+1,:,1+nzm),[n_cgnx,1]),symbol_zm_{1+nzm},'MarkerSize',msize_zm_(1+nzm));
end;%for nzm=0:n_zm-1;
hold off;
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
xlabel('data normalization+correction');
ylabel('log(p)');
set(gca,'TickLength',[0,0]);
%legend(prefix_zm_);
title(sprintf('log(p)'));
set(gca,'FontSize',6);
subplot(2,3,4); cla;
hold on;
for nzm=0:n_zm-1;
plot(1:n_cgnx,reshape(n_cluster_cgnx_zm__(:,1+nzm),[n_cgnx,1]),symbol_zm_{1+nzm},'MarkerSize',msize_zm_(1+nzm));
end;%for nzm=0:n_zm-1;
hold off;
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
xlabel('data normalization+correction');
ylabel('# clusters');
set(gca,'TickLength',[0,0]);
%legend(prefix_zm_);
title(sprintf('# clusters'));
set(gca,'FontSize',6);
subplot(2,3,[2,3,5,6]); cla;
hold on;
for nzm=0:n_zm-1;
tmp_index_ = efind(~cellfun(@isempty,strfind(prefix_cgnx_,tmp_infix_{1})));
tmp_n_cgnx = numel(tmp_index_);
plot( ...
 reshape(n_cluster_cgnx_zm__(1+tmp_index_,1+nzm),[tmp_n_cgnx,1]) ...
,reshape(nlpv_A_vs_B_Aa_cgnx_zm___(n_label_A+1,1+tmp_index_,1+nzm),[tmp_n_cgnx,1]) ...
,symbol_zm_{1+nzm},'MarkerSize',msize_zm_(1+nzm) ...
);
end;%for nzm=0:n_zm-1;
hold off;
xlabel('# clusters');
ylabel('log(p)');
%legend(prefix_zm_);
title(sprintf('scatterplot (liXf)'));
set(gca,'FontSize',6);
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for tmp_infix_ = {'logn','liXf'};

for tmp_infix_ = {'logn','liXf'};
fname_fig = sprintf('%s/test_loader_cluster_collect_17p_E_%s_I_%s_n_cluster__',dir_jpg,tmp_infix_{1},tmp_infix_{1});
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%%%%%;
figure(1);clf;
subplot(1,2,1); cla;
hold on;
nzm_use_ = [ ...
      efind(strcmp(prefix_zm_,'spectral isosplit5 (r1)')) ...
     ,efind(strcmp(prefix_zm_,'spectral isosplit5 (r2)')) ...
     ,efind(strcmp(prefix_zm_,'spectral isosplit5 (r3)')) ...
     ,efind(strcmp(prefix_zm_,'spectral isosplit5 (r4)')) ...
     ,efind(strcmp(prefix_zm_,'spectral isosplit5 (r5)')) ...
     ,efind(strcmp(prefix_zm_,'spectral isosplit5 (r6)')) ...
     ,efind(strcmp(prefix_zm_,'tsne00 isosplit5 (r2)')) ...
     ,efind(strcmp(prefix_zm_,'tsne50 isosplit5 (r2)')) ...
     ,efind(strcmp(prefix_zm_,'umap00 default (r1)')) ...
     ,efind(strcmp(prefix_zm_,'umap00 isosplit5 (r1)')) ...
     ,efind(strcmp(prefix_zm_,'umap00 hdbscan (r1)')) ...
     ,efind(strcmp(prefix_zm_,'louvain00 default (r1)')) ...
     ,efind(strcmp(prefix_zm_,'loop (pre)')) ...
     ,efind(strcmp(prefix_zm_,'loop (pos)')) ...
	     ];
legend_use_ = prefix_zm_(1+nzm_use_);
%%%%%%%%;
subplot(1,2,1); cla;
hold on;
for nzm=nzm_use_;
tmp_index_ = min(efind(strcmp(prefix_cgnx_,sprintf('E_%s',tmp_infix_{1}))));
plot( ...
 reshape(n_cluster_cgnx_zm__(1+tmp_index_,1+nzm),[1,1]) ...
,reshape(nlpv_A_vs_B_Aa_cgnx_zm___(n_label_A+1,1+tmp_index_,1+nzm),[1,1]) ...
,symbol_zm_{1+nzm},'MarkerSize',2*msize_zm_(1+nzm) ...
,'LineWidth',4 ...
);
end;%for nzm=0:n_zm-1;
hold off;
xlim([0,55]);ylim([0,1500]);
%legend(legend_use_,'location','SouthEast');
xlabel('# clusters');
ylabel('-log(p)');
title(sprintf('scatterplot (E_%s)',tmp_infix_{1}),'Interpreter','none');
%set(gca,'FontSize',8);
%%%%%%%%;
subplot(1,2,2); cla;
hold on;
for nzm=nzm_use_;
tmp_index_ = min(efind(strcmp(prefix_cgnx_,sprintf('I_%s',tmp_infix_{1}))));
plot( ...
 reshape(n_cluster_cgnx_zm__(1+tmp_index_,1+nzm),[1,1]) ...
,reshape(nlpv_A_vs_B_Aa_cgnx_zm___(n_label_A+1,1+tmp_index_,1+nzm),[1,1]) ...
,symbol_zm_{1+nzm},'MarkerSize',2*msize_zm_(1+nzm) ...
,'LineWidth',4 ...
);
end;%for nzm=0:n_zm-1;
hold off;
xlim([0,55]);ylim([0,1500]);
legend(legend_use_,'location','NorthWest');
xlabel('# clusters');
ylabel('-log(p)');
title(sprintf('scatterplot (I_%s)',tmp_infix_{1}),'Interpreter','none');
%set(gca,'FontSize',8);
%%%%%%%%;
figbig;
%%%%%%%%;
if (verbose); disp(sprintf(' %% writing %s',fname_fig)); end;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for tmp_infix_ = {'logn','liXf'};

%%%%%%%%;
% Now delve into E_logn and E_liXf
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% plot A_each_vs_B_each for E_liXf_ and E_liXf_G_call_. ;
%%%%%%%%;
n_E_GENE = n_X_GENE___{1,1,1}(1); n_I_GENE = n_X_GENE___{2,1,1}(1);
for tmp_prefix_cgnx_ = {'E_logn','E_liXf','I_logn','I_liXf'};
%%%%%%%%;
prefix_cgnx = tmp_prefix_cgnx_{1};
if (~isempty(strfind(prefix_cgnx,'E_')));
if (~isempty(strfind(prefix_cgnx,'logn'))); X_ = E_logn_; end;
if (~isempty(strfind(prefix_cgnx,'liXf'))); X_ = E_liXf_; end;
n_X_GENE = 1*n_E_GENE + 0*n_I_GENE; rank_estimate_X = rank_estimate_E;
end;%if (~isempty(strfind(prefix_cgnx,'E_')));
if (~isempty(strfind(prefix_cgnx,'I_'))); 
if (~isempty(strfind(prefix_cgnx,'logn'))); X_ = I_logn_; end;
if (~isempty(strfind(prefix_cgnx,'liXf'))); X_ = I_liXf_; end;
n_X_GENE = 0*n_E_GENE + 1*n_I_GENE; rank_estimate_X = rank_estimate_I;
end;%if (~isempty(strfind(prefix_cgnx,'I_'))); 
ncgnx = min(efind(strcmp(prefix_cgnx_,prefix_cgnx)));
assert(strcmp(prefix_cgnx_{1+ncgnx},prefix_cgnx));
tmp_dir_cluster = sprintf('%s/dir_%s_cluster',dir_cluster,prefix_cgnx);
%%%%%%%%;
% loading local hnbr0tZRgumb, hnbrtZRgumb and hnbtZRgumb results. ;
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for str_local_method_ = {'hnbr0tZRgumb','hnbrtZRgumb','hnbtZRgumb'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
str_local_method = str_local_method_{1};
if (verbose); disp(sprintf(' %% postprocessing for %s',str_local_method)); end;
%%%%;
if (strcmp(str_local_method,'hnbr0tZRgumb'));
tmp_infix_0 = sprintf('preproj_hnbr0tZRgumb_g000_p050_nml3'); 
tmp_infix_1 = sprintf('%s_n03',tmp_infix_0);
str_infix = tmp_infix_1;
fname_louvain_mat = sprintf('%s/preproj_louvain00_default.mat',tmp_dir_cluster);
fname_umap_mat = sprintf('%s/preproj_umap00_default.mat',tmp_dir_cluster);
end;%if (strcmp(str_local_method,'hnbr0tZRgumb'));
%%%%;
if (strcmp(str_local_method,'hnbrtZRgumb'));
tmp_infix_0 = sprintf('preproj_hnbtZRgumb_g000_p050_nml3'); 
tmp_infix_1 = sprintf('%s_n03',tmp_infix_0);
str_infix = tmp_infix_1;
fname_louvain_mat = sprintf('%s/preproj_louvain00_default.mat',tmp_dir_cluster);
fname_umap_mat = sprintf('%s/preproj_umap00_default.mat',tmp_dir_cluster);
end;%if (strcmp(str_local_method,'hnbrtZRgumb'));
%%%%;
if (strcmp(str_local_method,'hnbtZRgumb'));
tmp_infix_0 = sprintf('hnbtZRgumb_g010_p050_nml3'); 
tmp_infix_1 = sprintf('%s_g010_n03',tmp_infix_0);
str_infix = tmp_infix_1;
fname_louvain_mat = sprintf('%s/louvain00_default.mat',tmp_dir_cluster);
fname_umap_mat = sprintf('%s/umap00_default.mat',tmp_dir_cluster);
end;%if (strcmp(str_local_method,'hnbtZRgumb'));
%%%%;
tmp_dir_cluster_sub = sprintf('%s/dir_tmp_%s/dir_tmp_%s',tmp_dir_cluster,tmp_infix_0,tmp_infix_1);
local_fname_output_label = sprintf('%s/output_label__.txt',tmp_dir_cluster_sub);
local_fname_nlpbra_label = sprintf('%s/nlpbra_label__.txt',tmp_dir_cluster_sub);
local_fname_nlpnex_label = sprintf('%s/nlpnex_label__.txt',tmp_dir_cluster_sub);
str_sgtitle = sprintf('%s',tmp_dir_cluster_sub);
%%%%;
parameter = struct('type','parameter');
parameter.flag_force_create_mat = 0;
parameter.flag_force_create_tmp = 0;
parameter.flag_replot = flag_replot;
label_tree_compilation_plot_0( ...
 parameter ...
,tmp_dir_cluster ...
,str_infix ...
,str_sgtitle ...
,local_fname_output_label ...
,local_fname_nlpbra_label ...
,local_fname_nlpnex_label ...
,X_ ...
,rank_estimate_X ...
,label_A_ ...
,fname_louvain_mat ...
,fname_umap_mat ...
);
%%%%;

parameter = struct('type','parameter');
parameter.flag_force_create_mat = 0;
parameter.flag_force_create_tmp = 0;
parameter.flag_replot = flag_replot;
parameter.SVD_discat_alpha_over_beta = 1;
parameter.label_tree_SVD_discat_depth_upb = 1;
label_tree_SVD_discat_wrap_0( ...
 parameter ...
,tmp_dir_cluster ...
,str_infix ...
,str_sgtitle ...
,local_fname_output_label ...
,local_fname_nlpbra_label ...
,local_fname_nlpnex_label ...
,X_ ...
,rank_estimate_X ...
,label_A_ ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for str_local_method_ = {'hnbr0tZRgumb','hnbrtZRgumb','hnbtZRgumb'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for tmp_prefix_cgnx_ = {'E_logn','E_liXf'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

end;%if flag_plot_make;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

disp('returning'); return;

%%%%%%%%;
% quick text dump of output. ;
%%%%%%%%;

for nX=0:n_X-1;
str_X = str_X_{1+nX};
for nnormalization=0:n_normalization-1;
prefix_normalization = prefix_normalization_{1+nnormalization};
for ncovariate=0;%for ncovariate=0:n_covariate-1;
prefix_covariate = prefix_covariate_{1+ncovariate};
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
prefix_correction_ = prefix_correction___{1+nX,1+nnormalization,1+ncovariate};
ncorrection=0;
prefix_correction = prefix_correction_{1+ncorrection};
for nl=0:numel(index_method_)-1;
if (mod(nl,4)==0); disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')); end;
prefix_method = prefix_method_{1+index_method_(1+nl)};
[tmp_lpv,tmp_ij] = max(lpv_A_vs_B___{1+nX,1+nnormalization,1+ncovariate}{1,1+index_method_(1+nl)});
tmp_n_cluster = n_cluster_stage0___{1+nX,1+nnormalization,1+ncovariate}{1,1+index_method_(1+nl)}(tmp_ij);
disp(sprintf(' %% %2s %10s %1s %1s %64s: lpv %16.2f n_cluster %16d',str_X,prefix_normalization,prefix_covariate,prefix_correction,prefix_method,tmp_lpv,tmp_n_cluster));
end;%for nl=0:numel(index_method_)-1;
end;%for ncovariate=0:n_covariate-1;
end;%for nnormalization=0:n_normalization-1;
end;%for nX=0:n_X-1;
