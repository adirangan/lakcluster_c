%function test_loader_cluster_collect_17();
% intended for use with test_loader_38.m ;

clear;

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
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_dex = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
%%%%%%%%;
str_X_ = {'E','I','EI'}; n_X = numel(str_X_);
prefix_normalization_ = {'logn','liXf','liXf_clogr','fill'}; n_normalization = numel(prefix_normalization_);
prefix_covariate_ = {'C','G'}; n_covariate = numel(prefix_covariate_);
prefix_method_ = {'spectral_isosplit5','tsne00_isosplit5','tsne50_isosplit5','umap00_default','umap00_isosplit5','umap00_hdbscan','louvain00_default',str_dex}; n_method = numel(prefix_method_);
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
disp(sprintf(' %% %s n_QCluster_original %d',prefix_covariate,n_QCluster_original));
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
disp(sprintf(' %% %s_%s_%s n_QCluster_mahalanobis %d',str_X,prefix_normalization,prefix_covariate,n_QCluster_mahalanobis));
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
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
prefix_correction_ = prefix_correction___{1+nX,1+nnormalization,1+ncovariate};
for ncorrection=0:n_correction-1;
prefix_correction = prefix_correction_{1+ncorrection};
if 0;
elseif (numel(prefix_correction)==0); tmp_dir_cluster = sprintf('%s/dir_%s_%s_cluster',dir_cluster,str_X,prefix_normalization); 
elseif (numel(prefix_correction)> 0); tmp_dir_cluster = sprintf('%s/dir_%s_%s_%s_cluster',dir_cluster,str_X,prefix_normalization,prefix_correction);
end;%if;
disp(sprintf(' %% dir: %s',tmp_dir_cluster));
tmp_fname_mat = sprintf('%s/markergene_A__.mat',tmp_dir_cluster);
if (~exist(tmp_fname_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',tmp_fname_mat)); end;
end;%if (~exist(tmp_fname_mat,'file'));
if ( exist(tmp_fname_mat,'file'));
tmp_ = load(tmp_fname_mat);
n_X_GENE___{1+nX,1+nnormalization,1+ncovariate}(1+ncorrection) = tmp_.n_X_GENE;
markergene_A___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection} = tmp_.markergene_auz_A__;
clear tmp_;
end;%if ( exist(tmp_fname_mat,'file'));
for nmethod=0:n_method-1;
prefix_method = prefix_method_{1+nmethod};
tmp_fname_mat = sprintf('%s/%s.mat',tmp_dir_cluster,prefix_method);
if (~exist(tmp_fname_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',tmp_fname_mat)); end;
end;%if (~exist(tmp_fname_mat,'file'));
if ( exist(tmp_fname_mat,'file'));
flag_stage0___{1+nX,1+nnormalization,1+ncovariate}(1+ncorrection,1+nmethod) = 1;
tmp_ = load(tmp_fname_mat);
n_Z_index = numel(tmp_.label_B__);
n_cluster_stage0___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection,1+nmethod} = zeros(n_Z_index,1);
lpv_A_vs_B___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection,1+nmethod} = tmp_.lpv_;
lP0_A_vs_B___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection,1+nmethod} = tmp_.lP0_;
for nZ_index=0:n_Z_index-1;
n_cluster_stage0___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection,1+nmethod}(1+nZ_index) = numel(unique(tmp_.label_B__{1+nZ_index}));
end;%for nZ_index=0:n_Z_index-1;
clear tmp_;
end;%if ( exist(tmp_fname_mat,'file'));
tmp_fname_mat = sprintf('%s/%s_markergene_B___.mat',tmp_dir_cluster,prefix_method);
if (~exist(tmp_fname_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',tmp_fname_mat)); end;
end;%if (~exist(tmp_fname_mat,'file'));
if ( exist(tmp_fname_mat,'file'));
flag_stage1___{1+nX,1+nnormalization,1+ncovariate}(1+ncorrection,1+nmethod) = 1;
tmp_ = load(tmp_fname_mat);
n_Z_index___{1+nX,1+nnormalization,1+ncovariate}(1+ncorrection,1+nmethod) = tmp_.n_Z_index;
lpv_A_each_vs_B___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection,1+nmethod} = tmp_.lpv_A_each_vs_B__;
nlpv_A_each_vs_B_each___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection,1+nmethod} = tmp_.nlpv_each___;
clear tmp_;
end;%if ( exist(tmp_fname_mat,'file'));
end;%for nmethod=0:n_method-1;
end;%for ncorrection=0:n_correction-1;
end;%for ncovariate=0:n_covariate-1;
end;%for nnormalization=0:n_normalization-1;
end;%for nX=0:n_X-1;
%%%%%%%%;

save(sprintf('%s/tlcc17.mat',dir_mat),'-v7.3');

end;%if flag_compute;

if flag_plot_make;

%%%%%%%%;
% copy to OptiPlex.; 
%%%%%%%%;
% scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_jamison/dir_mat/tlcc17.mat /home/rangan/dir_bcc/dir_jamison/dir_mat/tlcc17.mat;

load(sprintf('/home/rangan/dir_bcc/dir_jamison/dir_mat/tlcc17.mat'));

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
% Plot size of AIBS clusters. ;
%%%%%%%%;
n_u_label_A_ = zeros(n_label_A,1);
for nlabel_A=0:n_label_A-1;
tmp_index_ = efind(strcmp(str_label_A_{1},u_label_A_{1+nlabel_A}));
n_u_label_A_(1+nlabel_A) = numel(tmp_index_);
end;%for nlabel_A=0:n_label_A-1;
fname_fig = sprintf('%s/n_u_label_A_',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
print('-depsc','~/dir_bcc/dir_jamison/dir_presentations/rank_estiamte_sample_1_G_gaus_n1024.eps');
print('-djpeg','~/dir_bcc/dir_jamison/dir_presentations/rank_estiamte_sample_1_G_gaus_n1024.jpg');
%%%%%%%%;
A__ = [ones(n_u,1) , mean_center_0(C_rank_)];
[rank_estimate_sample,s_B_nlp_,svd_sample__,eig_tw__,eig_B_,kta_opt__,h_x__,h_eig__,h_opt__] = rank_estimate_sample_1(A__,0.01,1024);
figure(1); set(gca,'Ytick',1:rank_estimate_sample,'YTickLabel',1:rank_estimate_sample); set(gcf,'Position',1+[0,0,512*1.25,1024]);
print('-depsc','~/dir_bcc/dir_jamison/dir_presentations/rank_estiamte_sample_1_C_rank_n1024.eps');
print('-djpeg','~/dir_bcc/dir_jamison/dir_presentations/rank_estiamte_sample_1_C_rank_n1024.jpg');
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
disp(sprintf('%s --> %d',prefix_covariate,n_QCluster_original));
if (prefix_covariate=='C'); n_HCOV = n_CCOV; H_VariableName_ = C_VariableName_; end;
if (prefix_covariate=='G'); n_HCOV = n_GCOV; H_VariableName_ = G_VariableName_; end;
str_tmp = sprintf('%s',prefix_covariate);
fname_fig = sprintf('%s/test_loader_cluster_collect_17_%s_covariate_tree',dir_jpg,str_tmp);
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
fname_fig = sprintf('%s/test_loader_cluster_collect_17_%s_covariate_tree',dir_jpg,str_tmp);
if (flag_replot | ~exist(sprintf('%s_FIGA.jpg',fname_fig)));
tmp_fname_tsv = sprintf('%s/n_QCluster_U%s_%s_%s_mahalanobis.tsv',dir_mat,str_X,prefix_normalization,prefix_covariate);
n_QCluster_mahalanobis = textread(tmp_fname_tsv); n_QCluster_mahalanobis = n_QCluster_mahalanobis(1);
disp(sprintf('%s_%s_%s --> %d',str_X,prefix_normalization,prefix_covariate,n_QCluster_mahalanobis));
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
disp(sprintf(' %% n_cgnx %d',n_cgnx));
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

fname_fig = sprintf('%s/test_loader_cluster_collect_17_markergene_A__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
disp(sprintf(' %% writing %s',fname_fig));
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
if (~isempty(strfind(prefix_method_{1+nmethod},'dexcluster'))); 
if (nZ_index==0); prefix_zm_{1+nzm} = sprintf('loop (pre)'); end;
if (nZ_index==1); prefix_zm_{1+nzm} = sprintf('loop (pos)'); end;
end;%if (~isempty(strfind(prefix_method_{1+nmethod},'dexcluster'))); 
nzm = nzm+1;
end;%for nZ_index=0:n_Z_index-1;
end;%for nmethod=0:n_method-1;
assert(nzm==n_zm);
%%%%%%%%;

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

fname_fig = sprintf('%s/test_loader_cluster_collect_17_lpv_A_each_vs_B__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(1);clf;
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
%c_ = colormap('hot'); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
np=0;n_prows=4;n_pcols=4;
for nzm=0:n_zm-1;
subplot(n_prows,n_pcols,1+np); cla;
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,:,1+nzm);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_cgnx,'XTickLabel',ticklabel_cgnx_); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
axis image;
np = np+1;
end;%for nzm=0:n_zm-1;
figbig;
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17_lpv_A_each_vs_B_max__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17_lpv_A_each_vs_B_umap__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/test_loader_cluster_collect_17_lpv_A_each_vs_B_dex__',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
disp(sprintf(' %% writing %s',fname_fig));
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
assert(nzm==n_zm);

for tmp_infix_ = {'logn','liXf'};
fname_fig = sprintf('%s/test_loader_cluster_collect_17_n_cluster_%s__',dir_jpg,tmp_infix_{1});
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for tmp_infix_ = {'logn','liXf'};

for tmp_infix_ = {'logn','liXf'};
fname_fig = sprintf('%s/test_loader_cluster_collect_17_E_%s_I_%s_n_cluster__',dir_jpg,tmp_infix_{1},tmp_infix_{1});
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for tmp_infix_ = {'logn','liXf'};

%%%%%%%%;
% plot markergene cumulative intersection. ;
%%%%%%%%;
for nX=0:n_X-1;
str_X = str_X_{1+nX};
for nnormalization=0:n_normalization-1;
prefix_normalization = prefix_normalization_{1+nnormalization};
for ncovariate=0:n_covariate-1;
prefix_covariate = prefix_covariate_{1+ncovariate};
n_correction = n_correction___(1+nX,1+nnormalization,1+ncovariate);
disp(sprintf('%s_%s_%s n_correction %d',str_X,prefix_normalization,prefix_covariate,n_correction));
prefix_correction_ = prefix_correction___{1+nX,1+nnormalization,1+ncovariate};
ncorrection_0 = 0;
ncorrection_1 = efind(strcmp(prefix_correction_,sprintf('%s_call',prefix_covariate)));
tmp_markergene_A_auz_0__ = markergene_A___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection_0};
tmp_markergene_A_auz_1__ = markergene_A___{1+nX,1+nnormalization,1+ncovariate}{1+ncorrection_1};
if ((~isempty(tmp_markergene_A_auz_0__)) & (~isempty(tmp_markergene_A_auz_1__)));
fname_fig = sprintf('%s/test_loader_cluster_collect_17_markergene_A_vs_%s_%s_%s_call__',dir_jpg,str_X,prefix_normalization,prefix_covariate);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
n_X_GENE = size(tmp_markergene_A_auz_0__,1);
n_top01 = ceil(n_X_GENE*0.01);
tmp_r_rem_ = transpose(flip(unique([n_top01:n_top01:n_X_GENE,n_X_GENE])));
fr_ = tmp_r_rem_/n_X_GENE; n_fr = numel(fr_);
r_fr__ = zeros(n_label_A,n_fr);
r_nlpv__ = zeros(n_label_A,n_fr);
for nlabel_A=0:n_label_A-1;
[~,tmp_rdrop_0_] = sort(abs(tmp_markergene_A_auz_0__(:,1+nlabel_A)),'ascend');
[~,tmp_rdrop_1_] = sort(abs(tmp_markergene_A_auz_1__(:,1+nlabel_A)),'ascend');
[tmp_r_fr__,tmp_r_pv__,tmp_r_rtn__,tmp_r_mu__,tmp_r_sg__] = cumulative_intersect_general(tmp_r_rem_,tmp_rdrop_0_,tmp_r_rem_,tmp_rdrop_1_);
tmp_r_nlpv__ = log(2) - erfcln((tmp_r_rtn__ - tmp_r_mu__)./(sqrt(2) * tmp_r_sg__));
r_fr__(1+nlabel_A,:) = diag(tmp_r_fr__);
r_nlpv__(1+nlabel_A,:) = diag(tmp_r_nlpv__);
end;%for nlabel_A=0:n_label_A-1;
%%%%%%%%;
figure(1);clf;
nlplim_ = [0,27];
%%%%%%%%;
subplot_(1) = subplot(1,9,1);
frlim_ = [0,1];
c_ = colormap_nlpv(64,0.85,0.15,0.85); n_c = size(c_,1);
tmp_x_ = zeros(4,n_label_A);
tmp_y_ = zeros(4,n_label_A);
tmp_c_ = zeros(1,n_label_A,3);
for nlabel_A=0:n_label_A-1;
tmp_val = r_fr__(1+index_label_A_(1+nlabel_A),end); %<-- not pre-sorted via index_label_A_. ;
nc = max(0,min(n_c-1,floor(n_c*(tmp_val-min(frlim_))/diff(frlim_))));
tmp_y_(:,1+nlabel_A) = -0.5 + 1 + nlabel_A + [0;1;1;0];
tmp_x_(:,1+nlabel_A) = [0;0;tmp_val;tmp_val];
tmp_c_(1,1+nlabel_A,:) = c_(1+nc,:);
end;%for nlabel_A=0:n_label_A-1;
patch(tmp_x_,tmp_y_,tmp_c_,'LineWidth',1);
xlabel('intersection fraction'); ylabel('AIBS cluster label');
xlim(frlim_);
set(gca,'XTick',[0:0.1:max(frlim_)],'XtickLabel',[0:0.1:max(frlim_)]); grid on;
ylim([+0.5,n_label_A+0.5]);
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gca,'Ydir','reverse');
title('intersection fraction (top 1%)');
%%%%%%%%;
subplot_(2) = subplot(1,9,[2,3]);
imagesc(r_fr__(1+index_label_A_,:),[0,1]);
set(gca,'XTick',[0:10:100],'XTickLabel',[0:10:100]); xlabel('percentile');
set(gca,'Ytick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); ylabel('AIBS cluster label');
title('intersection fraction');
colorbar;
%%%%%%%%;
subplot_(3) = subplot(1,9,[4,5]);
imagesc(r_fr__(1+index_label_A_,:) - ones(n_label_A,1)*transpose(fr_.^2),[-1,+1]);
set(gca,'XTick',[0:10:100],'XTickLabel',[0:10:100]); xlabel('percentile');
set(gca,'Ytick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); ylabel('AIBS cluster label');
title('excess intersection fraction');
colorbar;
%%%%%%%%;
subplot_(4) = subplot(1,9,[6,7]);
imagesc(r_nlpv__(1+index_label_A_,:),nlplim_);
set(gca,'XTick',[0:10:100],'XTickLabel',[0:10:100]); xlabel('percentile');
set(gca,'Ytick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); ylabel('AIBS cluster label');
title('-log(p)');
tmp_c_ = colorbar; set(tmp_c_,'YTick',[0:3:max(nlplim_)],'YTickLabel',[0:3:max(nlplim_)]);
%%%%%%%%;
subplot_(5) = subplot(1,9,[8,9]);
imagesc(r_nlpv__(1+index_label_A_,:) - log(n_label_A),nlplim_);
set(gca,'XTick',[0:10:100],'XTickLabel',[0:10:100]); xlabel('percentile');
set(gca,'Ytick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); ylabel('AIBS cluster label');
title('-log(p) (bonf)');
tmp_c_ = colorbar; set(tmp_c_,'YTick',[0:3:max(nlplim_)],'YTickLabel',[0:3:max(nlplim_)]);
%%%%%%%%;
figbig;
colormap(subplot_(2),colormap_nlpv(64,0.85,0.15,0.85));
colormap(subplot_(3),colormap_beach(64));
colormap(subplot_(4),colormap_nlpvt(64));
colormap(subplot_(5),colormap_nlpvt(64));
sgtitle(sprintf('%s_%s_%s',str_X,prefix_normalization,prefix_covariate),'Interpreter','none');
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if ((~isempty(tmp_markergene_A_auz_0__)) & (~isempty(tmp_markergene_A_auz_1__)));
end;%for ncovariate=0:n_covariate-1;
end;%for nnormalization=0:n_normalization-1;
end;%for nX=0:n_X-1;

%%%%%%%%;
% plot A_each_vs_B_each for E_liXf_ and E_liXf_G_call_. ;
%%%%%%%%;
n_E_GENE = n_X_GENE___{1,1,1}(1); n_I_GENE = n_X_GENE___{2,1,1}(1);
for tmp_prefix_cgnx_ = {'E_logn','E_liXf','I_logn','I_liXf'};
%for tmp_prefix_cgnx_ = {'E_liXf','E_liXf_G_call','E_liXf_clogr','E_fill','E_fill_G_call'};
%for tmp_prefix_cgnx_ = {'I_liXf','I_liXf_G_call','I_liXf_clogr','I_fill','I_fill_G_call'};
%for tmp_prefix_cgnx_ = {'E_liXf','E_liXf_G_call','E_liXf_clogr','E_fill','E_fill_G_call','I_liXf','I_liXf_G_call','I_liXf_clogr','I_fill','I_fill_G_call'};
%%%%%%%%;
prefix_cgnx = tmp_prefix_cgnx_{1};
if (~isempty(strfind(prefix_cgnx,'E_'))); n_X_GENE = 1*n_E_GENE + 0*n_I_GENE; end;
if (~isempty(strfind(prefix_cgnx,'I_'))); n_X_GENE = 0*n_E_GENE + 1*n_I_GENE; end;
if (~isempty(strfind(prefix_cgnx,'EI_'))); n_X_GENE = 1*n_E_GENE + 1*n_I_GENE; end;
ncgnx = min(efind(strcmp(prefix_cgnx_,prefix_cgnx)));
assert(strcmp(prefix_cgnx_{1+ncgnx},prefix_cgnx));
%%%%;
nmethod = efind(~cellfun(@isempty,strfind(prefix_method_,'dex'))); %<-- loop. ;
nZ_index = 1; %<-- pos. ;
nzm = efind(strcmp(prefix_zm_,'loop (pos)')); %<-- loop;
str_tmp = sprintf('%s',prefix_cgnx_{1+ncgnx});
fname_fig = sprintf('%s/test_loader_cluster_collect_17_%s',dir_jpg,str_tmp);
tmp_dir_cluster = sprintf('%s/dir_%s_cluster',dir_cluster,prefix_cgnx);
%%%%%%%%;
if (flag_replot | ~exist(sprintf('%s_FIGA.jpg',fname_fig)) | ~exist(sprintf('%s_FIGB.jpg',fname_fig)));
gamma = 0.01; n_shuffle = 64; p_set = 0.05; n_member_lob = 3; flag_force_create_mat = 0;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('tmp_%s',str_xfix); 
dir_out = []; E_array_base_ = zeros(n_u,n_X_GENE); E_array_r_ij_ = []; E_array_c_ij_ = [];
[ ...
 ZRimax_output_label_ ...
,ZRimax_lpFmax_label_ ...
,ZRimax_lpnext_label_ ...
] = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_4( ...
 tmp_dir_cluster ...
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
label_plot_recursive_1(ZRimax_output_label_,ZRimax_lpFmax_label_,ZRimax_lpnext_label_,[]);
figbig;
print('-depsc',sprintf('%s_FIGA.eps',fname_fig));
print('-djpeg',sprintf('%s_FIGA.jpg',fname_fig));
close(gcf);
%%%%%%%%;
u_label_B_ = unique(label_B_);
n_label_B = length(u_label_B_);
label_B_each__ = zeros(n_u,n_label_B);
for nlabel_B=0:n_label_B-1;
label_B_each__(:,1+nlabel_B) = zeros(n_u,1);
tmp_index_ = efind(label_B_==u_label_B_(1+nlabel_B));
label_B_each__(1+tmp_index_,1+nlabel_B) = 1;
end;%for nlabel_B=0:n_label_B-1;
%%%%%%%%;
[tmp_lP0,tmp_cap_] = label_to_label_enrichment_lP0(label_A_,label_B_);
tmp_cap_ = tmp_cap_/n_u;
tmp_cap_A_ = sum(tmp_cap_,2);
tmp_cap_B_ = sum(tmp_cap_,1);
tmp_cap_AB_ = tmp_cap_A_*tmp_cap_B_;
%tmp_cap_fr_ = (tmp_cap_ - tmp_cap_AB_)./tmp_cap_AB_; %<-- excess fraction. ;
tmp_cap_fr_ = tmp_cap_./repmat(tmp_cap_A_,[1,n_label_B]);
%%%%%%%%;
nlpv_each__ = zeros(n_label_A,n_label_B);
for nlabel_A=0:n_label_A-1;
for nlabel_B=0:n_label_B-1;
nlpv_each__(1+nlabel_A,1+nlabel_B) = -label_pair_enrichment(label_A_each__(:,1+nlabel_A),label_B_each__(:,1+nlabel_B));
end;%for nlabel_B=0:n_label_B-1;
end;%for nlabel_A=0:n_label_A-1;
nlpv_each_bonf__ = max(0,nlpv_each__ - log(n_label_A) - log(n_label_B));
%%%%%%%%;
tmp_tab__ = nlpv_each_bonf__(1+index_label_A_,:)>=3;
tmp_index_ = efind(sum(tmp_tab__,1));
[~,tmp_max_ij_] = max(nlpv_each_bonf__(1+index_label_A_,1+tmp_index_),[],1);
[~,index_label_B_] = sort(tmp_max_ij_,'ascend');
index_label_B_ = index_label_B_ - 1;
index_label_B_ = tmp_index_(1+index_label_B_);
index_label_B_ = [index_label_B_,setdiff([0:n_label_B-1],index_label_B_)];
%%%%%%%%;
figure(2);clf;
nlplim_ = [0,27];
%%%%;
subplot_(1) = subplot(1,8,[1,2]);
%imagesc(tmp_cap_fr_(1+index_label_A_,1+index_label_B_),0.25*[-1,+1]);
imagesc(tmp_cap_fr_(1+index_label_A_,1+index_label_B_),[0,1]);
xlabel('unsupervised cluster label'); ylabel('AIBS cluster label');
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gca,'XTick',1:n_label_B,'XTickLabel',u_label_B_(1+index_label_B_)); xtickangle(90);
title('overlap fraction (#/AIBS)');
colorbar;
%%%%;
subplot_(2) = subplot(1,8,[3,4]);
%imagesc(nlpv_each_bonf__(1+index_label_A_,1+index_label_B_),nlplim_);
imagesc(nlpv_each__(1+index_label_A_,1+index_label_B_),nlplim_);
xlabel('unsupervised cluster label'); ylabel('AIBS cluster label');
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gca,'XTick',1:n_label_B,'XTickLabel',u_label_B_(1+index_label_B_)); xtickangle(90);
tmp_c_ = colorbar; set(tmp_c_,'XTick',[0:3:max(nlplim_)]);
title('nlpv');
subplot_(3) = subplot(1,8,[5,6]);
imagesc(nlpv_each_bonf__(1+index_label_A_,1+index_label_B_),nlplim_);
%imagesc(nlpv_each__(1+index_label_A_,1+index_label_B_),nlplim_);
xlabel('unsupervised cluster label'); ylabel('AIBS cluster label');
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gca,'XTick',1:n_label_B,'XTickLabel',u_label_B_(1+index_label_B_)); xtickangle(90);
tmp_c_ = colorbar; set(tmp_c_,'XTick',[0:3:max(nlplim_)]);
title('nlpv (bonf)');
%%%%;
subplot_(4) = subplot(1,8,[7]); hold on;
c_ = colormap_nlpvt(64); n_c = size(c_,1);
tmp_x_ = zeros(4,n_label_A);
tmp_y_ = zeros(4,n_label_A);
tmp_c_ = zeros(1,n_label_A,3);
for nlabel_A=0:n_label_A-1;
tmp_val = nlpv_A_vs_B_Aa_cgnx_zm___(1+nlabel_A,1+ncgnx,1+nzm); %<-- pre sorted via index_label_A_. ;
nc = max(0,min(n_c-1,floor(n_c*(tmp_val-min(nlplim_))/diff(nlplim_))));
tmp_y_(:,1+nlabel_A) = -0.5 + 1 + nlabel_A + [0;1;1;0];
tmp_x_(:,1+nlabel_A) = [0;0;tmp_val;tmp_val];
tmp_c_(1,1+nlabel_A,:) = c_(1+nc,:);
end;%for nlabel_A=0:n_label_A-1;
patch(tmp_x_,tmp_y_,tmp_c_,'LineWidth',1);
xlabel('-log(p)'); ylabel('AIBS cluster label');
xlim(nlplim_);
set(gca,'XTick',[0:3:max(nlplim_)],'XtickLabel',[0:3:max(nlplim_)]); grid on;
ylim([+0.5,n_label_A+0.5]);
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gcf,'Position',1+[0,0,512*4,512*2.0]);
set(gca,'Ydir','reverse');
title('aggregate nlpv');
%%%%;
subplot_(5) = subplot(1,8,[8]); hold on;
c_ = colormap_nlpvt(64); n_c = size(c_,1);
tmp_x_ = zeros(4,n_label_A);
tmp_y_ = zeros(4,n_label_A);
tmp_c_ = zeros(1,n_label_A,3);
for nlabel_A=0:n_label_A-1;
tmp_val = nlpv_A_vs_B_Aa_cgnx_zm___(1+nlabel_A,1+ncgnx,1+nzm) - log(n_label_A); %<-- pre sorted via index_label_A_. ;
nc = max(0,min(n_c-1,floor(n_c*(tmp_val-min(nlplim_))/diff(nlplim_))));
tmp_y_(:,1+nlabel_A) = -0.5 + 1 + nlabel_A + [0;1;1;0];
tmp_x_(:,1+nlabel_A) = [0;0;tmp_val;tmp_val];
tmp_c_(1,1+nlabel_A,:) = c_(1+nc,:);
end;%for nlabel_A=0:n_label_A-1;
patch(tmp_x_,tmp_y_,tmp_c_,'LineWidth',1);
xlabel('-log(p)'); ylabel('AIBS cluster label');
xlim(nlplim_);
set(gca,'XTick',[0:3:max(nlplim_)],'XtickLabel',[0:3:max(nlplim_)]); grid on;
ylim([+0.5,n_label_A+0.5]);
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gcf,'Position',1+[0,0,512*4,512*2.0]);
set(gca,'Ydir','reverse');
title('aggregate nlpv (bonf)');
%%%%;
sgtitle(sprintf('%s',str_tmp),'Interpreter','none');
colormap(subplot_(1),colormap_nlpv(64,0.85,0.15,0.85));
colormap(subplot_(2),colormap_nlpvt(64));
colormap(subplot_(3),colormap_nlpvt(64));
print('-depsc',sprintf('%s_FIGB.eps',fname_fig));
print('-djpeg',sprintf('%s_FIGB.jpg',fname_fig));
close(gcf);
%%%%%%%%;
figure(3);clf;figbig;
nlplim_ = [0,27];
%%%%;
subplot_(1) = subplot(1,4,[1,2,3]);
imagesc(nlpv_each_bonf__(1+index_label_A_,1+index_label_B_),nlplim_);
%imagesc(nlpv_each__(1+index_label_A_,1+index_label_B_),nlplim_);
xlabel('unsupervised cluster label'); ylabel('AIBS cluster label');
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gca,'XTick',1:n_label_B,'XTickLabel',u_label_B_(1+index_label_B_)); xtickangle(90);
tmp_c_ = colorbar; set(tmp_c_,'XTick',[0:3:max(nlplim_)]);
title('nlpv (bonf)');
%%%%;
subplot_(2) = subplot(1,4,[4]); hold on;
c_ = colormap_nlpvt(64); n_c = size(c_,1);
tmp_x_ = zeros(4,n_label_A);
tmp_y_ = zeros(4,n_label_A);
tmp_c_ = zeros(1,n_label_A,3);
for nlabel_A=0:n_label_A-1;
tmp_val = nlpv_A_vs_B_Aa_cgnx_zm___(1+nlabel_A,1+ncgnx,1+nzm) - log(n_label_A); %<-- pre sorted via index_label_A_. ;
nc = max(0,min(n_c-1,floor(n_c*(tmp_val-min(nlplim_))/diff(nlplim_))));
tmp_y_(:,1+nlabel_A) = -0.5 + 1 + nlabel_A + [0;1;1;0];
tmp_x_(:,1+nlabel_A) = [0;0;tmp_val;tmp_val];
tmp_c_(1,1+nlabel_A,:) = c_(1+nc,:);
end;%for nlabel_A=0:n_label_A-1;
patch(tmp_x_,tmp_y_,tmp_c_,'LineWidth',1);
xlabel('-log(p)'); ylabel('AIBS cluster label');
xlim(nlplim_);
set(gca,'XTick',[0:3:max(nlplim_)],'XtickLabel',[0:3:max(nlplim_)]); grid on;
ylim([+0.5,n_label_A+0.5]);
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gcf,'Position',1+[0,0,512*4,512*2.0]);
set(gca,'Ydir','reverse');
title('aggregate nlpv (bonf)');
%%%%;
sgtitle(sprintf('%s',str_tmp),'Interpreter','none');
colormap(subplot_(1),colormap_nlpvt(64));
print('-depsc',sprintf('%s_FIGC.eps',fname_fig));
print('-djpeg',sprintf('%s_FIGC.jpg',fname_fig));
close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s_FIGA.jpg',fname_fig)) | ~exist(sprintf('%s_FIGB.jpg',fname_fig)));
end;%for tmp_prefix_cgnx_ = {'E_liXf','E_liXf_G_call','E_liXf_clogr','E_liXf_clogr_G_call','E_fill','E_fill_G_call'};

%%%%%%%%;
% Summary figure. ;
%%%%%%%%;
fname_fig = sprintf('%s/test_loader_cluster_collect_17_Summary_FIG0',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
nzm = efind(strcmp(prefix_zm_,'loop (pos)'));
%%%%%%%%;
figure(1);clf;
%%%%%%%%;
tmp_prefix_cgnx_ = { ...
 'E_logn' ...
,'E_logn_G_call' ...
,'E_logn_G_c1' ...
,'E_logn_G_c2' ...
,'E_liXf' ...
,'E_liXf_G_call' ...
,'E_liXf_G_c1' ...
,'E_liXf_G_c2' ...
,'E_liXf_clogr' ...
,'E_liXf_clogr_G_call' ...
,'E_liXf_clogr_G_c1' ...
,'E_liXf_clogr_G_c2' ...
,'E_fill' ...
,'E_fill_G_call' ...
,'E_fill_G_c1' ...
,'E_fill_G_c2' ...
};
n_tmp_prefix_cgnx = numel(tmp_prefix_cgnx_);
tmp_index_cgnx_ = zeros(n_tmp_prefix_cgnx,1);
for nl=0:n_tmp_prefix_cgnx-1;
tmp_index_cgnx_(1+nl) = min(efind(strcmp(prefix_cgnx_,tmp_prefix_cgnx_{1+nl})));
end;%for nl=0:n_tmp_prefix_cgnx-1;
%%%%
subplot_(1) = subplot(1,4,1);
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,1+tmp_index_cgnx_,1+nzm) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_tmp_prefix_cgnx,'XTickLabel',ticklabel_cgnx_(1+tmp_index_cgnx_)); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s -log(p) (bonf)',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
%%%%;
subplot_(2) = subplot(1,4,2);
tmp_lim_ = [0,1]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = markergene_top01_A_cgnx__(1+index_label_A_,1+tmp_index_cgnx_);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_));
set(gca,'XTick',1:n_tmp_prefix_cgnx,'XTickLabel',ticklabel_cgnx_(1+tmp_index_cgnx_)); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title('top01');
set(gca,'FontSize',6);
%%%%%%%%;
tmp_prefix_cgnx_ = { ...
,'I_liXf' ...
,'I_liXf_G_call' ...
,'I_liXf_G_c1' ...
,'I_liXf_G_c2' ...
,'I_liXf_clogr' ...
,'I_liXf_clogr_G_call' ...
,'I_liXf_clogr_G_c1' ...
,'I_liXf_clogr_G_c2' ...
,'I_fill' ...
,'I_fill_G_call' ...
,'I_fill_G_c1' ...
,'I_fill_G_c2' ...
};
n_tmp_prefix_cgnx = numel(tmp_prefix_cgnx_);
tmp_index_cgnx_ = zeros(n_tmp_prefix_cgnx,1);
for nl=0:n_tmp_prefix_cgnx-1;
tmp_index_cgnx_(1+nl) = min(efind(strcmp(prefix_cgnx_,tmp_prefix_cgnx_{1+nl})));
end;%for nl=0:n_tmp_prefix_cgnx-1;
%%%%;
subplot_(3) = subplot(1,4,3);
c_ = colormap_nlpvt(); n_c = size(c_,1); c_(1,:) = [1,1,1]; colormap(c_);
tmp_lim_ = [0,27]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = nlpv_A_vs_B_Aa_cgnx_zm___(:,1+tmp_index_cgnx_,1+nzm) - log(n_label_A);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A+1,'YTickLabel',ticklabel_Aa_);
set(gca,'XTick',1:n_tmp_prefix_cgnx,'XTickLabel',ticklabel_cgnx_(1+tmp_index_cgnx_)); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title(sprintf('%s -log(p) (bonf)',prefix_zm_{1+nzm}),'Interpreter','none');
set(gca,'FontSize',6);
%%%%;
subplot_(4) = subplot(1,4,4);
tmp_lim_ = [0,1]; tmp_min = min(tmp_lim_) + diff(tmp_lim_)/n_c;
tmp_ = markergene_top01_A_cgnx__(1+index_label_A_,1+tmp_index_cgnx_);
tmp_(find( isfinite(tmp_))) = max(tmp_min,tmp_(find( isfinite(tmp_))));
tmp_(find(~isfinite(tmp_))) = min(tmp_lim_);
imagesc(tmp_,tmp_lim_);colorbar;
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_));
set(gca,'XTick',1:n_tmp_prefix_cgnx,'XTickLabel',ticklabel_cgnx_(1+tmp_index_cgnx_)); xtickangle(90);
set(gca,'TickLength',[0,0]);
ylabel('AIBS cluster label');
xlabel('data normalization+correction');
title('top01');
set(gca,'FontSize',6);
%%%%%%%%;
figbig;
colormap(subplot_(1),colormap_nlpvt(64));
colormap(subplot_(2),colormap_nlpv(64,0.85,0.15,0.85));
colormap(subplot_(3),colormap_nlpvt(64));
colormap(subplot_(4),colormap_nlpv(64,0.85,0.15,0.85));
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Attempt to visualize clusters drifting and merging as a result of covariate-correction. ;
% Using 7 clusters, since there are 7 distinct clusters in colormap('lines');
%%%%%%%%;
%str_label_A_sub_ = {'1','6','7'}; %<-- not bad ;
%str_label_A_sub_ = {'1','2','3'}; %<-- some covariate-dependent separation ;
%str_label_A_sub_ = {'4','5','6','7'}; %<-- some dependence on one covariate group (aob 128) ;
%str_label_A_sub_ = {'4','5','10','11'}; %<-- some dependence on another covariate group. ;
str_label_A_sub_ = {'1','4','5'}; %<-- quite clear. ;
n_label_A_sub = numel(str_label_A_sub_);
str_label_A_sub_prefix_ = [];
index_label_A_sub_ = [];
for nlabel_A_sub=0:n_label_A_sub-1;
str_label_A_sub_prefix_ = sprintf('%s_%s',str_label_A_sub_prefix_,str_label_A_sub_{1+nlabel_A_sub});
index_label_A_sub_ = union(index_label_A_sub_,efind(strcmp(str_label_A_{1},str_label_A_sub_{1+nlabel_A_sub})));
end;%for nlabel_A_sub=0:n_label_A_sub-1;
disp(sprintf(' %% str_label_A_sub_: '));
disp(u_label_A_(unique(label_A_(1+index_label_A_sub_))));
label_A_sub_ = label_A_(1+index_label_A_sub_);
%%%%%%%%;
c__ = [ ...
 0.45,0.45,0.45 ...
;1,0.5,0.5 ...
;0,1,0 ...
;0,0,1 ...
;0,1,1 ...
;1,0,1 ...
;1,1,0 ...
;0.85,0.85,0.85 ...
];
%%%%%%%%;
str_X = 'E';
if (strcmp(str_X,'E')); rank_estimate_X = rank_estimate_E; end;
if (strcmp(str_X,'I')); rank_estimate_X = rank_estimate_I; end;
prefix_normalization = 'liXf_clogr';
prefix_covariate = 'G';
tmp_fname_tsv = sprintf('%s/n_QCluster_%s_original.tsv',dir_mat,prefix_covariate);
n_QCluster = textread(tmp_fname_tsv); n_QCluster = n_QCluster(1);
disp(sprintf(' %% %s n_QCluster %d',prefix_covariate,n_QCluster));
%%%%%%%%;
if (strcmp(str_X,'E') & ~isempty(strfind(prefix_normalization,'fill'))); X_ = E_fill_; end;
if (strcmp(str_X,'E') & ~isempty(strfind(prefix_normalization,'li16f'))); X_ = E_li16f_; end;
if (strcmp(str_X,'E') & ~isempty(strfind(prefix_normalization,'liXf'))); X_ = E_liXf_; end;
if (strcmp(str_X,'I') & ~isempty(strfind(prefix_normalization,'fill'))); X_ = I_fill_; end;
if (strcmp(str_X,'I') & ~isempty(strfind(prefix_normalization,'li16f'))); X_ = I_li16f_; end;
if (strcmp(str_X,'I') & ~isempty(strfind(prefix_normalization,'liXf'))); X_ = I_liXf_; end;
if (strcmp(str_X,'EI') & ~isempty(strfind(prefix_normalization,'fill'))); X_ = [E_fill_ , I_fill_]; end;
if (strcmp(str_X,'EI') & ~isempty(strfind(prefix_normalization,'li16f'))); X_ = [E_li16f_ , I_li16f_]; end;
if (strcmp(str_X,'EI') & ~isempty(strfind(prefix_normalization,'liXf'))); X_ = [E_liXf_ , I_liXf_]; end;
if (~isempty(strfind(prefix_normalization,'clogr'))); X_ = mean_center_0(X_,'col'); end;
Y00_ = X_; %<-- uncorrected backup. ;
[tmp_UY00_,tmp_SY00_,tmp_VY00_] = svds([Y00_(1+index_label_A_sub_,:)],rank_estimate_X);
Y00_sub_ = Y00_(1+index_label_A_sub_,:)*tmp_VY00_;
%%%%%%%%;
nQCluster = -1; %<-- -1 corresponds to all, nQCluster range in 0:n_QCluster-1;
tmp_H_ = [];
if ((nQCluster>=0) & (nQCluster< n_QCluster));
tmp_fname_mat = sprintf('%s/S_c%d_%s_original.mat',dir_mat,1+nQCluster,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate_mcH = tmp_.rank_estimate_sample; clear tmp_;
tmp_fname_tsv = sprintf('%s/c%d_%s__.tsv',dir_mat,1+nQCluster,prefix_covariate);
tmp_H_ = textread(tmp_fname_tsv); tmp_H_ = tmp_H_(:,1:end-1);
str_c = sprintf('%d',1+nQCluster);
end;%if ((nQCluster>=0) & (nQCluster< n_QCluster));
if (((nQCluster< 0) | (nQCluster>=n_QCluster)) & (n_QCluster> 1));
tmp_fname_mat = sprintf('%s/S_call_%s_original.mat',dir_mat,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate_mcH = tmp_.rank_estimate_sample; clear tmp_;
if (strcmp(prefix_covariate,'C')); tmp_H_ = C_rank_; end;
if (strcmp(prefix_covariate,'G')); tmp_H_ = G_gaus_; end;
str_c = 'all';
end;%if (((nQCluster< 0) | (nQCluster>=n_QCluster)) & (n_QCluster> 1));
%%%%;
rank_estimate_mcH_use = min(size(tmp_H_,2),max(1,rank_estimate_mcH)); %<-- use at least one rank. ;
[tmp_UmcH_,tmp_SmcH_,tmp_VmcH_] = svds(mean_center_0(tmp_H_),rank_estimate_mcH_use);
tmp_mcH_ = tmp_UmcH_*tmp_SmcH_*transpose(tmp_VmcH_);
tmp_1mcH_ = [ ones(n_u,1) , tmp_mcH_ ];
[tmp_U1mcH_,tmp_S1mcH_,tmp_V1mcH_] = svds(tmp_1mcH_,1+rank_estimate_mcH_use); %<-- allow for ones-vector to add a rank. ;
Y11_ = Y00_ - tmp_U1mcH_*transpose(tmp_U1mcH_)*Y00_; %<-- residual. ;
Y11_sub_ = Y11_(1+index_label_A_sub_,:)*tmp_VY00_; clear Y11_;
%%%%%%%%;
nQCluster = 0; %<-- -1 corresponds to all, nQCluster range in 0:n_QCluster-1;
tmp_H_ = [];
if ((nQCluster>=0) & (nQCluster< n_QCluster));
tmp_fname_mat = sprintf('%s/S_c%d_%s_original.mat',dir_mat,1+nQCluster,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate_mcH = tmp_.rank_estimate_sample; clear tmp_;
tmp_fname_tsv = sprintf('%s/c%d_%s__.tsv',dir_mat,1+nQCluster,prefix_covariate);
tmp_H_ = textread(tmp_fname_tsv); tmp_H_ = tmp_H_(:,1:end-1);
str_c = sprintf('%d',1+nQCluster);
end;%if ((nQCluster>=0) & (nQCluster< n_QCluster));
if (((nQCluster< 0) | (nQCluster>=n_QCluster)) & (n_QCluster> 1));
tmp_fname_mat = sprintf('%s/S_call_%s_original.mat',dir_mat,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate_mcH = tmp_.rank_estimate_sample; clear tmp_;
if (strcmp(prefix_covariate,'C')); tmp_H_ = C_rank_; end;
if (strcmp(prefix_covariate,'G')); tmp_H_ = G_gaus_; end;
str_c = 'all';
end;%if (((nQCluster< 0) | (nQCluster>=n_QCluster)) & (n_QCluster> 1));
%%%%;
rank_estimate_mcH_use = min(size(tmp_H_,2),max(1,rank_estimate_mcH)); %<-- use at least one rank. ;
[tmp_UmcH_,tmp_SmcH_,tmp_VmcH_] = svds(mean_center_0(tmp_H_),rank_estimate_mcH_use);
tmp_mcH_ = tmp_UmcH_*tmp_SmcH_*transpose(tmp_VmcH_);
tmp_1mcH_ = [ ones(n_u,1) , tmp_mcH_ ];
[tmp_U1mcH_,tmp_S1mcH_,tmp_V1mcH_] = svds(tmp_1mcH_,1+rank_estimate_mcH_use); %<-- allow for ones-vector to add a rank. ;
Y01_ = Y00_ - tmp_U1mcH_*transpose(tmp_U1mcH_)*Y00_; %<-- residual. ;
Y01_sub_ = Y01_(1+index_label_A_sub_,:)*tmp_VY00_; clear Y01_;
%%%%%%%%;
nQCluster = 1; %<-- -1 corresponds to all, nQCluster range in 0:n_QCluster-1;
tmp_H_ = [];
if ((nQCluster>=0) & (nQCluster< n_QCluster));
tmp_fname_mat = sprintf('%s/S_c%d_%s_original.mat',dir_mat,1+nQCluster,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate_mcH = tmp_.rank_estimate_sample; clear tmp_;
tmp_fname_tsv = sprintf('%s/c%d_%s__.tsv',dir_mat,1+nQCluster,prefix_covariate);
tmp_H_ = textread(tmp_fname_tsv); tmp_H_ = tmp_H_(:,1:end-1);
str_c = sprintf('%d',1+nQCluster);
end;%if ((nQCluster>=0) & (nQCluster< n_QCluster));
if (((nQCluster< 0) | (nQCluster>=n_QCluster)) & (n_QCluster> 1));
tmp_fname_mat = sprintf('%s/S_call_%s_original.mat',dir_mat,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate_mcH = tmp_.rank_estimate_sample; clear tmp_;
if (strcmp(prefix_covariate,'C')); tmp_H_ = C_rank_; end;
if (strcmp(prefix_covariate,'G')); tmp_H_ = G_gaus_; end;
str_c = 'all';
end;%if (((nQCluster< 0) | (nQCluster>=n_QCluster)) & (n_QCluster> 1));
%%%%;
rank_estimate_mcH_use = min(size(tmp_H_,2),max(1,rank_estimate_mcH)); %<-- use at least one rank. ;
[tmp_UmcH_,tmp_SmcH_,tmp_VmcH_] = svds(mean_center_0(tmp_H_),rank_estimate_mcH_use);
tmp_mcH_ = tmp_UmcH_*tmp_SmcH_*transpose(tmp_VmcH_);
tmp_1mcH_ = [ ones(n_u,1) , tmp_mcH_ ];
[tmp_U1mcH_,tmp_S1mcH_,tmp_V1mcH_] = svds(tmp_1mcH_,1+rank_estimate_mcH_use); %<-- allow for ones-vector to add a rank. ;
Y10_ = Y00_ - tmp_U1mcH_*transpose(tmp_U1mcH_)*Y00_; %<-- residual. ;
Y10_sub_ = Y10_(1+index_label_A_sub_,:)*tmp_VY00_; clear Y10_;
%%%%%%%%;

%%%%%%%%;
flag_2d_vs_3d = 1;
for laob = 0:1:5;
if (laob>=0); str_laob = sprintf('_p%d',+laob); end;
if (laob< 0); str_laob = sprintf('_n%d',-laob); end;
fname_fig = sprintf('%s/test_loader_cluster_collect_17_path_A_label_sub%s_laob%s',dir_jpg,str_label_A_sub_prefix_,str_laob);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig)));
disp(sprintf(' %% %s not found, creating',sprintf('%s.jpg',fname_fig)));
figure(1);clf;figbig;colormap('lines');
aob = 2^laob;
[H00] = SVD_discat_0([Y00_sub_],label_A_sub_,3,aob);
[H01] = SVD_discat_0([Y01_sub_],label_A_sub_,3,aob);
[H10] = SVD_discat_0([Y10_sub_],label_A_sub_,3,aob);
[H11] = SVD_discat_0([Y11_sub_],label_A_sub_,3,aob);
[tmp_VL_,tmp_SL_] = eigs([H00],3,'largestreal'); %<-- sum the hessians if desired. ;
Y00VL_ = Y00_sub_*tmp_VL_;
Y01VL_ = Y01_sub_*tmp_VL_;
Y10VL_ = Y10_sub_*tmp_VL_;
Y11VL_ = Y11_sub_*tmp_VL_;
%%%%%%%%;
tmp_ = [Y00VL_(:,1+0);Y01VL_(:,1+0);Y10VL_(:,1+0);Y11VL_(:,1+0)]; Y_lim_0_ = mean(tmp_) + std(tmp_,1)*1.5*[-1,+1]; clear tmp_;
tmp_ = [Y00VL_(:,1+1);Y01VL_(:,1+1);Y10VL_(:,1+1);Y11VL_(:,1+1)]; Y_lim_1_ = mean(tmp_) + std(tmp_,1)*1.5*[-1,+1]; clear tmp_;
tmp_ = [Y00VL_(:,1+2);Y01VL_(:,1+2);Y10VL_(:,1+2);Y11VL_(:,1+2)]; Y_lim_2_ = mean(tmp_) + std(tmp_,1)*1.5*[-1,+1]; clear tmp_;
Y_lim_equal_ = [min([Y_lim_0_,Y_lim_1_,Y_lim_2_]),max([Y_lim_0_,Y_lim_1_,Y_lim_2_])];
%%%%%%%%;
tmp_ = [Y00VL_(:,1+0);Y01VL_(:,1+0);Y10VL_(:,1+0);Y11VL_(:,1+0)]; Y_lim_0_ = [min(tmp_),max(tmp_)];
tmp_ = [Y00VL_(:,1+1);Y01VL_(:,1+1);Y10VL_(:,1+1);Y11VL_(:,1+1)]; Y_lim_1_ = [min(tmp_),max(tmp_)];
tmp_ = [Y00VL_(:,1+2);Y01VL_(:,1+2);Y10VL_(:,1+2);Y11VL_(:,1+2)]; Y_lim_2_ = [min(tmp_),max(tmp_)];
Y_lim_equal_ = [min([Y_lim_0_,Y_lim_1_,Y_lim_2_]),max([Y_lim_0_,Y_lim_1_,Y_lim_2_])];
%%%%%%%%;
if (flag_2d_vs_3d==1);
markersize_use = 35;
subplot(2,2,1); ZVL_ = Y00VL_;
scatter(ZVL_(:,1),ZVL_(:,2),markersize_use,c__(label_num_to_enum_0(label_A_sub_),:),'filled','MarkerFaceAlpha',0.67,'MarkerEdgeColor','k');
%xlim(Y_lim_equal_); ylim(Y_lim_equal_); grid on; xlabel('pc1'); ylabel('pc2');
xlim(Y_lim_0_); ylim(Y_lim_1_); grid on; xlabel('pc1'); ylabel('pc2');
axis equal;
  title(sprintf('Original Data'));
subplot(2,2,2); ZVL_ = Y01VL_;
scatter(ZVL_(:,1),ZVL_(:,2),markersize_use,c__(label_num_to_enum_0(label_A_sub_),:),'filled','MarkerFaceAlpha',0.67,'MarkerEdgeColor','k');
%xlim(Y_lim_equal_); ylim(Y_lim_equal_); grid on; xlabel('pc1'); ylabel('pc2');
xlim(Y_lim_0_); ylim(Y_lim_1_); grid on; xlabel('pc1'); ylabel('pc2');
axis equal;
  title(sprintf('Corrected for QC-cluster-1'));
subplot(2,2,3); ZVL_ = Y10VL_;
scatter(ZVL_(:,1),ZVL_(:,2),markersize_use,c__(label_num_to_enum_0(label_A_sub_),:),'filled','MarkerFaceAlpha',0.67,'MarkerEdgeColor','k');
%xlim(Y_lim_equal_); ylim(Y_lim_equal_); grid on; xlabel('pc1'); ylabel('pc2');
xlim(Y_lim_0_); ylim(Y_lim_1_); grid on; xlabel('pc1'); ylabel('pc2');
axis equal;
  title(sprintf('Corrected for QC-cluster-2'));
subplot(2,2,4); ZVL_ = Y11VL_;
scatter(ZVL_(:,1),ZVL_(:,2),markersize_use,c__(label_num_to_enum_0(label_A_sub_),:),'filled','MarkerFaceAlpha',0.67,'MarkerEdgeColor','k');
%xlim(Y_lim_equal_); ylim(Y_lim_equal_); grid on; xlabel('pc1'); ylabel('pc2');
xlim(Y_lim_0_); ylim(Y_lim_1_); grid on; xlabel('pc1'); ylabel('pc2');
axis equal;
  title(sprintf('Corrected for both QC-clusters'));
end;%if (flag_2d_vs_3d==0);
%%%%%%%%;
if (flag_2d_vs_3d==0);
markersize_use = 15;
subplot(2,2,1); ZVL_ = Y00VL_;
scatter3(ZVL_(:,1),ZVL_(:,2),ZVL_(:,3),markersize_use,c__(label_num_to_enum_0(label_A_sub_),:),'filled');
xlim(Y_lim_equal_); ylim(Y_lim_equal_); zlim(Y_lim_equal_); grid on;
axis vis3d; axis equal;  title(sprintf('Original Data'));
subplot(2,2,2); ZVL_ = Y01VL_;
scatter3(ZVL_(:,1),ZVL_(:,2),ZVL_(:,3),markersize_use,c__(label_num_to_enum_0(label_A_sub_),:),'filled');
xlim(Y_lim_equal_); ylim(Y_lim_equal_); zlim(Y_lim_equal_); grid on;
axis vis3d; axis equal;  title(sprintf('Corrected for QC-cluster-1'));
subplot(2,2,3); ZVL_ = Y10VL_;
scatter3(ZVL_(:,1),ZVL_(:,2),ZVL_(:,3),markersize_use,c__(label_num_to_enum_0(label_A_sub_),:),'filled');
xlim(Y_lim_equal_); ylim(Y_lim_equal_); zlim(Y_lim_equal_); grid on;
axis vis3d; axis equal;  title(sprintf('Corrected for QC-cluster-2'));
subplot(2,2,4); ZVL_ = Y11VL_;
scatter3(ZVL_(:,1),ZVL_(:,2),ZVL_(:,3),markersize_use,c__(label_num_to_enum_0(label_A_sub_),:),'filled');
xlim(Y_lim_equal_); ylim(Y_lim_equal_); zlim(Y_lim_equal_); grid on;
axis vis3d; axis equal;  title(sprintf('Corrected for both QC-clusters'));
end;%if (flag_2d_vs_3d==0);
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig)));
end;%for laob = -3:+3;

end;%if flag_plot_make;
