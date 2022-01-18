function test_dexcluster_multi_4(M,N,snr,n_cluster,n_rank,n_iteration,n_shuffle,flag_rerun);
% tries to find up to 8 clusters. ;
% single rank. ;
% see test_dexcluster_multi_4_collect.m ;

if (nargin<1);
for nl=0:5;
na=0;    
if (nl==na); M =  178; N =  2e3; n_cluster = 1; n_rank =  2; n_iteration = 8; n_shuffle = 256; end; na=na+1;
if (nl==na); M =  178; N =  2e3; n_cluster = 3; n_rank =  4; n_iteration = 8; n_shuffle = 256; end; na=na+1;
if (nl==na); M =  563; N = 6325; n_cluster = 1; n_rank =  3; n_iteration = 4; n_shuffle =  64; end; na=na+1;
if (nl==na); M =  563; N = 6325; n_cluster = 6; n_rank =  9; n_iteration = 4; n_shuffle =  64; end; na=na+1;
if (nl==na); M = 1781; N =  2e4; n_cluster = 1; n_rank =  4; n_iteration = 4; n_shuffle =  64; end; na=na+1;
if (nl==na); M = 1781; N =  2e4; n_cluster = 8; n_rank = 12; n_iteration = 4; n_shuffle =  64; end; na=na+1;
flag_rerun = 0;
n_step = 21;
snr_ = linspace(0.5,1.0,n_step);
for nstep=1:n_step;
snr = snr_(nstep);
test_dexcluster_multi_4(M,N,snr,n_cluster,n_rank,n_iteration,n_shuffle,flag_rerun);
end;%for nstep=1:n_step;
end;%for nl=0:5;
disp('returning'); return;
end;%if (nargin<1);

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jamison',string_root);
dir_base = sprintf('%s/dir_dexcluster_multi_4',dir_trunk);
n_rank_sub = 1;

if ( libisloaded('halfloop_lib'));
unloadlibrary('halfloop_lib');
end;%if ( libisloaded('halfloop_lib'));
if (~libisloaded('halfloop_lib'));
tmp_dir = pwd;
cd(sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root));
system('make -f halfloop_OptiPlex.make halfloop_lib;');
loadlibrary('halfloop_lib','halfloop_lib.h','includepath',sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root));
cd(tmp_dir);
end;%if (~libisloaded('halfloop_lib'));

if (~exist(dir_base,'dir')); disp(sprintf(' %% mkdir %s',dir_base)); mkdir(dir_base); end;
%%%%%%%%;
for niteration=1:n_iteration;
disp(sprintf(' %% niteration %d/%d',niteration,n_iteration));
%%%%%%%%;
rng(niteration);
[A_n_,label_A_,n_label_A_,pf_,pi_,snr_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);
[tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank);
UVS_n_ = tmp_U_*tmp_S_*transpose(tmp_V_);
US_n_ = tmp_U_*tmp_S_;
clear tmp_U_ tmp_S_ tmp_V_; 
str_xfix = test_dexcluster_multi_xfix_4('base',M,N,snr,n_cluster,n_rank,niteration);
str_base = sprintf('%s/%s.mat',dir_base,str_xfix);
if (~exist(str_base,'file')); 
disp(sprintf(' %% %s not found, creating',str_base)); 
save(str_base,'n_rank','M','N','snr','label_A_','n_label_A_','pf_','pi_','n_cluster');
end;%if (~exist(str_base,'file')); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% no preliminary projection. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne00 --> isosplit5 ;
str_xfix = test_dexcluster_multi_xfix_4('tsne00_isosplit5',M,N,snr,n_cluster,[],niteration);
str_tsne00_isosplit5_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_tsne00_isosplit5_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_tsne00_isosplit5_mat,'file') & ~exist(str_tsne00_isosplit5_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_tsne00_isosplit5_mat)); 
save(str_tsne00_isosplit5_tmp,'str_tsne00_isosplit5_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
A_n_sub_ = fast_tsne_dr_0(A_n_,struct('rand_seed',1,'no_dims',2,'theta',0));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_tsne00_isosplit5_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_tsne00_isosplit5_mat)); end;%try;
delete(str_tsne00_isosplit5_tmp);
end;%if (~exist(str_tsne00_isosplit5,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne50 --> isosplit5 ;
str_xfix = test_dexcluster_multi_xfix_4('tsne50_isosplit5',M,N,snr,n_cluster,[],niteration);
str_tsne50_isosplit5_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_tsne50_isosplit5_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_tsne50_isosplit5_mat,'file') & ~exist(str_tsne50_isosplit5_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_tsne50_isosplit5_mat)); 
save(str_tsne50_isosplit5_tmp,'str_tsne50_isosplit5_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
A_n_sub_ = fast_tsne_dr_0(A_n_,struct('rand_seed',1,'no_dims',2,'theta',0.5));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_tsne50_isosplit5_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_tsne50_isosplit5_mat)); end;%try;
delete(str_tsne50_isosplit5_tmp);
end;%if (~exist(str_tsne50_isosplit5,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> default ;
str_xfix = test_dexcluster_multi_xfix_4('umap00_default',M,N,snr,n_cluster,[],niteration);
str_umap00_default_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_umap00_default_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_umap00_default_mat,'file') & ~exist(str_umap00_default_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_umap00_default_mat)); 
save(str_umap00_default_tmp,'str_umap00_default_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
[A_n_sub_, ~, label_B__{nrank_sub}]=run_umap(A_n_,'verbose','none');
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_umap00_default_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_umap00_default_mat)); end;%try;
delete(str_umap00_default_tmp);
end;%if (~exist(str_umap00_default,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> isosplit5 ;
str_xfix = test_dexcluster_multi_xfix_4('umap00_isosplit5',M,N,snr,n_cluster,[],niteration);
str_umap00_isosplit5_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_umap00_isosplit5_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_umap00_isosplit5_mat,'file') & ~exist(str_umap00_isosplit5_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_umap00_isosplit5_mat)); 
save(str_umap00_isosplit5_tmp,'str_umap00_isosplit5_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
[A_n_sub_, ~, ~]=run_umap(A_n_,'verbose','none');
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_umap00_isosplit5_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_umap00_isosplit5_mat)); end;%try;
delete(str_umap00_isosplit5_tmp);
end;%if (~exist(str_umap00_isosplit5,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> hdbscan ;
str_xfix = test_dexcluster_multi_xfix_4('umap00_hdbscan',M,N,snr,n_cluster,[],niteration);
str_umap00_hdbscan_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_umap00_hdbscan_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_umap00_hdbscan_mat,'file') & ~exist(str_umap00_hdbscan_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_umap00_hdbscan_mat)); 
save(str_umap00_hdbscan_tmp,'str_umap00_hdbscan_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
[A_n_sub_, ~, ~]=run_umap(A_n_,'verbose','none');
hdbscan_ = HDBSCAN(A_n_sub_); 
hdbscan_.minpts = 10; hdbscan_.fit_model(); hdbscan_.get_best_clusters(); hdbscan_.get_membership();
label_B__{nrank_sub} = hdbscan_.labels;
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_umap00_hdbscan_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_umap00_hdbscan_mat)); end;%try;
delete(str_umap00_hdbscan_tmp);
end;%if (~exist(str_umap00_hdbscan,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% louvain00 --> default;
str_xfix = test_dexcluster_multi_xfix_4('louvain00_default',M,N,snr,n_cluster,[],niteration);
str_louvain00_default_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_louvain00_default_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_louvain00_default_mat,'file') & ~exist(str_louvain00_default_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_louvain00_default_mat)); 
save(str_louvain00_default_tmp,'str_louvain00_default_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
label_B__{nrank_sub} = GCModulMax1(A_n_*transpose(A_n_));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_louvain00_default_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_louvain00_default_mat)); end;%try;
delete(str_louvain00_default_tmp);
end;%if (~exist(str_louvain00_default,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%{
% dexcluster --> ZRimax ;
str_xfix = test_dexcluster_multi_xfix_4('dexcluster_nonbinary_trace_ZRimax',M,N,snr,n_cluster,[],niteration);
str_dexcluster_nonbinary_trace_ZRimax_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_dexcluster_nonbinary_trace_ZRimax_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_dexcluster_nonbinary_trace_ZRimax_mat,'file') & ~exist(str_dexcluster_nonbinary_trace_ZRimax_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_dexcluster_nonbinary_trace_ZRimax_mat)); 
save(str_dexcluster_nonbinary_trace_ZRimax_tmp,'str_dexcluster_nonbinary_trace_ZRimax_mat');
%try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
prefix_base = sprintf('tmp_%s_r%d',str_xfix,nrank_sub); 
dir_out = []; E_array_base_ = A_n_; E_array_r_ij_ = []; E_array_c_ij_ = []; gamma = 0.01; p_set = []; n_member_lob = 2;
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_4( ...
 dir_base ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
);
ZRimax_label_ = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_4( ...
 dir_base ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
);
label_B_1_ = label_str_to_enum_0(ZRimax_label_); 
label_B_2_ = label_rearrange_0(label_B_1_,E_array_base_);
[lpv_2,lP0_2,fla_2] = label_to_label_enrichment_quad_4(label_A_,label_B_2_);
lpv_(nrank_sub) = lpv_2; lP0_(nrank_sub) = lP0_2; fla_(nrank_sub) = fla_2; label_B__{nrank_sub} = label_B_2_;
try; rmdir(sprintf('%s/dir_%s',dir_base,prefix_base),'s'); catch; disp(sprintf(' %% could not remove %s',sprintf('%s/dir_%s',dir_base,prefix_base))); end;%try;
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_dexcluster_nonbinary_trace_ZRimax_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
%catch; disp(sprintf(' %% WARNING: error generating %s',str_dexcluster_nonbinary_trace_ZRimax_mat)); end;%try;
delete(str_dexcluster_nonbinary_trace_ZRimax_tmp);
end;%if (~exist(str_dexcluster_nonbinary_trace_ZRimax,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% dexnbtZRgumb ;
str_xfix = test_dexcluster_multi_xfix_4('dexnbtZRgumb',M,N,snr,n_cluster,[],niteration);
str_dexnbtZRgumb_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_dexnbtZRgumb_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_dexnbtZRgumb_mat,'file') & ~exist(str_dexnbtZRgumb_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_dexnbtZRgumb_mat)); 
save(str_dexnbtZRgumb_tmp,'str_dexnbtZRgumb_mat');
%try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
prefix_base = sprintf('tmp_%s_r%d',str_xfix,nrank_sub); 
dir_out = []; E_array_base_ = A_n_; E_array_r_ij_ = []; E_array_c_ij_ = []; gamma = 0.01; p_set = []; n_member_lob = 2;
dexnbtZRgumb_recursive_5( ...
 dir_base ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
);
ZRimax_label_ = ...
dexnbtZRgumb_recursive_5( ...
 dir_base ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
);
label_B_1_ = label_str_to_enum_0(ZRimax_label_); 
label_B_2_ = label_rearrange_0(label_B_1_,E_array_base_);
[lpv_2,lP0_2,fla_2] = label_to_label_enrichment_quad_4(label_A_,label_B_2_);
lpv_(nrank_sub) = lpv_2; lP0_(nrank_sub) = lP0_2; fla_(nrank_sub) = fla_2; label_B__{nrank_sub} = label_B_2_;
try; rmdir(sprintf('%s/dir_%s',dir_base,prefix_base),'s'); catch; disp(sprintf(' %% could not remove %s',sprintf('%s/dir_%s',dir_base,prefix_base))); end;%try;
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_dexnbtZRgumb_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
%catch; disp(sprintf(' %% WARNING: error generating %s',str_dexnbtZRgumb_mat)); end;%try;
delete(str_dexnbtZRgumb_tmp);
end;%if (~exist(str_dexnbtZRgumb,'file')); 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% hnbtZRgumb
str_xfix = test_dexcluster_multi_xfix_4('hnbtZRgumb',M,N,snr,n_cluster,[],niteration);
str_hnbtZRgumb_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_hnbtZRgumb_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_hnbtZRgumb_mat,'file') & ~exist(str_hnbtZRgumb_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_hnbtZRgumb_mat)); 
save(str_hnbtZRgumb_tmp,'str_hnbtZRgumb_mat');
%try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
prefix_base = sprintf('tmp_%s_r%d',str_xfix,nrank_sub); 
E_array_base_ = A_n_; gamma = 0.01; p_set = 0.05; n_member_lob = 2;
verbose=0;
flag_r0drop_vs_rcdrop = 0;
flag_force_create = 0;
flag_omp_use = 1;
binary_label_ = cast(zeros(M,1),'uint64');
[~,~,~,binary_label_] = ...
calllib( ...
 'halfloop_lib' ...
,'halfloop_nonbinary_f_gateway_matlab' ...
,verbose ...
,size(E_array_base_,1) ...
,size(E_array_base_,2) ...
,cast(E_array_base_,'single') ...
,flag_r0drop_vs_rcdrop ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,dir_base ...
,prefix_base ...
,flag_force_create ...
,flag_omp_use ...
,binary_label_ ...
);
label_B_1_ = label_num_to_enum_0(cast(binary_label_,'double'));
label_B_2_ = label_rearrange_0(label_B_1_,E_array_base_);
[lpv_2,lP0_2,fla_2] = label_to_label_enrichment_quad_4(label_A_,label_B_2_);
lpv_(nrank_sub) = lpv_2; lP0_(nrank_sub) = lP0_2; fla_(nrank_sub) = fla_2; label_B__{nrank_sub} = label_B_2_;
try; rmdir(sprintf('%s/dir_%s',dir_base,prefix_base),'s'); catch; disp(sprintf(' %% could not remove %s',sprintf('%s/dir_%s',dir_base,prefix_base))); end;%try;
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_hnbtZRgumb_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
%catch; disp(sprintf(' %% WARNING: error generating %s',str_hnbtZRgumb_mat)); end;%try;
delete(str_hnbtZRgumb_tmp);
end;%if (~exist(str_hnbtZRgumb,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% preliminary projection. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% spectral --> isosplit5. ;
str_xfix = test_dexcluster_multi_xfix_4('spectral_isosplit5',M,N,snr,n_cluster,n_rank,niteration);
str_spectral_isosplit5_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_spectral_isosplit5_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_spectral_isosplit5_mat,'file') & ~exist(str_spectral_isosplit5_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_spectral_isosplit5_mat)); 
save(str_spectral_isosplit5_tmp,'str_spectral_isosplit5_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
US_n_sub_ = US_n_;
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(US_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_spectral_isosplit5_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_spectral_isosplit5_mat)); end;%try;
delete(str_spectral_isosplit5_tmp);
end;%if (~exist(str_spectral_isosplit5,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne00pr --> isosplit5 ;
str_xfix = test_dexcluster_multi_xfix_4('tsne00pr_isosplit5',M,N,snr,n_cluster,n_rank,niteration);
str_tsne00pr_isosplit5_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_tsne00pr_isosplit5_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_tsne00pr_isosplit5_mat,'file') & ~exist(str_tsne00pr_isosplit5_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_tsne00pr_isosplit5_mat)); 
save(str_tsne00pr_isosplit5_tmp,'str_tsne00pr_isosplit5_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
US_n_sub_ = fast_tsne_dr_0(US_n_,struct('rand_seed',1,'no_dims',2,'theta',0));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(US_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_tsne00pr_isosplit5_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_tsne00pr_isosplit5_mat)); end;%try;
delete(str_tsne00pr_isosplit5_tmp);
end;%if (~exist(str_tsne00pr_isosplit5,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne50pr --> isosplit5 ;
str_xfix = test_dexcluster_multi_xfix_4('tsne50pr_isosplit5',M,N,snr,n_cluster,n_rank,niteration);
str_tsne50pr_isosplit5_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_tsne50pr_isosplit5_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_tsne50pr_isosplit5_mat,'file') & ~exist(str_tsne50pr_isosplit5_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_tsne50pr_isosplit5_mat)); 
save(str_tsne50pr_isosplit5_tmp,'str_tsne50pr_isosplit5_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
US_n_sub_ = fast_tsne_dr_0(US_n_,struct('rand_seed',1,'no_dims',2,'theta',0.5));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(US_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_tsne50pr_isosplit5_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_tsne50pr_isosplit5_mat)); end;%try;
delete(str_tsne50pr_isosplit5_tmp);
end;%if (~exist(str_tsne50pr_isosplit5,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00pr --> default ;
str_xfix = test_dexcluster_multi_xfix_4('umap00pr_default',M,N,snr,n_cluster,n_rank,niteration);
str_umap00pr_default_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_umap00pr_default_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_umap00pr_default_mat,'file') & ~exist(str_umap00pr_default_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_umap00pr_default_mat)); 
save(str_umap00pr_default_tmp,'str_umap00pr_default_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
[US_n_sub_, ~, label_B__{nrank_sub}]=run_umap(US_n_,'verbose','none');
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_umap00pr_default_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_umap00pr_default_mat)); end;%try;
delete(str_umap00pr_default_tmp);
end;%if (~exist(str_umap00pr_default,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00pr --> isosplit5 ;
str_xfix = test_dexcluster_multi_xfix_4('umap00pr_isosplit5',M,N,snr,n_cluster,n_rank,niteration);
str_umap00pr_isosplit5_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_umap00pr_isosplit5_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_umap00pr_isosplit5_mat,'file') & ~exist(str_umap00pr_isosplit5_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_umap00pr_isosplit5_mat)); 
save(str_umap00pr_isosplit5_tmp,'str_umap00pr_isosplit5_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
[US_n_sub_, ~, ~]=run_umap(US_n_,'verbose','none');
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(US_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_umap00pr_isosplit5_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_umap00pr_isosplit5_mat)); end;%try;
delete(str_umap00pr_isosplit5_tmp);
end;%if (~exist(str_umap00pr_isosplit5,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00pr --> hdbscan ;
str_xfix = test_dexcluster_multi_xfix_4('umap00pr_hdbscan',M,N,snr,n_cluster,n_rank,niteration);
str_umap00pr_hdbscan_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_umap00pr_hdbscan_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_umap00pr_hdbscan_mat,'file') & ~exist(str_umap00pr_hdbscan_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_umap00pr_hdbscan_mat)); 
save(str_umap00pr_hdbscan_tmp,'str_umap00pr_hdbscan_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
[US_n_sub_, ~, ~]=run_umap(US_n_,'verbose','none');
hdbscan_ = HDBSCAN(US_n_sub_); 
hdbscan_.minpts = 10; hdbscan_.fit_model(); hdbscan_.get_best_clusters(); hdbscan_.get_membership();
label_B__{nrank_sub} = hdbscan_.labels;
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_umap00pr_hdbscan_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_umap00pr_hdbscan_mat)); end;%try;
delete(str_umap00pr_hdbscan_tmp);
end;%if (~exist(str_umap00pr_hdbscan,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% louvain00pr --> default ;
str_xfix = test_dexcluster_multi_xfix_4('louvain00pr_default',M,N,snr,n_cluster,n_rank,niteration);
str_louvain00pr_default_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_louvain00pr_default_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_louvain00pr_default_mat,'file') & ~exist(str_louvain00pr_default_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_louvain00pr_default_mat)); 
save(str_louvain00pr_default_tmp,'str_louvain00pr_default_mat');
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
label_B__{nrank_sub} = GCModulMax1(US_n_*transpose(US_n_));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_louvain00pr_default_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',str_louvain00pr_default_mat)); end;%try;
delete(str_louvain00pr_default_tmp);
end;%if (~exist(str_louvain00pr_default,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%{
% dexcluster_rdrop --> ZRimax ;
str_xfix = test_dexcluster_multi_xfix_4('dexcluster_nonbinary_rdrop_trace_ZRimax',M,N,snr,n_cluster,n_rank,niteration);
str_dexcluster_nonbinary_rdrop_trace_ZRimax_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_dexcluster_nonbinary_rdrop_trace_ZRimax_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_dexcluster_nonbinary_rdrop_trace_ZRimax_mat,'file') & ~exist(str_dexcluster_nonbinary_rdrop_trace_ZRimax_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_dexcluster_nonbinary_rdrop_trace_ZRimax_mat)); 
save(str_dexcluster_nonbinary_rdrop_trace_ZRimax_tmp,'str_dexcluster_nonbinary_rdrop_trace_ZRimax_mat');
%try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
prefix_base = sprintf('tmp_%s_r%d',str_xfix,nrank_sub); 
dir_out = []; E_array_base_ = US_n_; E_array_r_ij_ = []; E_array_c_ij_ = []; gamma = 0.00; p_set = []; n_member_lob = 2;
test_loader_dexcluster_nonbinary_rdrop_trace_ZRimax_recursive_4( ...
 dir_base ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
);
ZRimax_label_ = ...
test_loader_dexcluster_nonbinary_rdrop_trace_ZRimax_recursive_4( ...
 dir_base ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
);
label_B_1_ = label_str_to_enum_0(ZRimax_label_); 
label_B_2_ = label_rearrange_0(label_B_1_,E_array_base_);
%[lpv_1,lP0_1,fla_1] = label_to_label_enrichment_quad_4(label_A_,label_B_1_);
[lpv_2,lP0_2,fla_2] = label_to_label_enrichment_quad_4(label_A_,label_B_2_);
lpv_(nrank_sub) = lpv_2; lP0_(nrank_sub) = lP0_2; fla_(nrank_sub) = fla_2; label_B__{nrank_sub} = label_B_2_;
try; rmdir(sprintf('%s/dir_%s',dir_base,prefix_base),'s'); catch; disp(sprintf(' %% could not remove %s',sprintf('%s/dir_%s',dir_base,prefix_base))); end;%try;
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_dexcluster_nonbinary_rdrop_trace_ZRimax_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
%catch; disp(sprintf(' %% WARNING: error generating %s',str_dexcluster_nonbinary_rdrop_trace_ZRimax_mat)); end;%try;
delete(str_dexcluster_nonbinary_rdrop_trace_ZRimax_tmp);
end;%if (~exist(str_dexcluster_nonbinary_rdrop_trace_ZRimax,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% dexnbrtZRgumb ;
str_xfix = test_dexcluster_multi_xfix_4('dexnbrtZRgumb',M,N,snr,n_cluster,n_rank,niteration);
str_dexnbrtZRgumb_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_dexnbrtZRgumb_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_dexnbrtZRgumb_mat,'file') & ~exist(str_dexnbrtZRgumb_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_dexnbrtZRgumb_mat)); 
save(str_dexnbrtZRgumb_tmp,'str_dexnbrtZRgumb_mat');
%try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
prefix_base = sprintf('tmp_%s_r%d',str_xfix,nrank_sub); 
dir_out = []; E_array_base_ = US_n_; E_array_r_ij_ = []; E_array_c_ij_ = []; gamma = 0.00; p_set = []; n_member_lob = 2;
dexnbrtZRgumb_recursive_5( ...
 dir_base ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
);
ZRimax_label_ = ...
dexnbrtZRgumb_recursive_5( ...
 dir_base ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
);
label_B_1_ = label_str_to_enum_0(ZRimax_label_); 
label_B_2_ = label_rearrange_0(label_B_1_,E_array_base_);
%[lpv_1,lP0_1,fla_1] = label_to_label_enrichment_quad_4(label_A_,label_B_1_);
[lpv_2,lP0_2,fla_2] = label_to_label_enrichment_quad_4(label_A_,label_B_2_);
lpv_(nrank_sub) = lpv_2; lP0_(nrank_sub) = lP0_2; fla_(nrank_sub) = fla_2; label_B__{nrank_sub} = label_B_2_;
try; rmdir(sprintf('%s/dir_%s',dir_base,prefix_base),'s'); catch; disp(sprintf(' %% could not remove %s',sprintf('%s/dir_%s',dir_base,prefix_base))); end;%try;
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_dexnbrtZRgumb_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
%catch; disp(sprintf(' %% WARNING: error generating %s',str_dexnbrtZRgumb_mat)); end;%try;
delete(str_dexnbrtZRgumb_tmp);
end;%if (~exist(str_dexnbrtZRgumb,'file')); 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% hnbrtZRgumb
str_xfix = test_dexcluster_multi_xfix_4('hnbrtZRgumb',M,N,snr,n_cluster,n_rank,niteration);
str_hnbrtZRgumb_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
str_hnbrtZRgumb_tmp = sprintf('%s/%s.tmp',dir_base,str_xfix);
if (~exist(str_hnbrtZRgumb_mat,'file') & ~exist(str_hnbrtZRgumb_tmp,'file')); 
disp(sprintf(' %% %s not found, creating',str_hnbrtZRgumb_mat)); 
save(str_hnbrtZRgumb_tmp,'str_hnbrtZRgumb_mat');
%try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
prefix_base = sprintf('tmp_%s_r%d',str_xfix,nrank_sub); 
E_array_base_ = US_n_; gamma = 0.00; p_set = 0.05; n_member_lob = 2;
verbose=0;
flag_r0drop_vs_rcdrop = 1;
flag_force_create = 0;
flag_omp_use = 1;
binary_label_ = cast(zeros(M,1),'uint64');
[~,~,~,binary_label_] = ...
calllib( ...
 'halfloop_lib' ...
,'halfloop_nonbinary_f_gateway_matlab' ...
,verbose ...
,size(E_array_base_,1) ...
,size(E_array_base_,2) ...
,cast(E_array_base_,'single') ...
,flag_r0drop_vs_rcdrop ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,dir_base ...
,prefix_base ...
,flag_force_create ...
,flag_omp_use ...
,binary_label_ ...
);
label_B_1_ = label_num_to_enum_0(cast(binary_label_,'double'));
label_B_2_ = label_rearrange_0(label_B_1_,E_array_base_);
[lpv_2,lP0_2,fla_2] = label_to_label_enrichment_quad_4(label_A_,label_B_2_);
lpv_(nrank_sub) = lpv_2; lP0_(nrank_sub) = lP0_2; fla_(nrank_sub) = fla_2; label_B__{nrank_sub} = label_B_2_;
try; rmdir(sprintf('%s/dir_%s',dir_base,prefix_base),'s'); catch; disp(sprintf(' %% could not remove %s',sprintf('%s/dir_%s',dir_base,prefix_base))); end;%try;
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(str_hnbrtZRgumb_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
%catch; disp(sprintf(' %% WARNING: error generating %s',str_hnbrtZRgumb_mat)); end;%try;
delete(str_hnbrtZRgumb_tmp);
end;%if (~exist(str_hnbrtZRgumb,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for niteration=1:n_iteration;

