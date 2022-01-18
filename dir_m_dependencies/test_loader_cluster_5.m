function...
[] = ...
test_loader_cluster_5( ...
 str_code ...
,dir_base ...
,str_label_A_ ...
,A_n_ ...
,rank_estimate_A ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);

%%%%%%%%;
na=0;
if (nargin<1+na); dir_base=pwd; end; na=na+1;
if (nargin<1+na); str_label_A_=[]; end; na=na+1;
if (nargin<1+na); A_n_=[]; end; na=na+1;
if (nargin<1+na); rank_estimate_A=[]; end; na=na+1;
if (nargin<1+na); date_diff_threshold=[]; end; na=na+1;
if (nargin<1+na); flag_force_create_mat=[]; end; na=na+1;
if (nargin<1+na); flag_force_create_tmp=[]; end; na=na+1;
%%%%;
if ( isempty(rank_estimate_A)); rank_estimate_A = []; end;
if ( isempty(date_diff_threshold)); date_diff_threshold = 1.5; end;
if ( isempty(flag_force_create_mat)); flag_force_create_mat = 0; end;
if ( isempty(flag_force_create_tmp)); flag_force_create_tmp = 0; end;
%%%%%%%%;

label_A_ = label_str_to_enum_0(str_label_A_{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% spectral --> isosplit5. ;
n_rank_svd = 6;
str_xfix = 'spectral_isosplit5';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
[tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank_svd);
lpv_ = zeros(n_rank_svd,1); lP0_ = zeros(n_rank_svd,1); fla_ = zeros(n_rank_svd,1); label_B__ = cell(n_rank_svd,1);
for nrank=1:n_rank_svd;
A_n_sub_ = tmp_U_(:,1:nrank);
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(nrank),lP0_(nrank),fla_(nrank)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank});
end;%for nrank=1:n_rank_svd;
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne00 --> isosplit5 ;
n_rank_tsne = 2;
str_xfix = 'tsne00_isosplit5';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_tsne,1); lP0_ = zeros(n_rank_tsne,1); fla_ = zeros(n_rank_tsne,1); label_B__ = cell(n_rank_tsne,1);
for nrank=1:n_rank_tsne;
A_n_sub_ = fast_tsne_dr_0(A_n_,struct('rand_seed',1,'no_dims',nrank,'theta',0));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(nrank),lP0_(nrank),fla_(nrank)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank});
end;%for nrank=1:n_rank_tsne;
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne50 --> isosplit5 ;
n_rank_tsne = 2;
str_xfix = 'tsne50_isosplit5';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_tsne,1); lP0_ = zeros(n_rank_tsne,1); fla_ = zeros(n_rank_tsne,1); 
label_B__ = cell(n_rank_tsne,1);
for nrank=1:n_rank_tsne;
A_n_sub_ = fast_tsne_dr_0(A_n_,struct('rand_seed',1,'no_dims',nrank,'theta',0.5));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(nrank),lP0_(nrank),fla_(nrank)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank});
end;%for nrank=1:n_rank_tsne;
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> default ;
str_xfix = 'umap00_default';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(1,1); lP0_ = zeros(1,1); fla_ = zeros(1,1); label_B__ = cell(1,1);
[A_n_sub_, ~, label_B__{1}]=run_umap(A_n_,'verbose','none');
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> isosplit5 ;
str_xfix = 'umap00_isosplit5';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(1,1); lP0_ = zeros(1,1); fla_ = zeros(1,1); label_B__ = cell(1,1);
[A_n_sub_, ~, ~]=run_umap(A_n_,'verbose','none');
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{1} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> hdbscan ;
str_xfix = 'umap00_hdbscan';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(1,1); lP0_ = zeros(1,1); fla_ = zeros(1,1); label_B__ = cell(1,1);
[A_n_sub_, ~, ~]=run_umap(A_n_,'verbose','none');
hdbscan_ = HDBSCAN(A_n_sub_); 
hdbscan_.minpts = 10; hdbscan_.fit_model(); hdbscan_.get_best_clusters(); hdbscan_.get_membership();
label_B__{1} = hdbscan_.labels;
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% louvain00 --> default ;
str_xfix = 'louvain00_default';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(1,1); lP0_ = zeros(1,1); fla_ = zeros(1,1); label_B__ = cell(1,1);
label_B__{1} = GCModulMax1(A_n_*transpose(A_n_));
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% dexcluster --> ZRimax ;
gamma = 0.01; n_shuffle = 64; p_set = 0.05; n_member_lob = 3;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
prefix_base = sprintf('tmp_%s',str_xfix); 
dir_out = []; E_array_base_ = A_n_; E_array_r_ij_ = []; E_array_c_ij_ = [];
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
,[] ...
,flag_force_create_mat ...
);
[ ...
 ZRimax_output_label_ ...
,ZRimax_lpFmax_label_ ...
,ZRimax_lpnext_label_ ...
] = ...
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
,[] ...
,flag_force_create_mat ...
);
lpv_ = zeros(2,1); lP0_ = zeros(2,1); fla_ = zeros(2,1); label_B__ = cell(2,1);
label_B__{1} = label_str_to_enum_0(ZRimax_output_label_); 
label_B__{2} = label_rearrange_0(label_B__{1},E_array_base_);
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
[lpv_(2),lP0_(2),fla_(2)] = label_to_label_enrichment_quad_4(label_A_,label_B__{2});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% hnbtZRgumb ;
gamma = 0.01; n_shuffle = 64; p_set = 0.05; n_member_lob = 3;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('hnbtZRgumb%s%s%s',str_g,str_p,str_nml);
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
prefix_base = sprintf('tmp_%s',str_xfix); 
dir_out = []; E_array_base_ = A_n_; E_array_r_ij_ = []; E_array_c_ij_ = [];
parameter = struct('type','parameter');
parameter.gamma = gamma;
parameter.n_shuffle = n_shuffle;
parameter.p_set = p_set;
parameter.n_member_lob = n_member_lob;
[ ...
parameter ...
,output_label_ ...
] = ...
halfloop_nonbinary_recursive_gateway_c( ...
 parameter ...
,str_code ...
,dir_base ...
,E_array_base_ ...
,prefix_base ...
);
lpv_ = zeros(2,1); lP0_ = zeros(2,1); fla_ = zeros(2,1); label_B__ = cell(2,1);
label_B__{1} = label_str_to_enum_0(output_label_); 
label_B__{2} = label_rearrange_0(label_B__{1},E_array_base_);
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
[lpv_(2),lP0_(2),fla_(2)] = label_to_label_enrichment_quad_4(label_A_,label_B__{2});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);

if ~isempty(rank_estimate_A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% preliminary projection. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
[tmp_UA_,tmp_SA_,tmp_VA_] = svds(A_n_,rank_estimate_A);
UASA_n_ = tmp_UA_*tmp_SA_;
clear tmp_UA_ tmp_SA_ tmp_VA_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% spectral --> isosplit5. ;
n_rank_svd = 1;
str_xfix = 'preproj_spectral_isosplit5';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
[tmp_U_,tmp_S_,tmp_V_] = svds(UASA_n_,n_rank_svd);
lpv_ = zeros(n_rank_svd,1); lP0_ = zeros(n_rank_svd,1); fla_ = zeros(n_rank_svd,1); label_B__ = cell(n_rank_svd,1);
for nrank=1:n_rank_svd;
A_n_sub_ = tmp_U_*tmp_S_;
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(nrank),lP0_(nrank),fla_(nrank)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank});
end;%for nrank=1:n_rank_svd;
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne00 --> isosplit5 ;
n_rank_tsne = 2;
str_xfix = 'preproj_tsne00_isosplit5';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_tsne,1); lP0_ = zeros(n_rank_tsne,1); fla_ = zeros(n_rank_tsne,1); label_B__ = cell(n_rank_tsne,1);
for nrank=1:n_rank_tsne;
A_n_sub_ = fast_tsne_dr_0(UASA_n_,struct('rand_seed',1,'no_dims',nrank,'theta',0));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(nrank),lP0_(nrank),fla_(nrank)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank});
end;%for nrank=1:n_rank_tsne;
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne50 --> isosplit5 ;
n_rank_tsne = 2;
str_xfix = 'preproj_tsne50_isosplit5';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_tsne,1); lP0_ = zeros(n_rank_tsne,1); fla_ = zeros(n_rank_tsne,1); 
label_B__ = cell(n_rank_tsne,1);
for nrank=1:n_rank_tsne;
A_n_sub_ = fast_tsne_dr_0(UASA_n_,struct('rand_seed',1,'no_dims',nrank,'theta',0.5));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(nrank),lP0_(nrank),fla_(nrank)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank});
end;%for nrank=1:n_rank_tsne;
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> default ;
str_xfix = 'preproj_umap00_default';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(1,1); lP0_ = zeros(1,1); fla_ = zeros(1,1); label_B__ = cell(1,1);
[A_n_sub_, ~, label_B__{1}]=run_umap(UASA_n_,'verbose','none');
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> isosplit5 ;
str_xfix = 'preproj_umap00_isosplit5';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(1,1); lP0_ = zeros(1,1); fla_ = zeros(1,1); label_B__ = cell(1,1);
[A_n_sub_, ~, ~]=run_umap(UASA_n_,'verbose','none');
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{1} = transpose(isosplit5(transpose(A_n_sub_),opts_isosplit5));
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> hdbscan ;
str_xfix = 'preproj_umap00_hdbscan';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(1,1); lP0_ = zeros(1,1); fla_ = zeros(1,1); label_B__ = cell(1,1);
[A_n_sub_, ~, ~]=run_umap(UASA_n_,'verbose','none');
hdbscan_ = HDBSCAN(A_n_sub_); 
hdbscan_.minpts = 10; hdbscan_.fit_model(); hdbscan_.get_best_clusters(); hdbscan_.get_membership();
label_B__{1} = hdbscan_.labels;
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% louvain00 --> default ;
str_xfix = 'preproj_louvain00_default';
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(1,1); lP0_ = zeros(1,1); fla_ = zeros(1,1); label_B__ = cell(1,1);
label_B__{1} = GCModulMax1(UASA_n_*transpose(UASA_n_));
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% dexcluster --> ZRimax ;
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 3;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('preproj_dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
prefix_base = sprintf('tmp_%s',str_xfix); 
dir_out = []; E_array_base_ = UASA_n_; E_array_r_ij_ = []; E_array_c_ij_ = [];
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
,[] ...
,flag_force_create_mat ...
);
[ ...
 ZRimax_output_label_ ...
,ZRimax_lpFmax_label_ ...
,ZRimax_lpnext_label_ ...
] = ...
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
,[] ...
,flag_force_create_mat ...
);
lpv_ = zeros(2,1); lP0_ = zeros(2,1); fla_ = zeros(2,1); label_B__ = cell(2,1);
label_B__{1} = label_str_to_enum_0(ZRimax_output_label_); 
label_B__{2} = label_rearrange_0(label_B__{1},E_array_base_);
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
[lpv_(2),lP0_(2),fla_(2)] = label_to_label_enrichment_quad_4(label_A_,label_B__{2});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% hnbtZRgumb ;
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 3;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('preproj_hnbtZRgumb%s%s%s',str_g,str_p,str_nml);
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
prefix_base = sprintf('tmp_%s',str_xfix); 
dir_out = []; E_array_base_ = UASA_n_; E_array_r_ij_ = []; E_array_c_ij_ = [];
parameter = struct('type','parameter');
parameter.gamma = gamma;
parameter.n_shuffle = n_shuffle;
parameter.p_set = p_set;
parameter.n_member_lob = n_member_lob;
parameter.flag_r0drop_vs_rcdrop = 0;
[ ...
parameter ...
,output_label_ ...
] = ...
halfloop_nonbinary_recursive_gateway_c_1( ...
 parameter ...
,str_code ...
,dir_base ...
,E_array_base_ ...
,prefix_base ...
);
lpv_ = zeros(2,1); lP0_ = zeros(2,1); fla_ = zeros(2,1); label_B__ = cell(2,1);
label_B__{1} = label_str_to_enum_0(output_label_); 
label_B__{2} = label_rearrange_0(label_B__{1},E_array_base_);
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
[lpv_(2),lP0_(2),fla_(2)] = label_to_label_enrichment_quad_4(label_A_,label_B__{2});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% hnbr0tZRgumb ;
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 3;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('preproj_hnbr0tZRgumb%s%s%s',str_g,str_p,str_nml);
tmp_fname_pre = sprintf('%s/%s',dir_base,str_xfix);
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_flag_skip;
try;
%%%%%%%%;
prefix_base = sprintf('tmp_%s',str_xfix); 
dir_out = []; E_array_base_ = UASA_n_; E_array_r_ij_ = []; E_array_c_ij_ = [];
parameter = struct('type','parameter');
parameter.gamma = gamma;
parameter.n_shuffle = n_shuffle;
parameter.p_set = p_set;
parameter.n_member_lob = n_member_lob;
parameter.flag_r0drop_vs_rcdrop = 1;
[ ...
parameter ...
,output_label_ ...
] = ...
halfloop_nonbinary_recursive_gateway_c_1( ...
 parameter ...
,str_code ...
,dir_base ...
,E_array_base_ ...
,prefix_base ...
);
lpv_ = zeros(2,1); lP0_ = zeros(2,1); fla_ = zeros(2,1); label_B__ = cell(2,1);
label_B__{1} = label_str_to_enum_0(output_label_); 
label_B__{2} = label_rearrange_0(label_B__{1},E_array_base_);
[lpv_(1),lP0_(1),fla_(1)] = label_to_label_enrichment_quad_4(label_A_,label_B__{1});
[lpv_(2),lP0_(2),fla_(2)] = label_to_label_enrichment_quad_4(label_A_,label_B__{2});
%%%%%%%%;
save(tmp_fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_fname_mat)); end;%try;
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
clear tmp_fname_pre tmp_flag_skip tmp_fname_mat;
%%%%%%%%;
test_loader_cluster_postprocess_3( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

end;%if ~isempty(rank_estimate_A);