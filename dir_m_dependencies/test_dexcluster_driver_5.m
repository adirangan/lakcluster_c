function ...
test_dexcluster_driver_5( ...
 str_code ...
,A_n_ ...
,label_A_ ...
,dir_mat ...
,n_rank ...
,n_shuffle ...
,flag_rerun ...
);
% tests out dexcluster (as well as several other clustering algorithms) on data-matrix A_n_. ;

if (nargin<1);
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); str_code = []; end; na=na+1;
if (nargin<1+na); A_n_ = []; end; na=na+1;
if (nargin<1+na); label_A_ = []; end; na=na+1;
if (nargin<1+na); dir_mat = []; end; na=na+1;
if (nargin<1+na); n_rank = []; end; na=na+1;
if (nargin<1+na); n_shuffle = []; end; na=na+1;
if (nargin<1+na); flag_rerun = []; end; na=na+1;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); string_root = 'home'; end;
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jamison',string_root);

date_diff_threshold = 0.5;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;

if isempty(dir_mat); dir_mat = dir_trunk; end;
if isempty(n_rank); n_rank = rank_estimate_sample_1(A_n_); end;
if isempty(n_shuffle); n_shuffle = 64; end;
if isempty(flag_rerun); flag_rerun = 0; end;

str_n_rank = ''; if ~isempty(n_rank); str_n_rank = sprintf('_r%d',n_rank); end;
n_rank_sub = 1;

if (~exist(dir_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_mat)); mkdir(dir_mat); end;

str_xfix = sprintf('base%s',str_n_rank);
str_base = sprintf('%s/%s.mat',dir_mat,str_xfix);
if (~exist(str_base,'file')); 
disp(sprintf(' %% %s not found, creating',str_base)); 
save(str_base,'n_rank','label_A_');
end;%if (~exist(str_base,'file')); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% no preliminary projection. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne00 --> isosplit5 ;
str_xfix = 'tsne00_isosplit5';
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
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
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne50 --> isosplit5 ;
str_xfix = 'tsne50_isosplit5';
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
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
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> default ;
str_xfix = 'umap00_default';
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
[A_n_sub_, ~, label_B__{nrank_sub}]=run_umap(A_n_,'verbose','none');
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> isosplit5 ;
str_xfix = 'umap00_isosplit5';
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
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
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00 --> hdbscan ;
str_xfix = 'umap00_hdbscan';
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
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
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% louvain00 --> default;
str_xfix = 'louvain00_default';
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
label_B__{nrank_sub} = GCModulMax1(A_n_*transpose(A_n_));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% hnbtZRgumb
str_xfix = 'hnbtZRgumb';
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
%try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
prefix_base = sprintf('tmp_%s_r%d',str_xfix,nrank_sub); 
E_array_base_ = A_n_; 
parameter = struct('type','parameter');
[ ...
parameter ...
,output_label_ ...
] = ...
halfloop_nonbinary_recursive_gateway_c( ...
 parameter ...
,str_code ...
,dir_mat ...
,E_array_base_ ...
,prefix_base ...
);
label_B_1_ = label_str_to_enum_1(output_label_);
label_B_2_ = label_rearrange_0(label_B_1_,E_array_base_);
[lpv_2,lP0_2,fla_2] = label_to_label_enrichment_quad_4(label_A_,label_B_2_);
lpv_(nrank_sub) = lpv_2; lP0_(nrank_sub) = lP0_2; fla_(nrank_sub) = fla_2; label_B__{nrank_sub} = label_B_2_;
%try; rmdir(sprintf('%s/dir_%s',dir_mat,prefix_base),'s'); catch; disp(sprintf(' %% could not remove %s',sprintf('%s/dir_%s',dir_mat,prefix_base))); end;%try;
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
%catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

US_n_ = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% preliminary projection. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% spectral --> isosplit5. ;
str_xfix = sprintf('spectral_isosplit5%s',str_n_rank);
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
if isempty(US_n_); [tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank); US_n_ = tmp_U_*tmp_S_; clear tmp_U_ tmp_S_ tmp_V_; end;
US_n_sub_ = US_n_;
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(US_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne00pr --> isosplit5 ;
str_xfix = sprintf('tsne00pr_isosplit5%s',str_n_rank);
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
if isempty(US_n_); [tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank); US_n_ = tmp_U_*tmp_S_; clear tmp_U_ tmp_S_ tmp_V_; end;
US_n_sub_ = fast_tsne_dr_0(US_n_,struct('rand_seed',1,'no_dims',2,'theta',0));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(US_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% tsne50pr --> isosplit5 ;
str_xfix = sprintf('tsne50pr_isosplit5%s',str_n_rank);
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
if isempty(US_n_); [tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank); US_n_ = tmp_U_*tmp_S_; clear tmp_U_ tmp_S_ tmp_V_; end;
US_n_sub_ = fast_tsne_dr_0(US_n_,struct('rand_seed',1,'no_dims',2,'theta',0.5));
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(US_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00pr --> default ;
str_xfix = sprintf('umap00pr_default%s',str_n_rank);
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
if isempty(US_n_); [tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank); US_n_ = tmp_U_*tmp_S_; clear tmp_U_ tmp_S_ tmp_V_; end;
[US_n_sub_, ~, label_B__{nrank_sub}]=run_umap(US_n_,'verbose','none');
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00pr --> isosplit5 ;
str_xfix = sprintf('umap00pr_isosplit5%s',str_n_rank);
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
if isempty(US_n_); [tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank); US_n_ = tmp_U_*tmp_S_; clear tmp_U_ tmp_S_ tmp_V_; end;
[US_n_sub_, ~, ~]=run_umap(US_n_,'verbose','none');
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_B__{nrank_sub} = transpose(isosplit5(transpose(US_n_sub_),opts_isosplit5));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% umap00pr --> hdbscan ;
str_xfix = sprintf('umap00pr_hdbscan%s',str_n_rank);
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
if isempty(US_n_); [tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank); US_n_ = tmp_U_*tmp_S_; clear tmp_U_ tmp_S_ tmp_V_; end;
[US_n_sub_, ~, ~]=run_umap(US_n_,'verbose','none');
hdbscan_ = HDBSCAN(US_n_sub_); 
hdbscan_.minpts = 10; hdbscan_.fit_model(); hdbscan_.get_best_clusters(); hdbscan_.get_membership();
label_B__{nrank_sub} = hdbscan_.labels;
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% louvain00pr --> default ;
str_xfix = sprintf('louvain00pr_default%s',str_n_rank);
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
if isempty(US_n_); [tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank); US_n_ = tmp_U_*tmp_S_; clear tmp_U_ tmp_S_ tmp_V_; end;
label_B__{nrank_sub} = GCModulMax1(US_n_*transpose(US_n_));
[lpv_(nrank_sub),lP0_(nrank_sub),fla_(nrank_sub)] = label_to_label_enrichment_quad_4(label_A_,label_B__{nrank_sub});
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear US_n_sub_ lpv_ lP0_ fla_ label_B__ ;
catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% hnbrtZRgumb
str_xfix = sprintf('hnbrtZRgumb%s',str_n_rank);
fname_pre = sprintf('%s/%s',dir_mat,str_xfix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
%%%%
if ~flag_skip;
%try;
%%%%%%%%;
lpv_ = zeros(n_rank_sub,1); lP0_ = zeros(n_rank_sub,1); fla_ = zeros(n_rank_sub,1); 
label_B__ = cell(n_rank_sub,1);
for nrank_sub=1:n_rank_sub;
prefix_base = sprintf('tmp_%s_r%d',str_xfix,nrank_sub); 
if isempty(US_n_); [tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank); US_n_ = tmp_U_*tmp_S_; clear tmp_U_ tmp_S_ tmp_V_; end;
E_array_base_ = US_n_; gamma = 0.00; p_set = 0.05; n_member_lob = 2;
parameter = struct('type','parameter');
parameter.gamma = 0.00;
parameter.flag_r0drop_vs_rcdrop = 1;
[ ...
parameter ...
,output_label_ ...
] = ...
halfloop_nonbinary_recursive_gateway_c( ...
 parameter ...
,str_code ...
,dir_mat ...
,E_array_base_ ...
,prefix_base ...
);
label_B_1_ = label_str_to_enum_1(output_label_);
label_B_2_ = label_rearrange_0(label_B_1_,E_array_base_);
[lpv_2,lP0_2,fla_2] = label_to_label_enrichment_quad_4(label_A_,label_B_2_);
lpv_(nrank_sub) = lpv_2; lP0_(nrank_sub) = lP0_2; fla_(nrank_sub) = fla_2; label_B__{nrank_sub} = label_B_2_;
%try; rmdir(sprintf('%s/dir_%s',dir_mat,prefix_base),'s'); catch; disp(sprintf(' %% could not remove %s',sprintf('%s/dir_%s',dir_mat,prefix_base))); end;%try;
end;%for nrank_sub=1:n_rank_sub;
%%%%%%%%;
save(fname_mat,'lpv_','lP0_','fla_','label_B__');
clear A_n_sub_ lpv_ lP0_ fla_ label_B__ ;
%catch; disp(sprintf(' %% WARNING: error generating %s',fname_mat)); end;%try;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;



