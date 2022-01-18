% collect results from test_loader_cluster_wrap_0.m ;

setup;
%%%%%%%%;
dir_trunk = '/data/rangan/dir_bcc/dir_jamison';
fp_label_A_ = fopen(sprintf('%s/dir_mat/str_CLabel_sub_.nsv',dir_trunk),'r');
str_label_A_ = textscan(fp_label_A_,'%s'); fclose(fp_label_A_);
n_u = numel(str_label_A_{1});
label_A_ = label_str_to_enum_0(str_label_A_{1});

%{
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/spectral_isosplit5.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/spectral_isosplit5.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/spectral_isosplit5.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/spectral_isosplit5.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
 %}

n_rank_dexnb = 2;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/dexcluster_nonbinary_trace_ZRmax_g010_p150.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/dexcluster_nonbinary_trace_ZRmax_g010_p150.mat';
%str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/dexcluster_nonbinary_trace_ZRmax_g010_p150.mat';
%str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/dexcluster_nonbinary_trace_ZRmax_g010_p150.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
%tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
%tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
[lP0_li30p_,cap_li30p_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_.label_B__{n_rank_dexnb});
[lP0_li30p_cPhr_r3_,cap_li30p_cPhr_r3_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r3_.label_B__{n_rank_dexnb});
%[lP0_li30p_cPhr_r6_,cap_li30p_cPhr_r6_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r6_.label_B__{n_rank_dexnb});
%[lP0_li30p_cPhr_r9_,cap_li30p_cPhr_r9_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r9_.label_B__{n_rank_dexnb});

n_rank_tsne = 2;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/tsne50_isosplit5.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/tsne50_isosplit5.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/tsne50_isosplit5.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/tsne50_isosplit5.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
[lP0_li30p_,cap_li30p_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_.label_B__{n_rank_tsne});
[lP0_li30p_cPhr_r3_,cap_li30p_cPhr_r3_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r3_.label_B__{n_rank_tsne});
[lP0_li30p_cPhr_r6_,cap_li30p_cPhr_r6_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r6_.label_B__{n_rank_tsne});
[lP0_li30p_cPhr_r9_,cap_li30p_cPhr_r9_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r9_.label_B__{n_rank_tsne});

n_rank_umap = 1;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/umap00_isosplit5.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/umap00_isosplit5.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/umap00_isosplit5.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/umap00_isosplit5.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
[lP0_li30p_,cap_li30p_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_.label_B__{n_rank_umap});
[lP0_li30p_cPhr_r3_,cap_li30p_cPhr_r3_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r3_.label_B__{n_rank_umap});
[lP0_li30p_cPhr_r6_,cap_li30p_cPhr_r6_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r6_.label_B__{n_rank_umap});
[lP0_li30p_cPhr_r9_,cap_li30p_cPhr_r9_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r9_.label_B__{n_rank_umap});

n_rank_umap = 1;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/umap00_hdbscan.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/umap00_hdbscan.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/umap00_hdbscan.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/umap00_hdbscan.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
[lP0_li30p_,cap_li30p_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_.label_B__{n_rank_umap});
[lP0_li30p_cPhr_r3_,cap_li30p_cPhr_r3_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r3_.label_B__{n_rank_umap});
[lP0_li30p_cPhr_r6_,cap_li30p_cPhr_r6_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r6_.label_B__{n_rank_umap});
[lP0_li30p_cPhr_r9_,cap_li30p_cPhr_r9_] = label_to_label_enrichment_lP0(label_A_,tmp_li30p_cPhr_r9_.label_B__{n_rank_umap});





