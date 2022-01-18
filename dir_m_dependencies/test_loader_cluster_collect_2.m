% collect results from test_loader_cluster_wrap_0.m ;

setup;
%%%%%%%%;
dir_trunk = '/data/rangan/dir_bcc/dir_jamison';
fp_label_A_ = fopen(sprintf('%s/dir_mat/str_CLabel_sub_.nsv',dir_trunk),'r');
str_label_A_ = textscan(fp_label_A_,'%s'); fclose(fp_label_A_);
n_u = numel(str_label_A_{1});
prefix_normalization = 'li30p';
str_E_rank_ = sprintf('%s/dir_mat/E_%s_.tsv',dir_trunk,prefix_normalization);
fp = fopen(str_E_rank_,'r'); tmp_ = fgetl(fp); fclose(fp);
n_E_GENE = numel(str2num(tmp_)); clear tmp_;
str_I_rank_ = sprintf('%s/dir_mat/I_%s_.tsv',dir_trunk,prefix_normalization);
fp = fopen(str_I_rank_,'r'); tmp_ = fgetl(fp); fclose(fp);
n_I_GENE = numel(str2num(tmp_)); clear tmp_;
label_A_ = label_str_to_enum_0(str_label_A_{1});
%%%%%%%%;

gamma = 0.01; n_shuffle = 64; p_set = 0.05; n_member_lob = 3;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('tmp_%s',str_xfix);
dir_base = sprintf('%s/dir_loader_cluster/dir_E_li30p_cluster',dir_trunk);
dir_out = []; E_array_base_ = sparse([],[],[],n_u,n_E_GENE,0); E_array_r_ij_ = []; E_array_c_ij_ = [];
ZRimax_label_ = test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_2(dir_base,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
% Note that setting p_set=0.05 reduces the number of groups to 26 or so. ;

n_rank_dexnb = 2;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
disp(sprintf(' %% dexnb: li30p: %0.6f',tmp_li30p_.lpv_(n_rank_dexnb)));
disp(sprintf(' %% dexnb: li30p_cPhr_r3: %0.6f',tmp_li30p_cPhr_r3_.lpv_(n_rank_dexnb)));
disp(sprintf(' %% dexnb: li30p_cPhr_r6: %0.6f',tmp_li30p_cPhr_r6_.lpv_(n_rank_dexnb)));
disp(sprintf(' %% dexnb: li30p_cPhr_r9: %0.6f',tmp_li30p_cPhr_r9_.lpv_(n_rank_dexnb)));

n_rank_svd = 2;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/spectral_isosplit5.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/spectral_isosplit5.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/spectral_isosplit5.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/spectral_isosplit5.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
disp(sprintf(' %% spectral_isosplit5: li30p: %0.6f',tmp_li30p_.lpv_(n_rank_svd)));
disp(sprintf(' %% spectral_isosplit5: li30p_cPhr_r3: %0.6f',tmp_li30p_cPhr_r3_.lpv_(n_rank_svd)));
disp(sprintf(' %% spectral_isosplit5: li30p_cPhr_r6: %0.6f',tmp_li30p_cPhr_r6_.lpv_(n_rank_svd)));
disp(sprintf(' %% spectral_isosplit5: li30p_cPhr_r9: %0.6f',tmp_li30p_cPhr_r9_.lpv_(n_rank_svd)));

n_rank_tsne = 2;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/tsne50_isosplit5.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/tsne50_isosplit5.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/tsne50_isosplit5.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/tsne50_isosplit5.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
disp(sprintf(' %% tsne50_isosplit5: li30p: %0.6f',tmp_li30p_.lpv_(n_rank_tsne)));
disp(sprintf(' %% tsne50_isosplit5: li30p_cPhr_r3: %0.6f',tmp_li30p_cPhr_r3_.lpv_(n_rank_tsne)));
disp(sprintf(' %% tsne50_isosplit5: li30p_cPhr_r6: %0.6f',tmp_li30p_cPhr_r6_.lpv_(n_rank_tsne)));
disp(sprintf(' %% tsne50_isosplit5: li30p_cPhr_r9: %0.6f',tmp_li30p_cPhr_r9_.lpv_(n_rank_tsne)));

n_rank_tsne = 2;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/tsne00_isosplit5.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/tsne00_isosplit5.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/tsne00_isosplit5.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/tsne00_isosplit5.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
disp(sprintf(' %% tsne00_isosplit5: li30p: %0.6f',tmp_li30p_.lpv_(n_rank_tsne)));
disp(sprintf(' %% tsne00_isosplit5: li30p_cPhr_r3: %0.6f',tmp_li30p_cPhr_r3_.lpv_(n_rank_tsne)));
disp(sprintf(' %% tsne00_isosplit5: li30p_cPhr_r6: %0.6f',tmp_li30p_cPhr_r6_.lpv_(n_rank_tsne)));
disp(sprintf(' %% tsne00_isosplit5: li30p_cPhr_r9: %0.6f',tmp_li30p_cPhr_r9_.lpv_(n_rank_tsne)));

n_rank_umap = 1;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/umap00_isosplit5.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/umap00_isosplit5.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/umap00_isosplit5.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/umap00_isosplit5.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
disp(sprintf(' %% umap00_isosplit5: li30p: %0.6f',tmp_li30p_.lpv_(n_rank_umap)));
disp(sprintf(' %% umap00_isosplit5: li30p_cPhr_r3: %0.6f',tmp_li30p_cPhr_r3_.lpv_(n_rank_umap)));
disp(sprintf(' %% umap00_isosplit5: li30p_cPhr_r6: %0.6f',tmp_li30p_cPhr_r6_.lpv_(n_rank_umap)));
disp(sprintf(' %% umap00_isosplit5: li30p_cPhr_r9: %0.6f',tmp_li30p_cPhr_r9_.lpv_(n_rank_umap)));

n_rank_umap = 1;
str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/umap00_hdbscan.mat';
str_li30p_cPhr_r3 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r3_cluster/umap00_hdbscan.mat';
str_li30p_cPhr_r6 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r6_cluster/umap00_hdbscan.mat';
str_li30p_cPhr_r9 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r9_cluster/umap00_hdbscan.mat';
tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r3_ = load(str_li30p_cPhr_r3);
tmp_li30p_cPhr_r6_ = load(str_li30p_cPhr_r6);
tmp_li30p_cPhr_r9_ = load(str_li30p_cPhr_r9);
disp(sprintf(' %% umap00_hdbscan: li30p: %0.6f',tmp_li30p_.lpv_(n_rank_umap)));
disp(sprintf(' %% umap00_hdbscan: li30p_cPhr_r3: %0.6f',tmp_li30p_cPhr_r3_.lpv_(n_rank_umap)));
disp(sprintf(' %% umap00_hdbscan: li30p_cPhr_r6: %0.6f',tmp_li30p_cPhr_r6_.lpv_(n_rank_umap)));
disp(sprintf(' %% umap00_hdbscan: li30p_cPhr_r9: %0.6f',tmp_li30p_cPhr_r9_.lpv_(n_rank_umap)));






