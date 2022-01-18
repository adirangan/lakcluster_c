  %%%%%%%%;
  % cluster imputed data for UEx. ;
  %%%%%%%%;
  clear;
  setup;
  %%%%%%%%;
  dir_trunk = '/data/rangan/dir_bcc/dir_jamison';
  fp_label_A_ = fopen(sprintf('%s/dir_mat/str_CLabel_sub_.nsv',dir_trunk),'r');
  str_label_A_ = textscan(fp_label_A_,'%s'); fclose(fp_label_A_);
  n_u = numel(str_label_A_{1});
  label_A_ = label_str_to_enum_0(str_label_A_{1});
  %%%%%%%%;
  for prefix_normalization_ = {'li16f'};%for prefix_normalization_ = {'li16f','li16f_clogr','li16f_rankn'};
  E_rank_ = textread(sprintf('%s/dir_mat/E_%s_.tsv',dir_trunk,prefix_normalization_{1})); E_rank_ = E_rank_(:,1:end-1);
  I_rank_ = textread(sprintf('%s/dir_mat/I_%s_.tsv',dir_trunk,prefix_normalization_{1})); I_rank_ = I_rank_(:,1:end-1);
  %%%%%%%%;
  n_rank_xeta = 16;
  for prefix_covariate_ = {'UE1','UE2','UE3','UE4','UE5','UE6','127','012','RRR'};
  C_rank_ = textread(sprintf('%s/dir_mat/C_rank_c%s_.tsv',dir_trunk,prefix_covariate_{1})); C_rank_ = C_rank_(:,1:end-1);
  prefix_xeta = sprintf('%s_c%s',prefix_normalization_{1},prefix_covariate_{1});
  [xeta_] = test_loader_lm_RRR_1(C_rank_,E_rank_,I_rank_,n_rank_xeta,dir_trunk,prefix_xeta);
  %%%%%%%%;
  if (strcmp(prefix_covariate_{1},'UE1')); n_rank_xeta_sub_ = [1]; end;
  if (strcmp(prefix_covariate_{1},'UE2')); n_rank_xeta_sub_ = [1,5]; end;
  if (strcmp(prefix_covariate_{1},'UE3')); n_rank_xeta_sub_ = [1,8]; end;
  if (strcmp(prefix_covariate_{1},'UE4')); n_rank_xeta_sub_ = [1]; end;
  if (strcmp(prefix_covariate_{1},'UE5')); n_rank_xeta_sub_ = [1]; end;
  if (strcmp(prefix_covariate_{1},'UE6')); n_rank_xeta_sub_ = [1,2,8]; end;
  if (strcmp(prefix_covariate_{1},'127')); n_rank_xeta_sub_ = [12]; end;
  if (strcmp(prefix_covariate_{1},'012')); n_rank_xeta_sub_ = [12]; end;
  if (strcmp(prefix_covariate_{1},'RRR')); n_rank_xeta_sub_ = [12]; end;
  test_loader_cluster_wrap_2(dir_trunk,label_A_,E_rank_,I_rank_,prefix_normalization_{1},xeta_,n_rank_xeta_sub_);
  %%%%%%%%;
  end;%for prefix_covariate_ = {'Phr','127','012','RRR'};
  end;%for prefix_normalization_ = {'li16f','li16f_clogr','li16f_rankn'};
