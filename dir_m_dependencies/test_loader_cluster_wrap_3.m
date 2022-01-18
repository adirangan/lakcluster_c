function test_loader_cluster_wrap_3(dir_trunk,label_A_,E_rank_,I_rank_,prefix_normalization,xeta_,n_rank_xeta_sub_);
%{
  % try: ;

  %{
  %%%%%%%%;
  % cluster perforations themselves. ;
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
  E_rank_ = load(sprintf('%s/dir_mat/E_fill_.mat',dir_trunk),'E_fill_'); E_rank_ = E_rank_.E_fill_;
  I_rank_ = load(sprintf('%s/dir_mat/I_fill_.mat',dir_trunk),'I_fill_'); I_rank_ = I_rank_.I_fill_;
  %%%%%%%%;
  prefix_normalization = 'fill';
  test_loader_cluster_wrap_2(dir_trunk,label_A_,E_rank_,I_rank_,prefix_normalization,[],[]);
  %}

  %{
  %%%%%%%%;
  % cluster imputed data. ;
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
  for prefix_covariate_ = {'127'};%for prefix_covariate_ = {'127','012','RRR','UE1','UE2','UE3','UE4','UE5','UE6'};%for prefix_covariate_ = {'Phr','127','012','RRR'};
  C_rank_ = textread(sprintf('%s/dir_mat/C_rank_c%s_.tsv',dir_trunk,prefix_covariate_{1})); C_rank_ = C_rank_(:,1:end-1);
  prefix_xeta = sprintf('%s_c%s',prefix_normalization_{1},prefix_covariate_{1});
  [xeta_] = test_loader_lm_RRR_1(C_rank_,E_rank_,I_rank_,n_rank_xeta,dir_trunk,prefix_xeta);
  %%%%%%%%;
  n_rank_xeta_sub_ = [9,16];%n_rank_xeta_sub_ = [1];
  test_loader_cluster_wrap_2(dir_trunk,label_A_,E_rank_,I_rank_,prefix_normalization_{1},xeta_,n_rank_xeta_sub_);
  %%%%%%%%;
  end;%for prefix_covariate_ = {'Phr','127','012','RRR'};
  end;%for prefix_normalization_ = {'li16f','li16f_clogr','li16f_rankn'};
  %}

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
  for prefix_normalization_ = {'fill','li16f'};%for prefix_normalization_ = {'li16f','li16f_clogr','li16f_rankn'};
  E_rank_ = textread(sprintf('%s/dir_mat/E_%s_.tsv',dir_trunk,prefix_normalization_{1})); E_rank_ = E_rank_(:,1:end-1);
  I_rank_ = textread(sprintf('%s/dir_mat/I_%s_.tsv',dir_trunk,prefix_normalization_{1})); I_rank_ = I_rank_(:,1:end-1);
  %%%%%%%%;
  n_rank_xeta = 16;
  for prefix_covariate_ = {'127'};%for prefix_covariate_ = {'UE1','UE2','UE3','UE4','UE5','UE6','127','012','RRR'};
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

  %}


dir_base = sprintf('%s/dir_loader_cluster',dir_trunk);
if (~exist(dir_base,'dir')); disp(sprintf(' %% mkdir %s',dir_base)); mkdir(dir_base); end;
str_prefix = sprintf('E_%s',prefix_normalization); dir_cluster = sprintf('%s/dir_%s_cluster',dir_base,str_prefix);
if (~exist(dir_cluster,'dir')); disp(sprintf(' %% mkdir %s',dir_cluster)); mkdir(dir_cluster); end;
test_loader_cluster_0(dir_cluster,label_A_,E_rank_);
str_prefix = sprintf('I_%s',prefix_normalization); dir_cluster = sprintf('%s/dir_%s_cluster',dir_base,str_prefix);
if (~exist(dir_cluster,'dir')); disp(sprintf(' %% mkdir %s',dir_cluster)); mkdir(dir_cluster); end;
test_loader_cluster_0(dir_cluster,label_A_,I_rank_);
str_prefix = sprintf('EI_%s',prefix_normalization); dir_cluster = sprintf('%s/dir_%s_cluster',dir_base,str_prefix);
if (~exist(dir_cluster,'dir')); disp(sprintf(' %% mkdir %s',dir_cluster)); mkdir(dir_cluster); end;
test_loader_cluster_0(dir_cluster,label_A_,[E_rank_ , I_rank_]);

if (nargin>4 & ~isempty(xeta_));
if (nargin<6); n_rank_xeta_sub_=1; end;
n_u = size(E_rank_,1);
n_CCOV = size(xeta_.C_rank_,2);
n_I_GENE = xeta_.n_I_GENE;
n_E_GENE = xeta_.n_E_GENE;
n_rank_xeta_sub = length(n_rank_xeta_sub_);
for nrank_xeta_sub=1:n_rank_xeta_sub;
nrank = n_rank_xeta_sub_(nrank_xeta_sub);
disp(sprintf(' %% nrank %d <-- %d/%d',nrank,nrank_xeta_sub,n_rank_xeta_sub));
n_zeta_E = min(nrank,xeta_.n_zeta_E);
n_zeta_I = min(nrank,xeta_.n_zeta_I);
n_zeta_EI = min(nrank,xeta_.n_zeta_EI);
%%%%%%%%;
str_prefix = sprintf('E_%s_r%d',xeta_.infix,nrank); dir_cluster = sprintf('%s/dir_%s_cluster',dir_base,str_prefix);
if (~exist(dir_cluster,'dir')); disp(sprintf(' %% mkdir %s',dir_cluster)); mkdir(dir_cluster); end;
tmp_E_ = E_rank_ - (ones(n_u,1)*xeta_.zeta_E_un_(1,1:n_zeta_E) + xeta_.C_rank_(:,:)*xeta_.zeta_E_un_(2:end,1:n_zeta_E))*transpose(xeta_.zeta_E_vn_(:,1:n_zeta_E));
test_loader_cluster_0(dir_cluster,label_A_,tmp_E_);
clear tmp_E_;
%%%%%%%%;
str_prefix = sprintf('I_%s_r%d',xeta_.infix,nrank); dir_cluster = sprintf('%s/dir_%s_cluster',dir_base,str_prefix);
if (~exist(dir_cluster,'dir')); disp(sprintf(' %% mkdir %s',dir_cluster)); mkdir(dir_cluster); end;
tmp_I_ = I_rank_ - (ones(n_u,1)*xeta_.zeta_I_un_(1,1:n_zeta_I) + xeta_.C_rank_(:,:)*xeta_.zeta_I_un_(2:end,1:n_zeta_I))*transpose(xeta_.zeta_I_vn_(:,1:n_zeta_I));
test_loader_cluster_0(dir_cluster,label_A_,tmp_I_);
clear tmp_I_;
%%%%%%%%;
str_prefix = sprintf('EI_%s_r%d',xeta_.infix,nrank); dir_cluster = sprintf('%s/dir_%s_cluster',dir_base,str_prefix);
if (~exist(dir_cluster,'dir')); disp(sprintf(' %% mkdir %s',dir_cluster)); mkdir(dir_cluster); end;
tmp_EI_ = [E_rank_ , I_rank_] - (ones(n_u,1)*xeta_.zeta_EI_un_(1,1:n_zeta_EI) + xeta_.C_rank_(:,:)*xeta_.zeta_EI_un_(2:end,1:n_zeta_EI))*transpose(xeta_.zeta_EI_vn_(:,1:n_zeta_EI));
test_loader_cluster_0(dir_cluster,label_A_,tmp_EI_);
clear tmp_EI_;
%%%%%%%%;
end;%for nrank_xeta_sub=1:n_rank_xeta_sub;
end;%if (nargin>4 & ~isempty(xeta_));
