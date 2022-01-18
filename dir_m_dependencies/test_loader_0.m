% testing a loader for table data ;
clear;
setup;
dir_trunk = sprintf('/data/rangan/dir_bcc/dir_jamison');
dir_data = sprintf('%s/data_summary_20190730',dir_trunk);
flag_load=0;
if flag_load;
str_table_name = sprintf('%s/20161026_covariate_table.format.txt',dir_data);
disp(sprintf(' %% reading %s',str_table_name));
C_ = readtable(str_table_name);
disp(sprintf(' %% saving %s',str_table_name));
save(sprintf('%s/20161026_covariate_table.format.mat',dir_data),'C_');
str_table_name = sprintf('%s/EXON.Gx_50421.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.txt',dir_data);
disp(sprintf(' %% reading %s',str_table_name));
E_ = readtable(str_table_name);
disp(sprintf(' %% writing %s',str_table_name));
save(sprintf('%s/EXON.Gx_50421.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'E_');
str_table_name = sprintf('%s/INTRON_INTERGENIC.Gx_25998.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.txt',dir_data);
disp(sprintf(' %% reading %s',str_table_name));
I_ = readtable(str_table_name);
disp(sprintf(' %% writing %s',str_table_name));
save(sprintf('%s/INTRON_INTERGENIC.Gx_25998.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'I_');
end;%if flag_load;
load(sprintf('%s/EXON.Gx_50421.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'E_');
load(sprintf('%s/INTRON_INTERGENIC.Gx_25998.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'I_');
load(sprintf('%s/20161026_covariate_table.format.mat',dir_data),'C_');

Cluster_ID_ = unique(C_.Cluster_ID_20161007);
n_Cluster_ID = length(Cluster_ID_);
n_Cluster_ID_ = zeros(n_Cluster_ID,1);
for nCluster_ID = 1:n_Cluster_ID;
n_Cluster_ID_(nCluster_ID) = length(find(strcmp(C_.Cluster_ID_20161007,Cluster_ID_(nCluster_ID))));
end;%for nCluster_ID = 1:n_Cluster_ID;
flag_plot=0;
if flag_plot;
bar(1:n_Cluster_ID,n_Cluster_ID_);
set(gca,'XTick',1:n_Cluster_ID,'XTickLabel',Cluster_ID_); xtickangle(90);
xlabel('cluster label');
ylabel('number');
title('histogram of cluster label counts');
end;%if flag_plot;

C_ID_ = C_{:,1}; 
%%%%%%%%;
E_ID_ = E_{:,1};
[~,~,C_to_E_ID_] = intersect(E_ID_,C_ID_,'stable');
C_by_E_xref_ = sparse(C_to_E_ID_,1:length(E_ID_),1,length(C_ID_),length(E_ID_));
E_by_C_xref_ = sparse(1:length(E_ID_),C_to_E_ID_,1,length(E_ID_),length(C_ID_));
flag_time=0;
if flag_time;
n_i = 1024;
% find 2 is much faster than find 1;
tic; for ni=1:n_i; find(C_by_E_xref_(:,max(1,min(length(E_ID_),floor(length(E_ID_)*rand()))))); end;%for ni=1:n_i; 
disp(sprintf(' %% C_by_E_xref_ find 2: %f',toc));
tic; for ni=1:n_i; find(C_by_E_xref_(max(1,min(length(C_ID_),floor(length(C_ID_)*rand()))),:)); end;%for ni=1:n_i; 
disp(sprintf(' %% C_by_E_xref_ find 1: %f',toc));
tic; for ni=1:n_i; find(E_by_C_xref_(:,max(1,min(length(C_ID_),floor(length(C_ID_)*rand()))))); end;%for ni=1:n_i; 
disp(sprintf(' %% E_by_C_xref_ find 2: %f',toc));
tic; for ni=1:n_i; find(E_by_C_xref_(max(1,min(length(E_ID_),floor(length(E_ID_)*rand()))),:)); end;%for ni=1:n_i; 
disp(sprintf(' %% E_by_C_xref_ find 1: %f',toc));
end;%if flag_time;
E_to_C_ID_ = zeros(length(E_ID_),1);
for nE_ID=1:length(E_ID_);
E_to_C_ID_(nE_ID) = find(C_by_E_xref_(:,nE_ID));
end;%for nE_ID=1:length(E_ID_);
flag_test = 0;
if flag_test;
disp(sprintf(' %% testing E_to_C_ID_'));
n_i = 1024;
for ni=1:n_i;
nE_ID = max(1,min(length(E_ID_),floor(length(E_ID_)*rand())));
nC_ID = E_to_C_ID_(nE_ID);
assert(strcmp(C_ID_(nC_ID),E_ID_(nE_ID)));
end;%for ni=1:n_i;
disp(sprintf(' %% finished testing E_to_C_ID_'));
end;%if flag_test;
%%%%%%%%;
I_ID_ = I_{:,1};
[~,~,C_to_I_ID_] = intersect(I_ID_,C_ID_,'stable');
C_by_I_xref_ = sparse(C_to_I_ID_,1:length(I_ID_),1,length(C_ID_),length(I_ID_));
I_by_C_xref_ = sparse(1:length(I_ID_),C_to_I_ID_,1,length(I_ID_),length(C_ID_));
flag_time=0;
if flag_time;
n_i = 1024;
% find 2 is much faster than find 1;
tic; for ni=1:n_i; find(C_by_I_xref_(:,max(1,min(length(I_ID_),floor(length(I_ID_)*rand()))))); end;%for ni=1:n_i; 
disp(sprintf(' %% C_by_I_xref_ find 2: %f',toc));
tic; for ni=1:n_i; find(C_by_I_xref_(max(1,min(length(C_ID_),floor(length(C_ID_)*rand()))),:)); end;%for ni=1:n_i; 
disp(sprintf(' %% C_by_I_xref_ find 1: %f',toc));
tic; for ni=1:n_i; find(I_by_C_xref_(:,max(1,min(length(C_ID_),floor(length(C_ID_)*rand()))))); end;%for ni=1:n_i; 
disp(sprintf(' %% I_by_C_xref_ find 2: %f',toc));
tic; for ni=1:n_i; find(I_by_C_xref_(max(1,min(length(I_ID_),floor(length(I_ID_)*rand()))),:)); end;%for ni=1:n_i; 
disp(sprintf(' %% I_by_C_xref_ find 1: %f',toc));
end;%if flag_time;
I_to_C_ID_ = zeros(length(I_ID_),1);
for nI_ID=1:length(I_ID_);
I_to_C_ID_(nI_ID) = find(C_by_I_xref_(:,nI_ID));
end;%for nI_ID=1:length(I_ID_);
flag_test = 0;
if flag_test;
disp(sprintf(' %% testing I_to_C_ID_'));
n_i = 1024;
for ni=1:n_i;
nI_ID = max(1,min(length(I_ID_),floor(length(I_ID_)*rand())));
nC_ID = I_to_C_ID_(nI_ID);
assert(strcmp(C_ID_(nC_ID),I_ID_(nI_ID)));
end;%for ni=1:n_i;
disp(sprintf(' %% finished testing I_to_C_ID_'));
end;%if flag_test;
%%%%%%%%;

%%%%%%%%;
ij_I_ID_A_ = find(strcmp(C_.Cluster_ID_20161007(I_to_C_ID_),'28')); n_A = length(ij_I_ID_A_);
ij_I_ID_B_ = find(strcmp(C_.Cluster_ID_20161007(I_to_C_ID_),'33')); n_B = length(ij_I_ID_B_);
I_A_ = decostand_total_0(I_{ij_I_ID_A_,51:end-1},'col');
I_B_ = decostand_total_0(I_{ij_I_ID_B_,51:end-1},'col');
flag_test=0;
if flag_test;
disp(sprintf(' %% testing ij_I_ID_A_'));
for nA=1:n_A;
assert(strcmp(C_ID_(I_to_C_ID_(ij_I_ID_A_(nA))),I_ID_(ij_I_ID_A_(nA))));
end;%for nA=1:n_A;
for nB=1:n_B;
assert(strcmp(C_ID_(I_to_C_ID_(ij_I_ID_B_(nB))),I_ID_(ij_I_ID_B_(nB))));
end;%for nB=1:n_B;
disp(sprintf(' %% finished testing ij_I_ID_A_'));
end;%if flag_test;
%%%%%%%%;
% determine p-values for genes alone. ;
%%%%%%%%;
n_GENE = size(I_,2)-51;
I_auc0_ = zeros(n_GENE,1);
I_beta_ = zeros(n_GENE,1);
I_praw_ = zeros(n_GENE,1);
I_pfdr_ = zeros(n_GENE,1);
tmp_label_ = [1*ones(size(I_A_,1),1);2*ones(size(I_B_,1),1)];
for ng=1:n_GENE;
if (mod(ng,10000)==0); disp(sprintf(' %% ng %d/%d',ng,n_GENE)); end;
I_auc0_(ng) = auc_0(I_A_(:,ng),I_B_(:,ng));
[tmp_b,tmp_dev,tmp_stats] = glmfit(tmp_label_,[I_A_(:,ng);I_B_(:,ng)]);
I_beta_(ng) = tmp_stats.beta(2);
I_praw_(ng) = tmp_stats.p(2);
end;%for ng=1:n_GENE;
%%%%%%%%;
% by default R ignores NAN values. ;
%%%%%%%%;
tmp_ij_ = find(isfinite(I_praw_)); I_pfdr_(tmp_ij_) = p_adjust_fdr_0(I_praw_(tmp_ij_));
%%%%%%%%;
% determine p-values for continuous covariates alone. ;
%%%%%%%%;
nctncov_ = setdiff(24:size(C_,2),28);
n_ctncov = length(nctncov_);
C_VariableName_ = cell(n_ctncov,1);
C_auc0_ = zeros(n_ctncov,1);
C_beta_ = zeros(n_ctncov,1);
C_praw_ = zeros(n_ctncov,1);
C_pfdr_ = zeros(n_ctncov,1);
tmp_label_ = [1*ones(size(I_A_,1),1);2*ones(size(I_B_,1),1)];
C_A_ = zeros(n_A,n_ctncov);
C_B_ = zeros(n_B,n_ctncov);
for nc=1:n_ctncov;
nctncov = nctncov_(nc);
tmp_C_A_ = C_(I_to_C_ID_(ij_I_ID_A_),nctncov); 
disp(sprintf(' %% nc %d nctncov %d class %s',nc,nctncov,class(tmp_C_A_{1,1})));
end;%for nc=1:n_ctncov;
for nc=1:n_ctncov;
nctncov = nctncov_(nc);
C_VariableName_{nc} = C_.Properties.VariableNames{nctncov};
tmp_C_A_ = C_(I_to_C_ID_(ij_I_ID_A_),nctncov); 
if strcmp(class(tmp_C_A_{1,1}),'double'); C_A_(:,nc) = tmp_C_A_{1:n_A,1};
elseif strcmp(class(tmp_C_A_{1,1}),'cell'); C_A_(:,nc) = cellfun(@str2num,tmp_C_A_{1:n_A,1}); 
else disp(sprintf(' %% nc %d nctncov %d class %s',nc,nctncov,class(tmp_C_A_{1,1}))); end;
tmp_C_B_ = C_(I_to_C_ID_(ij_I_ID_B_),nctncov); 
if strcmp(class(tmp_C_B_{1,1}),'double'); C_B_(:,nc) = tmp_C_B_{1:n_B,1};
elseif strcmp(class(tmp_C_B_{1,1}),'cell'); C_B_(:,nc) = cellfun(@str2num,tmp_C_B_{1:n_B,1});
else disp(sprintf(' %% nc %d nctncov %d class %s',nc,nctncov,class(tmp_C_B_{1,1}))); end;
C_auc0_(nc) = auc_0(C_A_(:,nc),C_B_(:,nc));
[tmp_b,tmp_dev,tmp_stats] = glmfit(tmp_label_,[C_A_(:,nc);C_B_(:,nc)]);
C_beta_(nc) = tmp_stats.beta(2);
C_praw_(nc) = tmp_stats.p(2);
end;%for nc=1:n_ctncov;
%%%%%%%%;
% by default R ignores NAN values. ;
%%%%%%%%;
tmp_ij_ = find(isfinite(C_praw_)); C_pfdr_(tmp_ij_) = p_adjust_fdr_0(C_praw_(tmp_ij_));

%%%%%%%%;
% Now sort genes by p-value and plot heatmap: ;
%%%%%%%%;
[~,tmp_I_ij_] = sort(I_auc0_);
tmp_AB_ = [I_A_(:,tmp_I_ij_) ; I_B_(:,tmp_I_ij_)];
imagedisp_0(decostand_total_0(tmp_AB_,'row'));
%%%%%%%%;
% Now sort covariates by p-value and plot heatmap: ;
%%%%%%%%;
[~,tmp_C_ij_] = sort(C_auc0_);
tmp_AB_ = [C_A_(:,tmp_C_ij_) ; C_B_(:,tmp_C_ij_)];
imagedisp_0(decostand_total_0(tmp_AB_,'row'));

%%%%%%%%;
% Now set up a simple linear model to assess gene + covariate effects. ;
%%%%%%%%;
tmp_An_ = [ones(n_A,1) , I_A_ , C_A_ ; ones(n_B,1) , I_B_ , C_B_];
[xn_,nb] = conjgrad_0(tmp_An_,tmp_label_,zeros(1+n_GENE+n_ctncov,1));
xn_1_ = xn_(1);
xn_I_ = xn_(1 + (1:n_GENE));
xn_C_ = xn_(1 + n_GENE + (1:n_ctncov));

%%%%%%%%;
% what does linear model do if both genes and covariates explain labels? ;
%%%%%%%%;
tmp_mu = 10.5;
tmp_An_ =  [ones(n_A,1) , randn(n_A,n_GENE) - tmp_mu , randn(n_A,n_ctncov) - tmp_mu ; ...
            ones(n_B,1) , randn(n_B,n_GENE) + tmp_mu , randn(n_B,n_ctncov) + tmp_mu ];
[xn_,nb] = conjgrad_0(tmp_An_,tmp_label_,zeros(1+n_GENE+n_ctncov,1));
xn_1_ = xn_(1);
xn_I_ = xn_(1 + (1:n_GENE));
xn_C_ = xn_(1 + n_GENE + (1:n_ctncov));
disp(sprintf(' %% mean(xn_I_) %0.8f std(xn_I_) %0.8f mean(xn_C_) %0.8f std(xn_C_) %0.8f',mean(xn_I_),std(xn_I_),mean(xn_C_),std(xn_C_)));

%%%%%%%%;
% try an even simpler linear model. ;
%%%%%%%%;
tmp_mu = +1.5; n_test = 500;
tmp_label_balanced_ = tmp_label_ - mean(tmp_label_);
tmp_An_ =  [ones(n_A,1) , 1*randn(n_A,n_test) - tmp_mu ; ...
            ones(n_B,1) , 1*randn(n_B,n_test) + tmp_mu ];
[xn_,nb] = conjgrad_0(tmp_An_,tmp_label_balanced_,zeros(1+n_test,1));
xn_1_ = xn_(1);
xn_T_ = xn_(1 + (1:n_test));
disp(sprintf(' %% mean(xn_T_) %0.8f std(xn_T_) %0.8f ',mean(xn_T_),std(xn_T_)));
prctile(xn_T_,0:10:100),;


%{
tmp_m = 23; %<-- colnames of covariate matrix. ;
tmp_VariableName = C_.Properties.VariableNames{tmp_m}; %<-- Should be 'RunID'. ;
C_A_ = C_(I_to_C_ID_(ij_I_ID_A_),tmp_m); C_A_ = C_A_{1:n_A,1};
C_B_ = C_(I_to_C_ID_(ij_I_ID_B_),tmp_m); C_B_ = C_B_{1:n_B,1};
tmp_C_ = [C_A_ ; C_B_];
ng = 3;
[tmp_b,tmp_dev,tmp_stats] = glmfit([tmp_label_ , tmp_C_],[I_A_(:,ng);I_B_(:,ng)]);
% not sure what to do here; what does anova(lm) do in R?. ;
 %}

%{
tmp_m = 24; %<-- colnames of covariate matrix. ;
tmp_VariableName = C_.Properties.VariableNames{tmp_m}; %<-- Should be 'cDNA_Pico_AdjConc'. ;
C_A_ = C_(I_to_C_ID_(ij_I_ID_A_),tmp_m); C_A_ = cellfun(@str2num,C_A_{1:n_A,1});
C_B_ = C_(I_to_C_ID_(ij_I_ID_B_),tmp_m); C_B_ = cellfun(@str2num,C_B_{1:n_B,1});
tmp_C_ = [C_A_ ; C_B_];
for ng=1:40;
[tmp_b,tmp_dev,tmp_stats] = glmfit([tmp_label_ , tmp_C_],[I_A_(:,ng);I_B_(:,ng)]);
disp(sprintf(' %% pval: %f %f %f',tmp_stats.p));
tmp_11_(ng) = tmp_stats.p(2);
tmp_12_(ng) = tmp_stats.p(3);
end;%for ng=1:40;
 %}

%{
tmp_20_ = [1.158716e-02,  3.386000e-01,  9.660253e-01,  1.220678e-01 ... 
, 5.211601e-02,  9.403348e-06,  9.575741e-02,  3.789066e-02 ... 
, 1.244929e-03,  2.347152e-03,  3.358386e-02,  8.429909e-04 ... 
, 7.571563e-01,  2.081093e-02,  1.051360e-03,  2.883896e-01 ... 
, 1.623833e-01,  5.214384e-01,  3.948735e-01,  9.419286e-08 ... 
, 9.156190e-02,  4.645287e-02,  5.817189e-05,  2.579568e-05 ... 
, 9.660364e-02,  6.376687e-03,  1.019442e-01,  3.081057e-05 ... 
, 9.882653e-02,  2.322626e-03,  1.438826e-05,  8.092762e-02 ... 
, 9.276010e-02,  5.780512e-02,  4.947630e-03,  3.592598e-04 ... 
, 7.417675e-02,  1.126596e-03,  7.337254e-01,  1.245327e-01 ];... 
  %}


