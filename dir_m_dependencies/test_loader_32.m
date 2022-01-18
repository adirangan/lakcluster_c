% testing a loader for table data ;
% running loader_29 again, except with a different imputation seed. ;
% 20210110: checking the original (perforated) data for rank. ;
% 20210110: checking the missingness for structure. ;
% 20210115: switching from RRR to PCA for covariate-correction. ;

clear;

platform = 'access1';
%platform = 'OptiPlex';
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;

dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jamison',string_root);
dir_data = sprintf('%s/data_summary_20190730',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_cluster = sprintf('%s/dir_loader_cluster',dir_trunk);

%%%%%%%%;
flag_load=1;
%%%%%%%%;
if ~flag_load;
%%%%%%%%;
str_table_name = sprintf('%s/20161026_covariate_table.format.txt',dir_data);
disp(sprintf(' %% reading %s',str_table_name));
C_ = readtable(str_table_name);
disp(sprintf(' %% saving %s',str_table_name));
save(sprintf('%s/20161026_covariate_table.format.mat',dir_data),'C_');
%%%%%%%%;
str_table_name = sprintf('%s/EXON.Gx_50421.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.txt',dir_data);
disp(sprintf(' %% reading %s',str_table_name));
E_ = readtable(str_table_name);
disp(sprintf(' %% writing %s',str_table_name));
save(sprintf('%s/EXON.Gx_50421.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'E_');
%%%%%%%%;
str_table_name = sprintf('%s/INTRON_INTERGENIC.Gx_25998.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.txt',dir_data);
disp(sprintf(' %% reading %s',str_table_name));
I_ = readtable(str_table_name);
disp(sprintf(' %% writing %s',str_table_name));
save(sprintf('%s/INTRON_INTERGENIC.Gx_25998.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'I_');
%%%%%%%%;
end;%if ~flag_load;
%%%%%%%%;
if flag_load;
%%%%%%%%;
load(sprintf('%s/EXON.Gx_50421.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'E_');
load(sprintf('%s/INTRON_INTERGENIC.Gx_25998.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'I_');
load(sprintf('%s/20161026_covariate_table.format.mat',dir_data),'C_');
%%%%%%%%;
end;%if flag_load;

flag_E = 1;
flag_I = 1;

%%%%%%%%;
% Here we attempt to extract the gene-lengths. ;
%%%%%%%%;
str_table_name = sprintf('%s/EXON.Gx_50421.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.txt',dir_data);
fp = fopen(str_table_name); tmp_str = fgetl(fp); fclose(fp);
E_VariableName_ = strread(tmp_str,'%s','delimiter','\t');
clear tmp_str ;
%%%%%%%%;
str_table_name = sprintf('%s/INTRON_INTERGENIC.Gx_25998.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.txt',dir_data);
fp = fopen(str_table_name); tmp_str = fgetl(fp); fclose(fp);
I_VariableName_ = strread(tmp_str,'%s','delimiter','\t');
clear tmp_str ;
%%%%%%%%;
str_table_name = sprintf('%s/ensg_lengths.tsv',dir_data);
fp = fopen(str_table_name); gene_length_ = textscan(fp,'%s %d'); fclose(fp);
for ng=1:length(gene_length_{1});
end;%for ng=1:length(gene_length_{1});
%%%%%%%%;
[tmp_cap_,E_VariableName_to_cap_,gene_length_to_cap_] = intersect(E_VariableName_,gene_length_{1},'stable');
E_length_ = zeros(length(E_VariableName_),1);
E_length_(E_VariableName_to_cap_) = gene_length_{2}(gene_length_to_cap_);
flag_test=1;
if flag_test;
disp(sprintf(' %% testing E_length_ vs gene_length_ '));
n_i=1024;
for ni=1:n_i;
ng = max(1,min(length(tmp_cap_),floor(length(tmp_cap_)*rand())));
assert(strcmp(E_VariableName_(E_VariableName_to_cap_(ng)),gene_length_{1}(gene_length_to_cap_(ng))));
end;%for ni=1:n_i;
end;%if flag_test;
clear tmp_cap_;
%%%%%%%%;
[tmp_cap_,I_VariableName_to_cap_,gene_length_to_cap_] = intersect(I_VariableName_,gene_length_{1},'stable');
I_length_ = zeros(length(I_VariableName_),1);
I_length_(I_VariableName_to_cap_) = gene_length_{2}(gene_length_to_cap_);
flag_test=1;
if flag_test;
disp(sprintf(' %% testing I_length_ vs gene_length_ '));
n_i=1024;
for ni=1:n_i;
ng = max(1,min(length(tmp_cap_),floor(length(tmp_cap_)*rand())));
assert(strcmp(I_VariableName_(I_VariableName_to_cap_(ng)),gene_length_{1}(gene_length_to_cap_(ng))));
end;%for ni=1:n_i;
end;%if flag_test;
clear tmp_cap_;
%%%%%%%%;
% Note that only 6660 (13%) of the E_ are found, and only 5379 (21%) of the I_ are found. ;
%%%%%%%%;

Label_ID_ = unique(C_.Cluster_ID_20161007);
n_Label_ID = length(Label_ID_);
n_Label_ID_ = zeros(n_Label_ID,1);
for nLabel_ID = 1:n_Label_ID;
n_Label_ID_(nLabel_ID) = length(find(strcmp(C_.Cluster_ID_20161007,Label_ID_(nLabel_ID))));
end;%for nLabel_ID = 1:n_Label_ID;
flag_plot=0;
if flag_plot;
bar(1:n_Label_ID,n_Label_ID_);
set(gca,'XTick',1:n_Label_ID,'XTickLabel',Label_ID_); xtickangle(90);
xlabel('cluster label');
ylabel('number');
title('histogram of cluster label counts');
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_base = sprintf('%s/dir_jpg/label_count_C',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% First pass at defining E_val_ and I_val_ ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
E_ID_ = E_{:,1};
I_ID_ = E_{:,1};
C_ID_ = C_{:,1}; 
u_ID_ = intersectall({E_ID_,I_ID_,C_ID_});
n_u = length(u_ID_);
[~,~,E_to_u_ID_] = intersect(u_ID_,E_ID_,'stable');
[~,~,I_to_u_ID_] = intersect(u_ID_,I_ID_,'stable');
[~,~,C_to_u_ID_] = intersect(u_ID_,C_ID_,'stable');
%%%%%%%%;
flag_test=1;
if flag_test;
disp(sprintf(' %% testing E_to_C_ID_'));
n_i = 1024;
for ni=1:n_i;
nu = max(1,min(n_u,floor(n_u*rand())));
nE_ID = E_to_u_ID_(nu);
nI_ID = I_to_u_ID_(nu);
nC_ID = C_to_u_ID_(nu);
assert(strcmp(C_ID_(nC_ID),E_ID_(nE_ID)));
assert(strcmp(C_ID_(nC_ID),I_ID_(nE_ID)));
end;%for ni=1:n_i;
disp(sprintf(' %% finished testing E_to_C_ID_'));
end;%if flag_test;
%%%%%%%%;
% Extracting unique cluster labels. ;
%%%%%%%%;
u_Sample_Label_ = unique(C_.Cluster_ID_20161007(C_to_u_ID_));
n_Sample_Label = length(u_Sample_Label_);
n_Sample_Label_ = zeros(n_Sample_Label,1);
for nSample_Label = 1:n_Sample_Label;
n_Sample_Label_(nSample_Label) = length(find(strcmp(C_.Cluster_ID_20161007(C_to_u_ID_),u_Sample_Label_(nSample_Label))));
end;%for nSample_Label = 1:n_Sample_Label;
flag_plot=0;
if flag_plot;
bar(1:n_Sample_Label,n_Sample_Label_);
set(gca,'XTick',1:n_Sample_Label,'XTickLabel',u_Sample_Label_); xtickangle(90);
xlabel('cluster label');
ylabel('number');
title('histogram of cluster label counts');
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_base = sprintf('%s/dir_jpg/label_count_I',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;
%%%%%%%%;
% Generating ordered matrices of gene-values. ;
%%%%%%%%;
if flag_E;
E_col_val_ = 52:size(E_,2)-1; %<-- gene list begins after 'Fail Confidence Sum' ;
%E_val_ = decostand_total_0(E_{E_to_u_ID_,E_col_val_},'col');
E_val_ = E_{E_to_u_ID_,E_col_val_};
n_E_GENE = size(E_val_,2);
end;%if flag_E;
if flag_I;
I_col_val_ = 52:size(I_,2)-1;
%I_val_ = decostand_total_0(I_{I_to_u_ID_,I_col_val_},'col');
I_val_ = I_{I_to_u_ID_,I_col_val_};
n_I_GENE = size(I_val_,2);  %<-- gene list begins after 'Fail Confidence Sum' ;
end;%if flag_I;
%assert(size(I_val_,1)==size(E_val_,1)); %<-- fewer samples in E_. ;
%%%%%%%%;
% Now check for bimodality of gene expression across samples: ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
clf;
subplot(2,1,1); tmp_b_ = 0:0.25:6; tmp_h_ = hist(log10(1+E_{E_to_u_ID_,E_col_val_}),tmp_b_); colormap(colormap_beach()); imagesc(log10(1+tmp_h_),[0,3]); xlabel('gene index'); ylabel('log #'); title('E distribution');
subplot(2,1,2); tmp_b_ = 0:0.25:6; tmp_h_ = hist(log10(1+I_{I_to_u_ID_,I_col_val_}),tmp_b_); colormap(colormap_beach()); imagesc(log10(1+tmp_h_),[0,3]); xlabel('gene index'); ylabel('log #'); title('I distribution');
fname_base = sprintf('%s/dir_jpg/expression_histogram',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% Determining cutoffs for outliers. ;
%%%%%%%%;
E_val_m1_ = sum(E_val_<=0,1)./n_u; 
I_val_m1_ = sum(I_val_<=0,1)./n_u; 
m1_threshold = 0.90; %<-- 90 percent missing. ;
disp(sprintf(' %% s1 threshold %0.2f',m1_threshold));
clf;
subplot(1,2,1);plot(1:n_E_GENE,sort(E_val_m1_),'.',1:n_E_GENE,m1_threshold*ones(1,n_E_GENE),'k:'); xlabel('sorted gene index'); ylabel('fraction missing samples'); xlim([1,n_E_GENE]); ylim([0,1]);
subplot(1,2,2);plot(1:n_I_GENE,sort(I_val_m1_),'.',1:n_I_GENE,m1_threshold*ones(1,n_I_GENE),'k:'); xlabel('sorted gene index'); ylabel('fraction missing samples'); xlim([1,n_I_GENE]); ylim([0,1]);
flag_plot=0;
if flag_plot;
fname_base = sprintf('%s/dir_jpg/expression_m1_cdf',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;
E_m1_ij_ = find(E_val_m1_<m1_threshold); I_m1_ij_ = find(I_val_m1_<m1_threshold);
%E_m1_ij_ = 1:1e4; I_m1_ij_ = 1:1e4;
E_val_m2_ = sum(log10(1e-12 + E_val_(:,E_m1_ij_)),2);
I_val_m2_ = sum(log10(1e-12 + I_val_(:,I_m1_ij_)),2);
clf;
m2_threshold = so2g_mle_fminsearch(E_val_m2_ + I_val_m2_);
flag_plot=0;
if flag_plot;
fname_base = sprintf('%s/dir_jpg/expression_m2_histogram',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;
disp(sprintf(' %% m2 threshold %0.2f',m2_threshold));
X_m2_ij_ = find(E_val_m2_+I_val_m2_ > m2_threshold);
disp(sprintf(' %% %d/%d = %0.2f samples above m2_threshold',length(X_m2_ij_),n_u,length(X_m2_ij_)/n_u));
X_m2_ij_ = 1:n_u; disp(sprintf(' %% retaining all samples. '));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Second pass at defining E_val_ and I_val_ ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
u_ID_ = u_ID_(X_m2_ij_);
n_u = length(u_ID_);
[~,~,E_to_u_ID_] = intersect(u_ID_,E_ID_,'stable');
[~,~,I_to_u_ID_] = intersect(u_ID_,I_ID_,'stable');
[~,~,C_to_u_ID_] = intersect(u_ID_,C_ID_,'stable');
%%%%%%%%;
flag_test=1;
if flag_test;
disp(sprintf(' %% testing E_to_C_ID_'));
n_i = 1024;
for ni=1:n_i;
nu = max(1,min(n_u,floor(n_u*rand())));
nE_ID = E_to_u_ID_(nu);
nI_ID = I_to_u_ID_(nu);
nC_ID = C_to_u_ID_(nu);
assert(strcmp(C_ID_(nC_ID),E_ID_(nE_ID)));
assert(strcmp(C_ID_(nC_ID),I_ID_(nE_ID)));
end;%for ni=1:n_i;
disp(sprintf(' %% finished testing E_to_C_ID_'));
end;%if flag_test;
%%%%%%%%;
% Extracting unique cluster labels. ;
%%%%%%%%;
str_Sample_Label_ = C_.Cluster_ID_20161007(C_to_u_ID_);
u_Sample_Label_ = unique(C_.Cluster_ID_20161007(C_to_u_ID_));
n_Sample_Label = length(u_Sample_Label_);
n_Sample_Label_ = zeros(n_Sample_Label,1);
for nSample_Label = 1:n_Sample_Label;
n_Sample_Label_(nSample_Label) = length(find(strcmp(C_.Cluster_ID_20161007(C_to_u_ID_),u_Sample_Label_(nSample_Label))));
end;%for nSample_Label = 1:n_Sample_Label;
Sample_Label_sub_ = zeros(n_u,1);
for nSample_Label=1:n_Sample_Label;
tmp_ij_ = find(strcmp(str_Sample_Label_,u_Sample_Label_(nSample_Label)));
if (~isempty(str2num(u_Sample_Label_{nSample_Label})));
Sample_Label_sub_(tmp_ij_) = str2num(u_Sample_Label_{nSample_Label});
end;%if (~isempty(str2num(u_Sample_Label_{nSample_Label})));
end;%for nSample_Label=1:n_Sample_Label;
flag_plot=0;
if flag_plot;
bar(1:n_Sample_Label,n_Sample_Label_);
set(gca,'XTick',1:n_Sample_Label,'XTickLabel',u_Sample_Label_); xtickangle(90);
xlabel('cluster label');
ylabel('number');
title('histogram of cluster label counts');
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_base = sprintf('%s/dir_jpg/label_count_sub',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;
%%%%%%%%;
% Generating ordered matrices of gene-values. ;
%%%%%%%%;
if flag_E;
E_col_val_ = 52:size(E_,2)-1; %<-- gene list begins after 'Fail Confidence Sum' ;
E_col_val_ = E_col_val_(E_m1_ij_);
E_VariableName_sub_ = E_VariableName_(E_col_val_);
%E_val_ = decostand_total_0(E_{E_to_u_ID_,E_col_val_},'col');
E_val_ = E_{E_to_u_ID_,E_col_val_};
n_E_GENE = size(E_val_,2);
end;%if flag_E;
if flag_I;
I_col_val_ = 52:size(I_,2)-1;
I_col_val_ = I_col_val_(I_m1_ij_);
I_VariableName_sub_ = I_VariableName_(I_col_val_);
%I_val_ = decostand_total_0(I_{I_to_u_ID_,I_col_val_},'col');
I_val_ = I_{I_to_u_ID_,I_col_val_};
n_I_GENE = size(I_val_,2);  %<-- gene list begins after 'Fail Confidence Sum' ;
end;%if flag_I;
assert(size(I_val_,1)==size(E_val_,1)); %<-- same number of samples in I_ and E_. ;
%%%%%%%%;
% Now check for bimodality of gene expression across samples: ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
clf;
subplot(2,1,1); tmp_b_ = 0:0.25:6; tmp_h_ = hist(log10(1+E_{E_to_u_ID_,E_col_val_}),tmp_b_); colormap(colormap_beach()); imagesc(log10(1+tmp_h_),[0,3]); xlabel('gene index'); ylabel('log #'); title('E distribution');
subplot(2,1,2); tmp_b_ = 0:0.25:6; tmp_h_ = hist(log10(1+I_{I_to_u_ID_,I_col_val_}),tmp_b_); colormap(colormap_beach()); imagesc(log10(1+tmp_h_),[0,3]); xlabel('gene index'); ylabel('log #'); title('I distribution');
fname_base = sprintf('%s/dir_jpg/expression_histogram_sub',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;
%%%%%%%%;

%{
%%%%%%%%;
% look at bimodal distribution of values. ;
% Does not appear as though the zeros are ordered. ;
%%%%%%%%;
tmp_ij_ = find(E_val_~=0);
E_nz_ = log(E_val_(tmp_ij_));
[E_lval_threshold] = so2g_mle_fminsearch(E_nz_);
tmp_E_ = zeros(size(E_val_));
tmp_E_(tmp_ij_(find(E_nz_<E_lval_threshold)))=-1;
tmp_E_(tmp_ij_(find(E_nz_>E_lval_threshold)))=+1;
tmp_E_(find(E_val_==0))=0;
flag_plot=0;
if flag_plot;
subplot(4,3,[1,2,3,4,5,6]); imagesc(tmp_E_,[-1,+1]); colormap(colormap_beach()); title('E_lval thresholded');
subplot(4,3,7); plot(sum(tmp_E_==-1,2),'.'); title('m2 -1'); xlim([1,n_u]);
subplot(4,3,8); plot(sum(tmp_E_== 0,2),'.'); title('m2  0'); xlim([1,n_u]);
subplot(4,3,9); plot(sum(tmp_E_==+1,2),'.'); title('m2 +1'); xlim([1,n_u]);
subplot(4,3,10); plot(sum(tmp_E_==-1,1),'.'); title('m1 -1'); xlim([1,n_E_GENE]);
subplot(4,3,11); plot(sum(tmp_E_== 0,1),'.'); title('m1  0'); xlim([1,n_E_GENE]);
subplot(4,3,12); plot(sum(tmp_E_==+1,1),'.'); title('m1 +1'); xlim([1,n_E_GENE]);
figbig;
fname_base = sprintf('%s/dir_jpg/E_lval_so2g',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;
%%%%%%%%;
tmp_ij_ = find(I_val_~=0);
I_nz_ = log(I_val_(tmp_ij_));
[I_lval_threshold] = so2g_mle_fminsearch(I_nz_);
tmp_I_ = zeros(size(I_val_));
tmp_I_(tmp_ij_(find(I_nz_<I_lval_threshold)))=-1;
tmp_I_(tmp_ij_(find(I_nz_>I_lval_threshold)))=+1;
tmp_I_(find(I_val_==0))=0;
flag_plot=0;
if flag_plot;
subplot(4,3,[1,2,3,4,5,6]); imagesc(tmp_I_,[-1,+1]); colormap(colormap_beach()); title('I_lval thresholded');
subplot(4,3,7); plot(sum(tmp_I_==-1,2),'.'); title('m2 -1'); xlim([1,n_u]);
subplot(4,3,8); plot(sum(tmp_I_== 0,2),'.'); title('m2  0'); xlim([1,n_u]);
subplot(4,3,9); plot(sum(tmp_I_==+1,2),'.'); title('m2 +1'); xlim([1,n_u]);
subplot(4,3,10); plot(sum(tmp_I_==-1,1),'.'); title('m1 -1'); xlim([1,n_I_GENE]);
subplot(4,3,11); plot(sum(tmp_I_== 0,1),'.'); title('m1  0'); xlim([1,n_I_GENE]);
subplot(4,3,12); plot(sum(tmp_I_==+1,1),'.'); title('m1 +1'); xlim([1,n_I_GENE]);
figbig;
fname_base = sprintf('%s/dir_jpg/I_lval_so2g',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;
%%%%%%%%;
 %}

%%%%%%%%;
% determine initial rank of E_val_. ;
%%%%%%%%;
flag_recalc = 0;
fname_mat = sprintf('%s/E_svd_ori_.mat',dir_mat);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
E_index_missing_ = efind(E_val_<=0); E_index_filled_ = setdiff(0:n_u*n_E_GENE-1,E_index_missing_);
%hist(log(max(1e-12,E_val_(1+E_index_filled_))),128);
n_shuffle = 32; n_svd = 32;
E_svd_ori__ = zeros(n_svd,1+n_shuffle);
for nshuffle=0:n_shuffle-1+1;
if (mod(nshuffle,8)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
if (nshuffle==0);
tmp_E_val_ = log(max(1e-12,E_val_));
tmp_svd_ = svds(tmp_E_val_,n_svd);
end;%if (nshuffle==0);
if (nshuffle> 0); tmp_E_val_ = log(1e-12)*ones(n_u,n_E_GENE);
rng(1024*nshuffle);
tmp_p_ = randperm(numel(E_index_filled_));
tmp_E_val_(1+E_index_filled_) = log(max(1e-12,E_val_(1+E_index_filled_(tmp_p_))));
tmp_svd_ = svds(tmp_E_val_,n_svd);
end;%if (nshuffle> 0); 
E_svd_ori__(:,1+nshuffle) = tmp_svd_;
clear tmp_E_val_ tmp_svd_ tmp_p_;
end;%for nshuffle=0:n_shuffle-1+1;
save(fname_mat ...
     ,'n_shuffle','n_svd','E_svd_ori__' ...
     );
end;%if (~exist(fname_mat,'file'));
load(fname_mat);
%%%%%%%%;
fname_mat = sprintf('%s/I_svd_ori_.mat',dir_mat);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
I_index_missing_ = efind(I_val_<=0); I_index_filled_ = setdiff(0:n_u*n_I_GENE-1,I_index_missing_);
%hist(log(max(1e-12,I_val_(1+I_index_filled_))),128);
n_shuffle = 32; n_svd = 32;
I_svd_ori__ = zeros(n_svd,1+n_shuffle);
for nshuffle=0:n_shuffle-1+1;
if (mod(nshuffle,8)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
if (nshuffle==0);
tmp_I_val_ = log(max(1e-12,I_val_));
tmp_svd_ = svds(tmp_I_val_,n_svd);
end;%if (nshuffle==0);
if (nshuffle> 0); 
tmp_I_val_ = log(1e-12)*ones(n_u,n_I_GENE);
rng(1024*nshuffle); tmp_p_ = randperm(numel(I_index_filled_));
tmp_I_val_(1+I_index_filled_) = log(max(1e-12,I_val_(1+I_index_filled_(tmp_p_))));
tmp_svd_ = svds(tmp_I_val_,n_svd);
end;%if (nshuffle> 0); 
I_svd_ori__(:,1+nshuffle) = tmp_svd_;
clear tmp_I_val_ tmp_svd_ tmp_p_;
end;%for nshuffle=0:n_shuffle-1+1;
save(fname_mat ...
     ,'n_shuffle','n_svd','I_svd_ori__' ...
     );
end;%if (~exist(fname_mat,'file'));
load(fname_mat);
%%%%%%%%;
% Now examine the singular-values. ;
%%%%%%%%;
E_svd_avg_ = mean(E_svd_ori__(:,2:end),2);
E_svd_std_ = std(E_svd_ori__(:,2:end),1,2);
E_svd_Z__ = (E_svd_ori__ - repmat(E_svd_avg_,[1,1+n_shuffle]))./repmat(E_svd_std_,[1,1+n_shuffle]);
E_svd_nlpZ__ = -z_to_lp(E_svd_Z__);
%%%%;
I_svd_avg_ = mean(I_svd_ori__(:,2:end),2);
I_svd_std_ = std(I_svd_ori__(:,2:end),1,2);
I_svd_Z__ = (I_svd_ori__ - repmat(I_svd_avg_,[1,1+n_shuffle]))./repmat(I_svd_std_,[1,1+n_shuffle]);
I_svd_nlpZ__ = -z_to_lp(I_svd_Z__);
%%%%;
fname_fig = sprintf('%s/X_svd_ori_FIGA',dir_jpg);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;
subplot(3,2,1); hold on;
plot(1:n_svd,E_svd_avg_ - 2.0*E_svd_std_,'k-','LineWidth',0.5);
plot(1:n_svd,E_svd_avg_ - 1.0*E_svd_std_,'k-','LineWidth',1.0);
plot(1:n_svd,E_svd_avg_ + 0.0*E_svd_std_,'k-','LineWidth',2.0);
plot(1:n_svd,E_svd_avg_ + 1.0*E_svd_std_,'k-','LineWidth',1.0);
plot(1:n_svd,E_svd_avg_ + 2.0*E_svd_std_,'k-','LineWidth',0.5);
plot(1:n_svd,E_svd_ori__(:,1),'ro-','LineWidth',2); grid on;
xlabel('rank');ylabel('sigma'); title('E_svd_ori_','Interpreter','none');
subplot(3,2,2); hold on;
plot(1:n_svd,I_svd_avg_ - 2.0*I_svd_std_,'k-','LineWidth',0.5);
plot(1:n_svd,I_svd_avg_ - 1.0*I_svd_std_,'k-','LineWidth',1.0);
plot(1:n_svd,I_svd_avg_ + 0.0*I_svd_std_,'k-','LineWidth',2.0);
plot(1:n_svd,I_svd_avg_ + 1.0*I_svd_std_,'k-','LineWidth',1.0);
plot(1:n_svd,I_svd_avg_ + 2.0*I_svd_std_,'k-','LineWidth',0.5);
plot(1:n_svd,I_svd_ori__(:,1),'ro-','LineWidth',2); grid on;
xlabel('rank');ylabel('sigma'); title('I_svd_ori_','Interpreter','none');
subplot(3,2,3); hold on;
plot(1:n_svd,-2.0*ones(1,n_svd),'k-','LineWidth',0.5);
plot(1:n_svd,-1.0*ones(1,n_svd),'k-','LineWidth',1.0);
plot(1:n_svd,+0.0*ones(1,n_svd),'k-','LineWidth',2.0);
plot(1:n_svd,+1.0*ones(1,n_svd),'k-','LineWidth',1.0);
plot(1:n_svd,+2.0*ones(1,n_svd),'k-','LineWidth',0.5);
plot(1:n_svd,E_svd_Z__(:,1),'ro-','LineWidth',2); grid on;
xlabel('rank');ylabel('Z-score'); title('E_svd_Z_','Interpreter','none');
subplot(3,2,4); hold on;
plot(1:n_svd,-2.0*ones(1,n_svd),'k-','LineWidth',0.5);
plot(1:n_svd,-1.0*ones(1,n_svd),'k-','LineWidth',1.0);
plot(1:n_svd,+0.0*ones(1,n_svd),'k-','LineWidth',2.0);
plot(1:n_svd,+1.0*ones(1,n_svd),'k-','LineWidth',1.0);
plot(1:n_svd,+2.0*ones(1,n_svd),'k-','LineWidth',0.5);
plot(1:n_svd,I_svd_Z__(:,1),'ro-','LineWidth',2); grid on;
xlabel('rank');ylabel('Z-score'); title('I_svd_Z_','Interpreter','none');
subplot(3,2,5); hold on;
plot(1:n_svd,ones(1,n_svd)*-log(0.0500),'k-','LineWidth',0.5);
plot(1:n_svd,ones(1,n_svd)*-log(0.0100),'k-','LineWidth',1.0);
plot(1:n_svd,ones(1,n_svd)*-log(0.0010),'k-','LineWidth',2.0);
plot(1:n_svd,ones(1,n_svd)*-log(0.0001),'k-','LineWidth',4.0);
plot(1:n_svd,E_svd_nlpZ__(:,1),'ro-','LineWidth',2); grid on;
legend({'0.05','0.01','0.001','0.0001','E'});
ylim([0,128]); xlabel('rank');ylabel('nlpZ'); title('E_svd_nlpZ_','Interpreter','none');
subplot(3,2,6); hold on;
plot(1:n_svd,ones(1,n_svd)*-log(0.0500),'k-','LineWidth',0.5);
plot(1:n_svd,ones(1,n_svd)*-log(0.0100),'k-','LineWidth',1.0);
plot(1:n_svd,ones(1,n_svd)*-log(0.0010),'k-','LineWidth',2.0);
plot(1:n_svd,ones(1,n_svd)*-log(0.0001),'k-','LineWidth',4.0);
plot(1:n_svd,I_svd_nlpZ__(:,1),'ro-','LineWidth',2); grid on;
legend({'0.05','0.01','0.001','0.0001','E'});
ylim([0,128]); xlabel('rank');ylabel('nlpz'); title('I_svd_nlpZ_','Interpreter','none');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Conclusion: estimated rank should be 16. ;
%%%%%%%%;
rank_estimate_E = max(find(E_svd_nlpZ__(:,1)> -log(0.01)));
rank_estimate_I = max(find(I_svd_nlpZ__(:,1)> -log(0.01)));
disp(sprintf(' %% rank_estimate_E %d rank_estimate_I %d',rank_estimate_E,rank_estimate_I));

%%%%%%%%;
% Is missingness itself correlated with the cell-types? ;
%%%%%%%%;
E_index_missing_ = efind(E_val_<=0); E_index_filled_ = setdiff(0:n_u*n_E_GENE-1,E_index_missing_);
E_fill_ = zeros(size(E_val_)); E_fill_(1+E_index_filled_)=1; E_fill_ = E_fill_ - mean(E_fill_,'all');
[tmp_U_E_fill__,tmp_S_E_fill__,tmp_V_E_fill__] = svds(E_fill_,rank_estimate_E);
%save(sprintf('%s/E_fill_.mat',dir_mat),'E_fill_');
I_index_missing_ = efind(I_val_<=0); I_index_filled_ = setdiff(0:n_u*n_I_GENE-1,I_index_missing_);
I_fill_ = zeros(size(I_val_)); I_fill_(1+I_index_filled_)=1; I_fill_ = I_fill_ - mean(I_fill_,'all');
[tmp_U_I_fill__,tmp_S_I_fill__,tmp_V_I_fill__] = svds(I_fill_,rank_estimate_I);
%save(sprintf('%s/I_fill_.mat',dir_mat),'I_fill_');
%%%%%%%%;
flag_plot=0;
if (flag_plot);
figure(2);clf;
colormap('lines');
subplot(1,2,1); cla;
scatter3(tmp_U_E_fill__(:,1),tmp_U_E_fill__(:,2),tmp_U_E_fill__(:,3),25,Sample_Label_sub_,'filled'); 
axis vis3d; title('E_fill','Interpreter','none');
subplot(1,2,2); cla;
scatter3(tmp_U_I_fill__(:,1),tmp_U_I_fill__(:,2),tmp_U_I_fill__(:,3),25,Sample_Label_sub_,'filled'); 
axis vis3d; title('I_fill','Interpreter','none');
figbig;
close(gcf);
end;%flag_plot=0;
%%%%%%%%;
% Conclusion: could be; will investigate later on. ;
%%%%%%%%; 
tmp_ = E_fill_; fname_tsv = sprintf('%s/dir_mat/E_fill_.tsv',dir_trunk); 
if (~exist(fname_tsv,'file')); disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs'); end;
clear tmp_;
tmp_ = I_fill_; fname_tsv = sprintf('%s/dir_mat/I_fill_.tsv',dir_trunk); 
if (~exist(fname_tsv,'file')); disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs'); end;
clear tmp_;

%%%%%%%%;
% attempt svd-based imputation. ;
%%%%%%%%;
tmp_rseed = 1; %<-- may be different from tmp_rseed = 1 or 2, which is used in test_loader_29.m and 29b, respectively ;
tmp_n_iteration = 256;
tmp_tolerance = 1e-3;
if flag_load==0;
%%%%%%%%;
E_ij_missing_ = find(E_val_<=0); E_ij_filled_ = setdiff(1:n_u*n_E_GENE,E_ij_missing_); 
E_liX_prefix_ = sprintf('%s/dir_jpg/E_liX_impute_fit',dir_trunk);
[~,E_liXf_] = svd_impute_fit_3(log(max(1e-12,E_val_)),E_ij_missing_,rank_estimate_E,E_liX_prefix_,tmp_rseed,tmp_n_iteration,tmp_tolerance);
%save(sprintf('%s/E_liX_.mat',dir_mat),'E_liX_');
save(sprintf('%s/E_liXf_.mat',dir_mat),'E_liXf_');
%tmp_ = E_liX_; fname_tsv = sprintf('%s/dir_mat/E_liX_.tsv',dir_trunk); 
%disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = E_liXf_; fname_tsv = sprintf('%s/dir_mat/E_liXf_.tsv',dir_trunk); 
disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs');
%%%%%%%%;
% demonstrate structure imposed by typical imputation. ;
%%%%%%%%;
flag_check=0;
if flag_check;
[tmp_M_] = sparse(1:n_u,randperm(n_u),1); %<-- random permutation. ;
tmp_E_ = impute_randperm(log(max(1e-12,E_val_)),E_ij_missing_);
tmp_sE_r_ = svd(tmp_E_); tmp_sME_r_ = svd(impute_randperm(tmp_M_*tmp_E_,E_ij_missing_));
tmp_E_ = impute_median(log(max(1e-12,E_val_)),E_ij_missing_);
tmp_sE_m_ = svd(tmp_E_); tmp_sME_m_ = svd(impute_median(tmp_M_*tmp_E_,E_ij_missing_));
tmp_E_ = impute_knn(0,log(max(1e-12,E_val_)),E_ij_missing_);
tmp_sE_n_ = svd(tmp_E_); tmp_sME_n_ = svd(impute_knn(0,tmp_M_*tmp_E_,E_ij_missing_));
flag_plot=1;
if flag_plot;
figure(4);clf;
subplot(2,3,1); plot(1:n_u,tmp_sE_r_,'ko',1:n_u,tmp_sME_r_,'rx'); xlim([1,n_u]);
xlabel('rank'); ylabel('singular value'); title('random imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthEast');
subplot(2,3,4); plot(1:n_u,log(tmp_sE_r_),'ko',1:n_u,log(tmp_sME_r_),'rx'); xlim([1,100]);
xlabel('rank'); ylabel('log(singular value)'); title('random imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthEast');
subplot(2,3,2); plot(1:n_u,tmp_sE_m_,'ko',1:n_u,tmp_sME_r_,'rx'); xlim([1,n_u]);
xlabel('rank'); ylabel('singular value'); title('median imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthEast');
subplot(2,3,5); plot(1:n_u,log(tmp_sE_m_),'ko',1:n_u,log(tmp_sME_r_),'rx'); xlim([1,100]);
xlabel('rank'); ylabel('log(singular value)'); title('median imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthEast');
subplot(2,3,3); plot(1:n_u,tmp_sE_n_,'ko',1:n_u,tmp_sME_r_,'rx'); xlim([1,n_u]);
xlabel('rank'); ylabel('singular value'); title('1nn imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthEast');
subplot(2,3,6); plot(1:n_u,log(tmp_sE_n_),'ko',1:n_u,log(tmp_sME_r_),'rx'); xlim([1,100]);
xlabel('rank'); ylabel('log(singular value)'); title('1nn imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthEast');
figbig;
fname_base = sprintf('%s/dir_jpg/E_knn_impute',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
close(gcf);
clear tmp_M_ tmp_E_ ;
end;%if flag_plot;
end;%if flag_check;
%%%%%%%%;
I_ij_missing_ = find(I_val_<=0); I_ij_filled_ = setdiff(1:n_u*n_I_GENE,I_ij_missing_); 
I_liX_prefix_ = sprintf('%s/dir_jpg/I_liX_impute_fit',dir_trunk);
[~,I_liXf_] = svd_impute_fit_3(log(max(1e-12,I_val_)),I_ij_missing_,rank_estimate_I,I_liX_prefix_,tmp_rseed,tmp_n_iteration,tmp_tolerance);
%save(sprintf('%s/I_liX_.mat',dir_mat),'I_liX_');
save(sprintf('%s/I_liXf_.mat',dir_mat),'I_liXf_');
%tmp_ = I_liX_; fname_tsv = sprintf('%s/dir_mat/I_liX_.tsv',dir_trunk); 
%disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = I_liXf_; fname_tsv = sprintf('%s/dir_mat/I_liXf_.tsv',dir_trunk); 
disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs');
%%%%%%%%;
% demonstrate structure imposed by typical imputation. ;
%%%%%%%%;
flag_check=0;
if flag_check;
[tmp_M_] = sparse(1:n_u,randperm(n_u),1); %<-- random permutation. ;
tmp_I_ = impute_randperm(log(max(1e-12,I_val_)),I_ij_missing_);
tmp_sI_r_ = svd(tmp_I_); tmp_sMI_r_ = svd(impute_randperm(tmp_M_*tmp_I_,I_ij_missing_));
tmp_I_ = impute_median(log(max(1e-12,I_val_)),I_ij_missing_);
tmp_sI_m_ = svd(tmp_I_); tmp_sMI_m_ = svd(impute_median(tmp_M_*tmp_I_,I_ij_missing_));
tmp_I_ = impute_knn(0,log(max(1e-12,I_val_)),I_ij_missing_);
tmp_sI_n_ = svd(tmp_I_); tmp_sMI_n_ = svd(impute_knn(0,tmp_M_*tmp_I_,I_ij_missing_));
flag_plot=1;
if flag_plot;
figure(4);clf;
subplot(2,3,1); plot(1:n_u,tmp_sI_r_,'ko',1:n_u,tmp_sMI_r_,'rx'); xlim([1,n_u]);
xlabel('rank'); ylabel('singular value'); title('random imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthIast');
subplot(2,3,4); plot(1:n_u,log(tmp_sI_r_),'ko',1:n_u,log(tmp_sMI_r_),'rx'); xlim([1,100]);
xlabel('rank'); ylabel('log(singular value)'); title('random imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthIast');
subplot(2,3,2); plot(1:n_u,tmp_sI_m_,'ko',1:n_u,tmp_sMI_r_,'rx'); xlim([1,n_u]);
xlabel('rank'); ylabel('singular value'); title('median imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthIast');
subplot(2,3,5); plot(1:n_u,log(tmp_sI_m_),'ko',1:n_u,log(tmp_sMI_r_),'rx'); xlim([1,100]);
xlabel('rank'); ylabel('log(singular value)'); title('median imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthIast');
subplot(2,3,3); plot(1:n_u,tmp_sI_n_,'ko',1:n_u,tmp_sMI_r_,'rx'); xlim([1,n_u]);
xlabel('rank'); ylabel('singular value'); title('1nn imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthIast');
subplot(2,3,6); plot(1:n_u,log(tmp_sI_n_),'ko',1:n_u,log(tmp_sMI_r_),'rx'); xlim([1,100]);
xlabel('rank'); ylabel('log(singular value)'); title('1nn imputation'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthIast');
figbig;
fname_base = sprintf('%s/dir_jpg/I_knn_impute',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
close(gcf);
clear tmp_M_ tmp_I_ ;
end;%if flag_plot;
end;%if flag_check;
%%%%%%%%;
end;%if flag_load==0;
if flag_load==1;
%load(sprintf('%s/E_liX_.mat',dir_mat),'E_liX_');
%load(sprintf('%s/I_liX_.mat',dir_mat),'I_liX_');
load(sprintf('%s/E_liXf_.mat',dir_mat),'E_liXf_');
load(sprintf('%s/I_liXf_.mat',dir_mat),'I_liXf_');
end;%if flag_load==1;

%%%%%%%%;
% Do not apply max/min cutoffs at 01 and 99 percentiles. ;
% Only apply centered log ratio. ;
%%%%%%%%;
if flag_E; 
tmp_99 = prctile(E_liXf_(:),99); tmp_01 = prctile(E_liXf_(:), 1);
%E_liXf_clogr_ = mean_center_0(max(tmp_01,min(tmp_99,E_liXf_)),'col');
E_liXf_clogr_ = mean_center_0(E_liXf_,'col');
%tmp_ = E_liXf_clogr_; fname_tsv = sprintf('%s/dir_mat/E_liXf_clogr_.tsv',dir_trunk); 
%if (~exist(fname_tsv,'file')); disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs'); end;
%clear tmp_;
end;%if flag_E; 
if flag_I; 
tmp_99 = prctile(I_liXf_(:),99); tmp_01 = prctile(I_liXf_(:), 1);
%I_liXf_clogr_ = mean_center_0(max(tmp_01,min(tmp_99,I_liXf_)),'col');
I_liXf_clogr_ = mean_center_0(I_liXf_,'col');
%tmp_ = I_liXf_clogr_; fname_tsv = sprintf('%s/dir_mat/I_liXf_clogr_.tsv',dir_trunk); 
%if (~exist(fname_tsv,'file')); disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs'); end;
clear tmp_;
end;%if flag_I; 

%%%%%%%%;
% Generating rank-normalized matrices. ;
%%%%%%%%;
if flag_E; 
%E_liXf_rankn_ = rank_normalize_0(E_liXf_clogr_,'row'); 
%tmp_ = E_liXf_rankn_; fname_tsv = sprintf('%s/dir_mat/E_liXf_rankn_.tsv',dir_trunk); 
%if (~exist(fname_tsv,'file')); disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs'); end;
%clear tmp_;
end;%if flag_E; 
if flag_I; 
%I_liXf_rankn_ = rank_normalize_0(I_liXf_clogr_,'row'); 
%tmp_ = I_liXf_rankn_; fname_tsv = sprintf('%s/dir_mat/I_liXf_rankn_.tsv',dir_trunk); 
%if (~exist(fname_tsv,'file')); disp(sprintf(' %% writing %s',fname_tsv)); save(fname_tsv,'tmp_','-ascii','-tabs'); end;
%clear tmp_;
end;%if flag_I;

%%%%%%%%;
% Extracting glossary. ;
% The vector 'LC' will list the LinearCorrelationWithPositiveQuality. ;
%%%%%%%%;
str_table_name = sprintf('%s/metric_glossary.txt',dir_data);
G_ = readtable(str_table_name); n_G = size(G_,1);
G_h_ = zeros(255,n_G); 
for nG=1:n_G;
G_h_(:,nG) = hist(cast(cast(G_{nG,1}{1},'uint8'),'double'),1:255);
end;%for nG=1:n_G;
G_h_mask_ = zeros(255,1);
G_h_mask_(cast('0','uint8'):cast('9','uint8')) = 1;
G_h_mask_(cast('A','uint8'):cast('Z','uint8')) = 1;
G_h_mask_(cast('a','uint8'):cast('z','uint8')) = 1;
G_h_ = G_h_.*repmat(G_h_mask_,1,n_G);
for nG=1:n_G;
G_h_(:,nG) = G_h_(:,nG)/norm(G_h_(:,nG),'fro');
end;%for nG=1:n_G;

%%%%%%%%;
% Extracting categorical covariates. ;
%%%%%%%%;
B_col_val_ = setdiff(7:23,[9,10,11,13,15,16,21,22,23]);
% 09 --> cell [1] x3ClassPrediction: 
% 10 <-- FailConfidencSum (double);
% 13 <-- SampleCount (double) ;
% 11 --> cell [47] max_leaf: 
% 15 --> cell [283] well: 
% 16 --> cell [9] batch: 
% 21 --> double [48] SequencingLane: 
% 22 --> double [3] FailState: 
% 23 --> double [12] RunID: 
n_BCOV = length(B_col_val_);
disp(sprintf(' %% %% %% %% '));
for nB=1:n_BCOV;
nc = B_col_val_(nB);
tmp_str = C_.Properties.VariableNames{nc};
tmp_ij = find(strcmp(G_{:,1},tmp_str));
tmp_LC = 0;
if ( isempty(tmp_ij)); disp(sprintf(' %% Warning! %s not found in glossary',tmp_str)); end;
if (~isempty(tmp_ij)); 
tmp_LC_str = G_{tmp_ij,2};
if strcmp(tmp_LC_str,'NA'); tmp_LC = 0; end;
if strcmp(tmp_LC_str,'Unknown'); tmp_LC = 0; end;
if strcmp(tmp_LC_str,'Positive'); tmp_LC = +1; end;
if strcmp(tmp_LC_str,'Negative'); tmp_LC = -1; end;
end;%if (~isempty(tmp_ij));
tmp_type = class(C_{C_to_u_ID_,nc});
tmp_u_ = unique(C_{C_to_u_ID_,nc});
disp(sprintf(' %% %0.2d --> %s [%d] %s [LC %d]: ',nc,tmp_type,length(tmp_u_),tmp_str,tmp_LC));
for nu=1:length(tmp_u_);
if strcmp(tmp_type,'cell');   disp(sprintf(' %% %% %s: %d',tmp_u_{nu},length(find(strcmp(C_{C_to_u_ID_,nc},tmp_u_{nu}))))); end;
if strcmp(tmp_type,'double'); disp(sprintf(' %% %% %d: %d',tmp_u_(nu),length(find(C_{C_to_u_ID_,nc} == tmp_u_(nu))))); end;
end;%for nu=1:length(tmp_u_);
disp(sprintf(' %% %% %% %% '));
end;%for nB=1:n_BCOV;
clear tmp_str tmp_ij tmp_LC tmp_LC_str tmp_type tmp_u_ ;
%%%%%%%%;
% Note that all LC are 0 for the categorical covariates. ;
%%%%%%%%;

%%%%%%%%;
n_FACTOR = 0; B_val_ = zeros(n_u,0); B_VariableName_ = cell(0);
%%%%%%%%;
% Cluster_Grp_20161007:
%%%%%%%%;
nc=7; tmp_str = 'exc'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=7; tmp_str = 'glia'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=7; tmp_str = 'inh'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
%%%%%%%%;
% IsOutlier_20161007: 
%%%%%%%%;
nc=8; tmp_str = 'no'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=8; tmp_str = 'yes'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
%%%%%%%%;
% cell_class: 
%%%%%%%%;
nc=12; tmp_str = 'GABAergic'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=12; tmp_str = 'Glutamatergic'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=12; tmp_str = 'Non-neuronal'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
%%%%%%%%;
% BatchCount: 
%%%%%%%%;
nc=14; tmp_d = 1;
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(C_{C_to_u_ID_,nc}==tmp_d); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=14; tmp_d = 2;
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(C_{C_to_u_ID_,nc}==tmp_d); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=14; tmp_d = 3;
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(C_{C_to_u_ID_,nc}==tmp_d); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
%%%%%%%%;
% neun:
%%%%%%%%;
nc=17; tmp_str = '1NeuNN'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=17; tmp_str = '1NeuNP'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
%%%%%%%%;
% patient_id: 
%%%%%%%%;
nc=18; tmp_str1 = 'H200-1025'; tmp_str2 = 'H200.1025'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str1) | strcmp(C_{C_to_u_ID_,nc},tmp_str2)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=18; tmp_str1 = 'H200-1030'; tmp_str2 = 'H200.1030'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str1) | strcmp(C_{C_to_u_ID_,nc},tmp_str2)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
%%%%%%%%;
% region:
%%%%%%%%;
nc=19; tmp_str = 'FI'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=19; tmp_str = 'MTG'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
%%%%%%%%;
% layer: 
%%%%%%%%;
nc=20; tmp_str = '1'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
nc=20; tmp_str = '5'; 
B_VariableName_{1+n_FACTOR} = sprintf('%s_%s',C_.Properties.VariableNames{nc},tmp_str); tmp_ij_ = find(strcmp(C_{C_to_u_ID_,nc},tmp_str)); B_val_(tmp_ij_,1+n_FACTOR)=1; n_FACTOR = n_FACTOR+1;
%%%%%%%%;
clear tmp_ij_ tmp_str tmp_str1 tmp_str2 ;

%%%%%%%%;
% Ranking categorical covariates. ;
%%%%%%%%;
B_rank_ = rank_normalize_0(B_val_); n_BCOV = size(B_rank_,2);
%%%%%%%%;
% Because LC = 0 for all categorical covariates, we can assign B_LC_ now. ;
%%%%%%%%;
B_LC_ = zeros(n_BCOV,1);

%%%%%%%%;
% Extracting continuous covariates. ;
%%%%%%%%;
D_col_val_ = [10,13,24:size(C_,2)];
str_exclude_ = {'fastQCResult','x_P_E__RawSeq__','x_Num_1_NumberOfInputCoreGenes_CoreGenes_','x_Num_2_NumberOfInputCoreGenes_CoreGenes_','x_MitoCore_NumberOfInputCoreGenes_CoreGenes_','x_MitoCore13_NumberOfInputCoreGenes_CoreGenes_','x_P_E__Pretrimmed_ERCCAligned_'};
for ne=1:length(str_exclude_);
nx = find(strcmp(C_.Properties.VariableNames,str_exclude_{ne}));
disp(sprintf(' %% excluding %s <-- %d <-- %s',str_exclude_{ne},nx,C_.Properties.VariableNames{nx}));
D_col_val_ = setdiff(D_col_val_,nx);
end;%for ne=1:length(str_exclude_);
%%%%%%%%;
n_D_col_val = length(D_col_val_); n_DCOV = n_D_col_val;
D_val_ = zeros(n_u,n_DCOV);
D_LC_ = zeros(n_DCOV,1);
D_VariableName_ = cell(n_D_col_val,1);
for nc=1:n_D_col_val;
D_col_val = D_col_val_(nc);
D_VariableName_{nc} = C_.Properties.VariableNames{D_col_val};
tmp_ij = find(strcmp(G_{:,1},D_VariableName_{nc}));
tmp_LC = 0;
if ( isempty(tmp_ij)); 
disp(sprintf(' %% %s no exact match found in glossary',D_VariableName_{nc})); 
tmp_h_ = hist(cast(cast(D_VariableName_{nc},'uint8'),'double'),1:255);
tmp_h_ = tmp_h_.*transpose(G_h_mask_);
tmp_h_ = tmp_h_/norm(tmp_h_,'fro');
[~,tmp_ij] = max(tmp_h_*G_h_);
disp(sprintf(' %% %s closest match: %s (%d)',D_VariableName_{nc},G_{tmp_ij,1}{1},tmp_ij));
end;%if ( isempty(tmp_ij)); 
if (~isempty(tmp_ij)); 
tmp_LC_str = G_{tmp_ij,2};
if strcmp(tmp_LC_str,'NA'); tmp_LC = 0; end;
if strcmp(tmp_LC_str,'Unknown'); tmp_LC = 0; end;
if strcmp(tmp_LC_str,'Positive'); tmp_LC = +1; end;
if strcmp(tmp_LC_str,'Negative'); tmp_LC = -1; end;
disp(sprintf(' %% %s found in glossary: %d',D_VariableName_{nc},tmp_LC));
end;%if (~isempty(tmp_ij)); 
D_LC_(nc) = tmp_LC;
tmp_D_val_ = C_(C_to_u_ID_,D_col_val); 
if strcmp(class(tmp_D_val_{1,1}),'double'); D_val_(:,nc) = tmp_D_val_{1:end,1};
elseif strcmp(class(tmp_D_val_{1,1}),'cell'); D_val_(:,nc) = cellfun(@str2num,tmp_D_val_{1:end,1}); 
else disp(sprintf(' %% nc %d D_col_val %d LC %d class %s',nc,D_col_val,D_LC_(nc),class(tmp_D_val_{1,1}))); end;
end;%for nc=1:n_D_col_val;
D_rank_ = rank_normalize_0(D_val_);
clear tmp_ij tmp_LC tmp_h_ tmp_LC_str tmp_D_val_ ;

%%%%%%%%;
% Stacking together B_rank_ and D_rank_. ;
%%%%%%%%;
n_CCOV = n_BCOV + n_DCOV;
C_VariableName_ = {B_VariableName_{:},D_VariableName_{:}};
C_val_ = [B_val_ , D_val_];
C_rank_ = [B_rank_ , D_rank_];
C_LC_ = [B_LC_ ; D_LC_];
C_VariableName_LC_ = C_VariableName_;
for nCCOV=1:n_CCOV;
if C_LC_(nCCOV) == -1; tmp_suffix = '[ - ]'; end;
if C_LC_(nCCOV) ==  0; tmp_suffix = '[    ]'; end;
if C_LC_(nCCOV) == +1; tmp_suffix = '[+]'; end;
C_VariableName_LC_{nCCOV} = sprintf('%s_%s',C_VariableName_{nCCOV},tmp_suffix);
end;%for nCCOV=1:n_CCOV;
clear tmp_suffix ;
%%%%%%%%;
% Note how clearly the missingness is correlated with some of the covariates. ;
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/E_fill_C_rank_bar',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
bar(corr(mean(E_fill_,2),C_rank_));
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_); xtickangle(90);
title('E_fill_ vs C_rank_','Interpreter','none');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/I_fill_C_rank_bar',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
bar(corr(mean(I_fill_,2),C_rank_));
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_); xtickangle(90);
title('I_fill_ vs C_rank_','Interpreter','none');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% estimate rank of C_rank_. ;
%%%%%%%%;
fname_mat = sprintf('%s/S_c127_original.mat',dir_mat);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
[ ...
 rank_estimate_sample ...
,s_B_nlp_ ...
,svd_sample__ ...
,eig_tw__ ...
,eig_B_ ...
,kta_opt__ ...
,h_x__ ...
,h_eig__ ...
,h_opt__ ...
] = ...
rank_estimate_sample_0(mean_center_0(C_rank_));
save(fname_mat ...
,'rank_estimate_sample' ...
,'s_B_nlp_' ...
,'svd_sample__' ...
,'eig_tw__' ...
,'eig_B_' ...
,'kta_opt__' ...
,'h_x__' ...
,'h_eig__' ...
,'h_opt__' ...
);
end;%if (~exist(fname_mat,'file'));
tmp_ = load(fname_mat);
rank_estimate_C = tmp_.rank_estimate_sample;
clear tmp_;

%%%%%%%%;
% Fully supervised clustering of covariates. ;
% i.e., direct measure of QC-differential-expression. ;
% Are any of these significant? ;
%%%%%%%%;
fname_mat = sprintf('%s/test_loader_cluster_C_rank_AB.mat',dir_mat);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
n_iteration = 1024;
[AB_] = test_loader_cluster_C_AB_1(n_iteration,str_Sample_Label_,C_rank_);
save(fname_mat,'n_iteration','str_Sample_Label_','AB_');
end;%if (~exist(fname_mat,'file'));
load(fname_mat);
%%%%%%%%;
flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dir_jpg/CcabsZ_shuffle_heatmap',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
Ccavg0_ = mean_center_0(AB_.C_rank_avg_);
[tmp_UCcavg0_,tmp_SCcavg0_,tmp_VCcavg0_] = svd(Ccavg0_);
CcabsZ_ = mean_center_0(AB_.AB_C_absZ_);
[tmp_UCcabsZ_,tmp_SCcabsZ_,tmp_VCcabsZ_] = svd(CcabsZ_);
CcrawZ_ = mean_center_0(AB_.AB_C_rawZ_);
[tmp_UCcrawZ_,tmp_SCcrawZ_,tmp_VCcrawZ_] = svd(CcrawZ_);
%subplot(1,3,1); scatter(tmp_VCcavg0_(:,1),tmp_VCcavg0_(:,2),15,'o','filled'); axisnotick; title('Ccavg0');
%subplot(1,3,2); scatter(tmp_VCcabsZ_(:,1),tmp_VCcabsZ_(:,2),15,'o','filled'); axisnotick; title('CcabsZ');
%subplot(1,3,3); scatter(tmp_VCcrawZ_(:,1),tmp_VCcrawZ_(:,2),15,'o','filled'); axisnotick; title('CcrawZ');
%%%%%%%%;
% quickly test statistical significance of CcabsZ_. ;
%%%%%%%%;
n_shuffle = 17; n_iteration = 256;
CcabsZ__ = zeros(n_Sample_Label^2,n_CCOV,1+n_shuffle);
for nshuffle=0:n_shuffle;
if (nshuffle==0); tmp_Q_ = eye(n_CCOV); end;
if (nshuffle> 0); rng(1024*nshuffle); [tmp_Q_,~] = qr(randn(n_CCOV)); end;
%disp(num2str(tmp_Q_(1,1:10)));
tmp_C_rank_ = mean_center_0(C_rank_)*tmp_Q_; 
tmp_AB_ = test_loader_cluster_C_AB_1(n_iteration,str_Sample_Label_,tmp_C_rank_);
CcabsZ__(:,:,1+nshuffle) = tmp_AB_.AB_C_absZ_;
end;%for nshuffle=0:n_shuffle;
lim_ = mean(CcabsZ__(:,:,1),'all') + std(CcabsZ__(:,:,1),1,'all')*3.5*[-1,+1];
for nshuffle=0:n_shuffle;
[tmp_UCcabsZ_,tmp_SCcabsZ_,tmp_VCcabsZ_] = svds(CcabsZ__(:,:,1+nshuffle),2);
[~,tmp_UCcabsZ_ij_] = sort(tmp_UCcabsZ_(:,1)); [~,tmp_VCcabsZ_ij_] = sort(tmp_VCcabsZ_(:,1));
figure(1); subplot(3,6,1+nshuffle); imagesc(CcabsZ__(tmp_UCcabsZ_ij_,tmp_VCcabsZ_ij_,1+nshuffle),lim_);
figure(2); subplot(3,6,1+nshuffle); scatter(tmp_VCcabsZ_(:,1),tmp_VCcabsZ_(:,2),'o','filled');
end;%for nshuffle=0:n_shuffle;
figure(1);colormap(colormap_beach());figbig;
figure(2);figbig;
figure(1);
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
figure(2);
fname_fig_2 = sprintf('%s/dir_jpg/CcabsZ_shuffle_scatter',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_fig_2));
print('-djpeg',sprintf('%s.jpg',fname_fig_2));
print('-depsc',sprintf('%s.eps',fname_fig_2));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Now retest CcabsZ_ using q3dcluster. ;
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/CcabsZ_shuffle_trace',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
n_shuffle = 17; n_iteration = 256;
QR__ = zeros(n_CCOV,1+n_shuffle);
QC__ = zeros(n_CCOV,1+n_shuffle);
r_rem_ = zeros(n_CCOV,1);
c_rem_ = zeros(n_CCOV,1);
for nshuffle=0:n_shuffle;
if (nshuffle==0); tmp_Q_ = eye(n_CCOV); end;
if (nshuffle> 0); rng(1024*nshuffle); [tmp_Q_,~] = qr(randn(n_CCOV)); end;
tmp_C_rank_ = mean_center_0(C_rank_)*tmp_Q_; 
tmp_AB_ = test_loader_cluster_C_AB_1(n_iteration,str_Sample_Label_,tmp_C_rank_);
[out_xdrop_,trace_] = q3dcluster_nonbinary_AAAA_ver0(transpose(tmp_AB_.AB_C_absZ_));
n_iter = size(trace_,1);
QR__(1:n_iter,1+nshuffle) = trace_(:,4);
QC__(1:n_iter,1+nshuffle) = trace_(:,5);
if (nshuffle==0); r_rem_(1:n_iter) = trace_(:,2); end;
if (nshuffle==0); c_rem_(1:n_iter) = trace_(:,3); end;
end;%for nshuffle=0:n_shuffle;
QR_avg_ = mean(QR__(:,2:end),2); QR_std_ = std(QR__(:,2:end),1,2);
ZR__ = (QR__ - repmat(QR_avg_,1,1+n_shuffle))./repmat(QR_std_,1,1+n_shuffle);
QC_avg_ = mean(QC__(:,2:end),2); QC_std_ = std(QC__(:,2:end),1,2);
ZC__ = (QC__ - repmat(QC_avg_,1,1+n_shuffle))./repmat(QC_std_,1,1+n_shuffle);
figure(1);
subplot(2,2,1);plot(1:n_iter,ZR__(1:n_iter,2:end),'k.-',1:n_iter,ZR__(1:n_iter,1),'ro-');
xlabel('iter'); ylabel('ZR');
subplot(2,2,2);plot(1:n_iter,ZC__(1:n_iter,2:end),'k.-',1:n_iter,ZC__(1:n_iter,1),'ro-');
xlabel('iter'); ylabel('ZC');
subplot(2,2,3);plot(1:n_iter,QR__(1:n_iter,2:end),'k.-',1:n_iter,QR__(1:n_iter,1),'ro-');
xlabel('iter'); ylabel('QR');
subplot(2,2,4);plot(1:n_iter,QC__(1:n_iter,2:end),'k.-',1:n_iter,QC__(1:n_iter,1),'ro-');
xlabel('iter'); ylabel('QC');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% dexnb. ;
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/Ccavg0_dexnb_scatter_FIGA',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 2;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('Ccavg0_%s',str_xfix); 
dir_out = []; E_array_base_ = transpose(Ccavg0_); E_array_r_ij_ = []; E_array_c_ij_ = [];
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
%ZRimax_label_ = test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
p_set_ = [0.05]*0.5.^[0:3];
for np=1:length(p_set_);
tmp_p_set = p_set_(np);
tmp_ZRimax_label_ = test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,tmp_p_set,n_member_lob);
subplot(2,2,np);
scatter(tmp_VCcavg0_(:,1),tmp_VCcavg0_(:,2),15,label_str_to_enum_0(tmp_ZRimax_label_));
colormap('lines');
title(sprintf('p %0.6f',tmp_p_set));
axisnotick;
end;%for np=1:length(p_set_);
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/Ccavg0_laknb_scatter_FIGA',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 2;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('lakcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('Ccavg0_%s',str_xfix); 
dir_out = []; E_array_base_ = transpose(Ccavg0_); E_array_r_ij_ = []; E_array_c_ij_ = [];
test_loader_lakcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
%ZRimax_label_ = test_loader_lakcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
p_set_ = [0.05]*0.5.^[0:3];
for np=1:length(p_set_);
tmp_p_set = p_set_(np);
tmp_ZRimax_label_ = test_loader_lakcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,tmp_p_set,n_member_lob);
subplot(2,2,np);
scatter(tmp_VCcavg0_(:,1),tmp_VCcavg0_(:,2),15,label_str_to_enum_0(tmp_ZRimax_label_));
colormap('lines');
title(sprintf('p %0.6f',tmp_p_set));
axisnotick;
end;%for np=1:length(p_set_);
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/CcabsZ_tsne_FIGA',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
tmp_Q_ = eye(n_CCOV); %[tmp_Q_,~] = qr(randn(n_CCOV));
n_rank_CC = 2;
CcabsZ_n_sub_ = fast_tsne_dr_0(transpose(CcabsZ_*tmp_Q_),struct('rand_seed',1,'no_dims',n_rank_CC,'theta',0.5));
opts_isosplit5 = struct('K_init',n_CCOV,'isocut_threshold',1.0);
label_CcabsZ_ = transpose(isosplit5(transpose(CcabsZ_n_sub_),opts_isosplit5));
subplot(1,1,1); scatter(CcabsZ_n_sub_(:,1),CcabsZ_n_sub_(:,2),15,label_CcabsZ_,'filled'); colormap('lines');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/CcabsZ_dexnb_scatter_FIGA',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 2;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('CcabsZ_%s',str_xfix); 
dir_out = []; E_array_base_ = transpose(CcabsZ_); E_array_r_ij_ = []; E_array_c_ij_ = [];
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
%ZRimax_label_ = test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
p_set_ = [0.05]*0.5.^[0:5];
for np=1:length(p_set_);
tmp_p_set = p_set_(np);
tmp_ZRimax_label_ = test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,tmp_p_set,n_member_lob);
subplot(2,3,np);
scatter(tmp_VCcabsZ_(:,1),tmp_VCcabsZ_(:,2),15,label_str_to_enum_0(tmp_ZRimax_label_),'filled');
colormap('lines');
title(sprintf('p %0.6f #=%d',tmp_p_set,length(unique(tmp_ZRimax_label_))));
axisnotick;
end;%for np=1:length(p_set_);
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Print output: ;
%%%%%%%%;
fname_base = sprintf('%s/dir_mat/CcabsZ_dexnb_p%.3d.txt',dir_trunk,floor(1000*p_set));
if (~exist(fname_base,'file'));
disp(sprintf(' %% %s not found, creating',fname_base));
p_set = 0.05; 
ZRimax_label_ = test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
ZRimax_enum_ = label_str_to_enum_0(ZRimax_label_);
n_QCluster = length(unique(ZRimax_enum_));
fp = fopen(fname_base,'w');
tmp_str = (sprintf(' n_QCluster: %d',n_QCluster)); fprintf(fp,'%s\n',tmp_str); disp(tmp_str);
tmp_str = (sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')); fprintf(fp,'%s\n',tmp_str); disp(tmp_str);
for nQCluster=0:n_QCluster-1;
tmp_str = (sprintf(' QCluster #%d:',1+nQCluster)); fprintf(fp,'%s\n',tmp_str); disp(tmp_str);
tmp_str = (sprintf('\t%s\n',C_VariableName_{find(ZRimax_enum_==1+nQCluster)})); fprintf(fp,'%s\n',tmp_str); disp(tmp_str);
tmp_str = (sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')); fprintf(fp,'%s\n',tmp_str); disp(tmp_str);
end;%for nQCluster=0:n_QCluster-1;
fclose(fp);
disp(sprintf(' %% writing %s',fname_base));
end;%if (~exist(fname_base,'file'));
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/Table_CcabsZ_dexnb_FIGB',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
Table_CcabsZ_ = zeros(n_QCluster,n_Sample_Label);
for nQCluster=0:n_QCluster-1;
QCluster_ij_ = find(ZRimax_enum_==1+nQCluster);
n_l = length(QCluster_ij_);
for nl=1:n_l;
QC_ij = QCluster_ij_(nl);
ns=0;
for nSample_Label_A=1:n_Sample_Label;
for nSample_Label_B=1:n_Sample_Label;
tmp_lp = z_to_p_0(AB_.AB_C_absZ_(1+ns,QC_ij)); tmp_lp = -tmp_lp/(n_l*n_Sample_Label);
Table_CcabsZ_(1+nQCluster,nSample_Label_A) = Table_CcabsZ_(1+nQCluster,nSample_Label_A) + tmp_lp;
Table_CcabsZ_(1+nQCluster,nSample_Label_B) = Table_CcabsZ_(1+nQCluster,nSample_Label_B) + tmp_lp;
ns=ns+1;
end;%for nSample_Label_B=1:n_Sample_Label;
end;%for nSample_Label_A=1:n_Sample_Label;
end;%for nl=1:n_l;
end;%for nQCluster=0:n_QCluster-1;
subplot(1,1,1);
colormap(colormap_beach());
imagesc(Table_CcabsZ_); colorbar; figbig;
set(gca,'XTick',1:n_Sample_Label,'XTickLabel',unique(str_Sample_Label_)); xtickangle(90);
set(gca,'YTick',1:n_QCluster,'YTickLabel',1:n_QCluster);
title('QC-cluster vs sample-cluster -log(p-value)');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/CcrawZ_laknb_scatter_FIGA',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 2;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('lakcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('CcrawZ_%s',str_xfix); 
dir_out = []; E_array_base_ = transpose(CcrawZ_); E_array_r_ij_ = []; E_array_c_ij_ = [];
test_loader_lakcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
%ZRimax_label_ = test_loader_lakcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
p_set_ = [0.05]*0.5.^[0:3];
for np=1:length(p_set_);
tmp_p_set = p_set_(np);
tmp_ZRimax_label_ = test_loader_lakcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,tmp_p_set,n_member_lob);
subplot(2,2,np);
scatter(tmp_VCcrawZ_(:,1),tmp_VCcrawZ_(:,2),15,label_str_to_enum_0(tmp_ZRimax_label_));
colormap('lines');
title(sprintf('p %0.6f',tmp_p_set));
axisnotick;
end;%for np=1:length(p_set_);
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if flag_plot;

%%%%%%%%;
% spectral clustering of covariates. ;
% Verdict: Covariates are not significantly clustered by themselves! ;
% Rather: certain covariates are significantly associated with the clusters found after analyzing transcript data. ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
%%%%%%%%;
tmp_Q_ = eye(n_CCOV); %[tmp_Q_,~] = qr(randn(n_CCOV));
CC_ = corr(mean_center_0(C_rank_)*tmp_Q_);
%%%%%%%%;
% spectral-->isosplit5;
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/C_rank_svd_scatter_FIGC',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
[tmp_U_,tmp_S_,tmp_V_] = svds(CC_,n_CCOV); %<-- around 11 large principal values. ;
n_rank_CC = 12;
for nrank_CC=1:n_rank_CC;
opts_isosplit5 = struct('K_init',n_CCOV,'isocut_threshold',1.0);
label_C__{nrank_CC} = transpose(isosplit5(transpose(tmp_V_(:,1:nrank_CC)*tmp_S_(1:nrank_CC,1:nrank_CC)),opts_isosplit5));
end;%for nrank_CC=1:n_rank_CC;
for nrank_CC=1:n_rank_CC;
subplot(3,4,nrank_CC); scatter(tmp_V_(:,1),tmp_V_(:,2),15,label_C__{nrank_CC}); colormap('lines');
end;%for nrank_CC=1:n_rank_CC;
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% tsne50-->isosplit5;
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/C_rank_tsne00_scatter_FIGC',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
tmp_Q_ = eye(n_CCOV); %[tmp_Q_,~] = qr(randn(n_CCOV));
n_rank_CC = 2;
C_n_sub_ = fast_tsne_dr_0(transpose(mean_center_0(C_rank_)*tmp_Q_),struct('rand_seed',1,'no_dims',n_rank_CC,'theta',0.0));
opts_isosplit5 = struct('K_init',n_CCOV,'isocut_threshold',1.0);
label_C__{nrank_CC} = transpose(isosplit5(transpose(C_n_sub_),opts_isosplit5));
subplot(1,1,1); scatter(C_n_sub_(:,1),C_n_sub_(:,2),15,label_C__{nrank_CC}); colormap('lines');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% dexnb; %<-- no statistical significance !? ;
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/C_rank_dexnb_scatter_FIGC',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 2;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('C_rank_%s',str_xfix); 
dir_out = []; E_array_base_ = transpose(C_rank_); E_array_r_ij_ = []; E_array_c_ij_ = [];
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
ZRimax_label_ = test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3(dir_cluster,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
%%%%%%%%;
[~,tmp_CC_ij_] = sort(tmp_V_(:,1));
%subplot(1,2,1); plot(1:n_CCOV,diag(tmp_S_),'ko'); xlabel('rank'); ylabel('sigma'); title('spectrum of CC');
subplot(1,1,1); imagesc(CC_(tmp_CC_ij_,tmp_CC_ij_),[-1,+1]); colormap(colormap_beach()); 
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(tmp_CC_ij_)); xtickangle(90);
set(gca,'YTick',1:n_CCOV,'YTickLabel',C_VariableName_LC_(tmp_CC_ij_)); 
set(gca,'FontSize',6);
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Is any of this statistically significant? ;
%%%%%%%%;
fname_fig = sprintf('%s/dir_jpg/C_rank_svd_scatter_FIGD',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
n_iteration = 54; rng(255);
for niteration=1:n_iteration;
if (niteration==1); tmp_Q_ = eye(n_CCOV); else; [tmp_Q_,~] = qr(randn(n_CCOV)); end;
[tmp_U_,tmp_S_,tmp_V_] = svds(mean_center_0(C_rank_)*tmp_Q_,n_CCOV); %<-- around 11 large principal values. ;
subplot(6,9,niteration); scatter(tmp_V_(:,1),tmp_V_(:,2),4,'ko','filled');
set(gca,'XTick',[],'YTick',[]);
%title(sprintf('ni%d',niteration));
end;%for niteration=1:n_iteration;
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if flag_plot;

%%%%%%%%;
% spectral clustering of covariates after mahalanobization via E_liXf_clogr_. ;
% Verdict: Covariates are significantly clustered with respect to transcript data. ;
%%%%%%%%;
flag_compute = strcmp(platform,'access1');

%{
if flag_compute;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% li16f ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
verbose=1;
E_li16f_ = textread(sprintf('%s/E_li16f_.tsv',dir_mat)); E_li16f_ = E_li16f_(:,1:end-1);
X_infix = 'E_li16f';
X_ = E_li16f_;
rank_estimate_X = n_u;
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
E_li16f_clogr_ = mean_center_0(E_li16f_,'col');
X_infix = 'E_li16f_clogr';
X_ = E_li16f_clogr_;
rank_estimate_X = n_u;
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
I_li16f_ = textread(sprintf('%s/I_li16f_.tsv',dir_mat)); I_li16f_ = I_li16f_(:,1:end-1);
X_infix = 'I_li16f';
X_ = I_li16f_;
rank_estimate_X = n_u;
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
I_li16f_clogr_ = mean_center_0(I_li16f_,'col');
X_infix = 'I_li16f_clogr';
X_ = I_li16f_clogr_;
rank_estimate_X = n_u;
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
X_infix = 'EI_li16f';
X_ = [E_li16f_ , I_li16f_];
rank_estimate_X = n_u;
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
X_infix = 'EI_li16f_clogr';
X_ = [E_li16f_clogr_ , I_li16f_clogr_];
rank_estimate_X = n_u;
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;%if flag_compute;
%}

if flag_compute;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% liXf ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
verbose=1;
X_infix = 'E_liXf_clogr';
X_ = E_liXf_clogr_;
rank_estimate_X = max(n_u,rank_estimate_E);
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
X_infix = 'E_liXf';
X_ = E_liXf_;
rank_estimate_X = max(n_u,rank_estimate_E);
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
X_infix = 'I_liXf_clogr';
X_ = I_liXf_clogr_;
rank_estimate_X = max(n_u,rank_estimate_I);
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
X_infix = 'I_liXf';
X_ = I_liXf_;
rank_estimate_X = max(n_u,rank_estimate_I);
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
X_infix = 'EI_liXf_clogr';
X_ = [E_liXf_clogr_ , I_liXf_clogr_];
rank_estimate_X = max(n_u,rank_estimate_I);
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
X_infix = 'EI_liXf';
X_ = [E_liXf_ , I_liXf_];
rank_estimate_X = max(n_u,rank_estimate_I);
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_compute;

if flag_compute;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% fill ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
verbose=1;
X_infix = 'E_fill';
X_ = E_fill_;
rank_estimate_X = n_u;
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
X_infix = 'I_fill';
X_ = I_fill_;
rank_estimate_X = n_u;
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%;
verbose=1;
X_infix = 'EI_fill';
X_ = [E_fill_ , I_fill_];
rank_estimate_X = n_u;
flag_force_create = 0;
test_loader_helper_covariate_cluster_0( ...
 verbose ...
,flag_force_create ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,n_CCOV ...
,C_rank_ ...
,C_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
,AB_ ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_compute;

%%%%%%%%;
fname_tsv = sprintf('%s/u_ID_sub_.nsv',dir_mat);
if (~exist(fname_tsv,'file'));
disp(sprintf(' %% %s not found, creating',fname_tsv));
fp=fopen(fname_tsv,'w');
fprintf(fp,'%s\n',u_ID_{1:end});
fclose(fp);
end;%if (~exist(fname_tsv,'file'));
%%%%%%%%;
fname_tsv = sprintf('%s/str_Sample_Label_sub_.nsv',dir_mat);
if (~exist(fname_tsv,'file'));
disp(sprintf(' %% %s not found, creating',fname_tsv));
fp=fopen(fname_tsv,'w');
fprintf(fp,'%s\n',str_Sample_Label_{1:end});
fclose(fp);
end;%if (~exist(fname_tsv,'file'));
fname_tsv = sprintf('%s/C_VariableName_.tsv',dir_mat);
if (~exist(fname_tsv,'file'));
disp(sprintf(' %% %s not found, creating',fname_tsv));
fp=fopen(fname_tsv,'w');
fprintf(fp,'%s\t',C_VariableName_{1:end});
fclose(fp);
end;%if (~exist(fname_tsv,'file'));
fname_tsv = sprintf('%s/C_VariableName_.nsv',dir_mat);
if (~exist(fname_tsv,'file'));
disp(sprintf(' %% %s not found, creating',fname_tsv));
fp=fopen(fname_tsv,'w');
fprintf(fp,'%s\n',C_VariableName_{1:end});
fclose(fp);
end;%if (~exist(fname_tsv,'file'));
%%%%%%%%;
fname_tsv = sprintf('%s/E_GeneName_sub_.nsv',dir_mat);
if (~exist(fname_tsv,'file'));
disp(sprintf(' %% %s not found, creating',fname_tsv));
fp=fopen(fname_tsv,'w');
fprintf(fp,'%s\n',E_VariableName_sub_{1:end});
fclose(fp);
end;%if (~exist(fname_tsv,'file'));
%%%%%%%%;
fname_tsv = sprintf('%s/I_GeneName_sub_.nsv',dir_mat);
if (~exist(fname_tsv,'file'));
disp(sprintf(' %% %s not found, creating',fname_tsv));
fp=fopen(fname_tsv,'w');
fprintf(fp,'%s\n',I_VariableName_sub_{1:end});
fclose(fp);
end;%if (~exist(fname_tsv,'file'));
%%%%%%%%;

flag_compute = strcmp(platform,'access1');

if flag_compute;
%%%%%%%%;
% Now original clustering (with no correction). ;
%%%%%%%%;
date_diff_threshold = 1.0;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;
fp_label_A_ = fopen(sprintf('%s/dir_mat/str_CLabel_sub_.nsv',dir_trunk),'r');
str_label_A_ = textscan(fp_label_A_,'%s'); fclose(fp_label_A_);
n_u = numel(str_label_A_{1});
label_A_ = label_str_to_enum_0(str_label_A_{1});
%%%%%%%%;
%for prefix_normalization_ = {'li16f','li16f_clogr','liXf','liXf_clogr','fill'};
for prefix_normalization_ = {'liXf','liXf_clogr','fill'};
prefix_normalization = prefix_normalization_{1};
if (~isempty(strfind(prefix_normalization,'fill'))); E_xxxx_ = E_fill_ ; I_xxxx_ = I_fill_; end;%if (~isempty(strfind(prefix_normalization,'fill')));
if (~isempty(strfind(prefix_normalization,'li16f'))); E_xxxx_ = E_li16f_ ; I_xxxx_ = I_li16f_; end;%if (~isempty(strfind(prefix_normalization,'li16f')));
if (~isempty(strfind(prefix_normalization,'liXf'))); E_xxxx_ = E_liXf_ ; I_xxxx_ = I_liXf_; end;%if (~isempty(strfind(prefix_normalization,'liXf')));
if (~isempty(strfind(prefix_normalization,'clogr'))); E_xxxx_ = mean_center_0(E_xxxx_,'col'); I_xxxx_ = mean_center_0(I_xxxx_,'col'); end;%if (~isempty(strfind(prefix_normalization,'clogr')));
str_infix = prefix_normalization;
test_loader_cluster_wrap_4(dir_trunk,label_A_,E_xxxx_,I_xxxx_,str_infix,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
clear E_xxxx_ I_xxxx_ ;
end;%for prefix_normalization_ = {'li16f','li16f_clogr','liXf','liXf_clogr','fill'};
%%%%%%%%;
end;%if flag_compute;

if flag_compute;
%%%%%%%%;
% Now correcting for covariates: original. ;
%%%%%%%%;
date_diff_threshold = 1.0;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;
fp_label_A_ = fopen(sprintf('%s/dir_mat/str_CLabel_sub_.nsv',dir_trunk),'r');
str_label_A_ = textscan(fp_label_A_,'%s'); fclose(fp_label_A_);
n_u = numel(str_label_A_{1});
label_A_ = label_str_to_enum_0(str_label_A_{1});
%%%%%%%%;
%for prefix_normalization_ = {'li16f','li16f_clogr','liXf','liXf_clogr','fill'};
for prefix_normalization_ = {'liXf','liXf_clogr','fill'};
prefix_normalization = prefix_normalization_{1};
for prefix_covariate_ = {'c127'};
prefix_covariate = prefix_covariate_{1};
if (~isempty(strfind(prefix_normalization,'fill'))); E_xxxx_ = E_fill_ ; I_xxxx_ = I_fill_; end;%if (~isempty(strfind(prefix_normalization,'fill')));
if (~isempty(strfind(prefix_normalization,'li16f'))); E_xxxx_ = E_li16f_ ; I_xxxx_ = I_li16f_; end;%if (~isempty(strfind(prefix_normalization,'li16f')));
if (~isempty(strfind(prefix_normalization,'liXf'))); E_xxxx_ = E_liXf_ ; I_xxxx_ = I_liXf_; end;%if (~isempty(strfind(prefix_normalization,'liXf')));
if (~isempty(strfind(prefix_normalization,'clogr'))); E_xxxx_ = mean_center_0(E_xxxx_,'col'); I_xxxx_ = mean_center_0(I_xxxx_,'col'); end;%if (~isempty(strfind(prefix_normalization,'clogr')));
if (~isempty(strfind(prefix_covariate,'c127'))); C_xxxx_ = mean_center_0(C_rank_); end;%if (~isempty(strfind(prefix_covariate,'c127')));
tmp_fname_mat = sprintf('%s/S_%s_original.mat',dir_mat,prefix_covariate);
tmp_ = load(tmp_fname_mat); rank_estimate_C_xxxx = tmp_.rank_estimate_sample; clear tmp_;
[U_C_xxxx_,S_C_xxxx_,V_C_xxxx_] = svds(C_xxxx_,rank_estimate_C_xxxx);
E_xxxx_ = E_xxxx_ - U_C_xxxx_*(transpose(U_C_xxxx_)*E_xxxx_);  %<-- residual. ;
I_xxxx_ = I_xxxx_ - U_C_xxxx_*(transpose(U_C_xxxx_)*I_xxxx_);  %<-- residual. ;
str_infix = sprintf('%s_%s',prefix_normalization,prefix_covariate);
test_loader_cluster_wrap_4(dir_trunk,label_A_,E_xxxx_,I_xxxx_,str_infix,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
clear E_xxxx_ I_xxxx_ C_xxxx_ U_C_xxxx_ S_C_xxxx_ V_C_xxxx_;
end;%for prefix_covariate_ = {'c127'};
end;%for prefix_normalization_ = {'li16f','li16f_clogr','liXf','liXf_clogr','fill'};
%%%%%%%%;
end;%if flag_compute;

if flag_compute;
%%%%%%%%;
% Now correcting for covariates: mahalanobis. ;
%%%%%%%%;
date_diff_threshold = 1.0;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;
fp_label_A_ = fopen(sprintf('%s/dir_mat/str_CLabel_sub_.nsv',dir_trunk),'r');
str_label_A_ = textscan(fp_label_A_,'%s'); fclose(fp_label_A_);
n_u = numel(str_label_A_{1});
label_A_ = label_str_to_enum_0(str_label_A_{1});
%%%%%%%%;
for str_X_ = {'E','I','EI'};
str_X = str_X_{1};
%for prefix_normalization_ = {'li16f','li16f_clogr','liXf','liXf_clogr','fill'};
for prefix_normalization_ = {'liXf','liXf_clogr','fill'};
prefix_normalization = prefix_normalization_{1};
tmp_fname_tsv = sprintf('%s/n_QCluster_U%s_%s_C_mahalanobis.tsv',dir_mat,str_X,prefix_normalization);
n_QCluster = textread(tmp_fname_tsv); n_QCluster = n_QCluster(1);
disp(sprintf(' %% %s_%s n_QCluster %d',str_X,prefix_normalization,n_QCluster));
for nQCluster=0:n_QCluster-1;
%%%%%%%%;
tmp_fname_mat = sprintf('%s/S_c%d_U%s_%s_C_mahalanobis.mat',dir_mat,1+nQCluster,str_X,prefix_normalization);
tmp_ = load(tmp_fname_mat); rank_estimate_UXC = tmp_.rank_estimate_sample; clear tmp_;
tmp_fname_tsv = sprintf('%s/c%d_U%s_%s__.tsv',dir_mat,1+nQCluster,str_X,prefix_normalization);
tmp_C_ = textread(tmp_fname_tsv); tmp_C_ = tmp_C_(:,1:end-1);
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
rank_estimate_X = n_u;
[tmp_UX_,tmp_SX_,tmp_VX_] = svds(X_,rank_estimate_X); %<-- use full rank of X_. ;
tmp_UXC_ = inv(tmp_SX_)*transpose(tmp_UX_)*mean_center_0(tmp_C_);
[tmp_UUXC_,tmp_SUXC_,tmp_VUXC_] = svds(tmp_UXC_,rank_estimate_UXC);
X_ = (tmp_UX_*tmp_SX_)*(transpose(tmp_VX_) - tmp_UUXC_*(transpose(tmp_UUXC_)*transpose(tmp_VX_))); %<-- transpose(VX_) = UXC_*Z_  --> residual. ;
str_infix = sprintf('%s_UX_c%d',prefix_normalization,1+nQCluster);
test_loader_cluster_wrap_5(dir_trunk,label_A_,X_,str_X,str_infix,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
clear X_ tmp_C_ tmp_UX_ tmp_SX_ tmp_VX_ tmp_UXC_ tmp_UUXC_ tmp_SUXC_ tmp_VUXC_ ;
%%%%%%%%;
end;%for nQCluster=0:n_QCluster-1;
end;%for prefix_normalization_ = {'li16f','li16f_clogr','liXf','liXf_clogr','fill'};
end;%for str_X_ = {'E','I'};
%%%%%%%%;
end;%if flag_compute;

disp('returning'); return;

%%%%%%%%;
% Collect output. ;
%%%%%%%%;
verbose=1;
fp_label_A_ = fopen(sprintf('%s/dir_mat/str_CLabel_sub_.nsv',dir_trunk),'r');
str_label_A_ = textscan(fp_label_A_,'%s'); fclose(fp_label_A_);
n_u = numel(str_label_A_{1});
label_A_ = label_str_to_enum_0(str_label_A_{1});
u_label_A_ = unique(str_label_A_{1});
n_label_A = length(u_label_A_);
label_A_each__ = zeros(n_u,n_label_A);
for nlabel_A=0:n_label_A-1;
label_A_each__(:,1+nlabel_A) = zeros(n_u,1);
tmp_index_ = efind(strcmp(str_label_A_{1},u_label_A_{1+nlabel_A}));
label_A_each__(1+tmp_index_,1+nlabel_A) = 1;
end;%for nlabel_A=0:n_label_A-1;
u_label_A_enum_ = zeros(n_label_A,1); 
u_label_A_enum_(1:n_label_A-1) = cellfun(@str2num,u_label_A_(1:n_label_A-1));
u_label_A_enum_(end) = n_label_A;
[~,index_label_A_] = sort(u_label_A_enum_,'ascend'); index_label_A_ = index_label_A_ - 1;
%%%%%%%%;
str_X = 'E';
prefix_normalization = 'liXf';
tmp_dir_cluster = sprintf('%s/dir_%s_%s_cluster',dir_cluster,str_X,prefix_normalization);
X_n_ = E_liXf_; n_X_GENE = size(X_n_,2);
assert(size(X_n_,1)==n_u);
%%%%;
tmp_fname_mat = sprintf('%s/markergene_A__.mat',tmp_dir_cluster);
if (~exist(tmp_fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
%%%%%%%%;
markergene_auc_A__ = zeros(n_X_GENE,n_label_A);
markergene_auz_A__ = zeros(n_X_GENE,n_label_A); %<-- bates-equivalent z-score. ;
markergene_nlp_A__ = zeros(n_X_GENE,n_label_A); %<-- bates-equivalent z-score. ;
for nlabel_A=0:n_label_A-1;
tmp_index_0on_ = efind(label_A_each__(:,1+nlabel_A)==1);
tmp_index_off_ = efind(label_A_each__(:,1+nlabel_A)==0);
tmp_N = min(numel(tmp_index_0on_),numel(tmp_index_off_));
if (verbose); disp(sprintf(' %% nlabel_A %d/%d --> N %d/%d',nlabel_A,n_label_A,tmp_N,n_u)); end;
for nX_GENE=0:n_X_GENE-1;
tmp_auc = auc_0(X_n_(1+tmp_index_off_,1+nX_GENE),X_n_(1+tmp_index_0on_,1+nX_GENE));
tmp_auz = (tmp_auc - 0.5)*sqrt(12*tmp_N);
tmp_nlp = -logp_auc_0(tmp_auc,tmp_N);
markergene_auc_A__(1+nX_GENE,1+nlabel_A) = tmp_auc;
markergene_auz_A__(1+nX_GENE,1+nlabel_A) = tmp_auz;
markergene_nlp_A__(1+nX_GENE,1+nlabel_A) = tmp_nlp;
end;%for nX_GENE=0:n_X_GENE-1;
end;%for nlabel_A=0:n_label_A-1;
%%%%%%%%;
save(tmp_fname_mat ...
,'fp_label_A_','str_label_A_','n_u','label_A_' ...
,'u_label_A_','n_label_A','label_A_each__' ...
,'u_label_A_enum_','index_label_A_' ...
,'n_X_GENE' ...
,'markergene_auc_A__' ...
,'markergene_auz_A__' ...
,'markergene_nlp_A__' ...
);
end;%if (~exist(tmp_fname_mat,'file'));
load(tmp_fname_mat);
%%%%%%%%;
gamma = 0.01; n_shuffle = 64; p_set = 0.05; n_member_lob = 3;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('tmp_%s',str_xfix); 
dir_out = []; X_array_base_ = X_n_; X_array_r_ij_ = []; X_array_c_ij_ = []; flag_force_create = 0;
[ ...
 ZRimax_output_label_ ...
,ZRimax_lpFmax_label_ ...
] = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3( ...
 tmp_dir_cluster ...
,dir_out ...
,prefix_base ...
,X_array_base_ ...
,X_array_r_ij_ ...
,X_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,flag_force_create ...
);
%%%%%%%%;
% Clustering looks okay. ;
%%%%%%%%;
tmp_fname_mat = sprintf('%s/%s.mat',tmp_dir_cluster,str_xfix);
tmp_ = load(tmp_fname_mat);
n_Z_index = numel(tmp_.label_B__);
label_B__ = tmp_.label_B__;
%%%%%%%%;
tmp_fname_mat = sprintf('%s/%s_markergene_B___.mat',tmp_dir_cluster,str_xfix);
if (~exist(tmp_fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
%%%%%%%%;
lp0_A_each_vs_B__ = zeros(n_label_A,n_Z_index);
lpv_A_each_vs_B__ = zeros(n_label_A,n_Z_index);
for nZ_index=0:n_Z_index-1;
for nlabel_A=0:n_label_A-1;
lp0_A_each_vs_B__(1+nlabel_A,1+nZ_index) = label_to_label_enrichment_lP0(label_A_each__(:,1+nlabel_A),label_B__{1+nZ_index});
lpv_A_each_vs_B__(1+nlabel_A,1+nZ_index) = label_to_label_enrichment_2(label_A_each__(:,1+nlabel_A),label_B__{1+nZ_index});
end;%for nlabel_A=0:n_label_A-1;
end;%for nZ_index=0:n_Z_index-1;
%%%%%%%%;
u_label_B__ = cell(n_Z_index,1);
n_label_B_ = zeros(n_Z_index,1);
label_B_each___ = cell(n_Z_index,1);
for nZ_index=0:n_Z_index-1;
%%%%%%%%;
u_label_B__{1+nZ_index} = unique(label_B__{1+nZ_index});
n_label_B_(1+nZ_index) = length(u_label_B__{1+nZ_index});
label_B_each___{1+nZ_index} = zeros(n_u,n_label_B_(1+nZ_index));
for nlabel_B=0:n_label_B_(1+nZ_index)-1;
label_B_each___{1+nZ_index}(:,1+nlabel_B) = zeros(n_u,1);
tmp_index_ = efind(label_B__{1+nZ_index}==u_label_B__{1+nZ_index}(1+nlabel_B));
label_B_each___{1+nZ_index}(1+tmp_index_,1+nlabel_B) = 1;
end;%for nlabel_B=0:n_label_B_(1+nZ_index)-1;
%%%%%%%%;
nlpv_each__ = zeros(n_label_A,n_label_B_(1+nZ_index));
for nlabel_A=0:n_label_A-1;
for nlabel_B=0:n_label_B_(1+nZ_index)-1;
nlpv_each__(1+nlabel_A,1+nlabel_B) = -label_pair_enrichment(label_A_each__(:,1+nlabel_A),label_B_each___{1+nZ_index}(:,1+nlabel_B));
end;%for nlabel_B=0:n_label_B_(1+nZ_index)-1;
end;%for nlabel_A=0:n_label_A-1;
%%%%%%%%;
markergene_auc_B___{1+nZ_index} = zeros(n_X_GENE,n_label_B_(1+nZ_index));
markergene_auz_B___{1+nZ_index} = zeros(n_X_GENE,n_label_B_(1+nZ_index)); %<-- bates-equivalent z-score. ;
markergene_nlp_B___{1+nZ_index} = zeros(n_X_GENE,n_label_B_(1+nZ_index)); %<-- bates-equivalent z-score. ;
for nlabel_B=0:n_label_B_(1+nZ_index)-1;
tmp_index_0on_ = efind(label_B_each___{1+nZ_index}(:,1+nlabel_B)==1);
tmp_index_off_ = efind(label_B_each___{1+nZ_index}(:,1+nlabel_B)==0);
tmp_N = min(numel(tmp_index_0on_),numel(tmp_index_off_));
if (verbose); disp(sprintf(' %% nlabel_B %d/%d --> N %d/%d',nlabel_B,n_label_B_(1+nZ_index),tmp_N,n_u)); end;
for nX_GENE=0:n_X_GENE-1;
tmp_auc = auc_0(X_n_(1+tmp_index_off_,1+nX_GENE),X_n_(1+tmp_index_0on_,1+nX_GENE));
tmp_auz = (tmp_auc - 0.5)*sqrt(12*tmp_N);
tmp_nlp = -logp_auc_0(tmp_auc,tmp_N);
markergene_auc_B___{1+nZ_index}(1+nX_GENE,1+nlabel_B) = tmp_auc;
markergene_auz_B___{1+nZ_index}(1+nX_GENE,1+nlabel_B) = tmp_auz;
markergene_nlp_B___{1+nZ_index}(1+nX_GENE,1+nlabel_B) = tmp_nlp;
end;%for nX_GENE=0:n_X_GENE-1;
end;%for nlabel_B=0:n_label_B_(1+nZ_index)-1;
%%%%%%%%;
end;%for nZ_index=0:n_Z_index-1;
%%%%%%%%;
save(tmp_fname_mat ...
,'fp_label_A_','str_label_A_','n_u','label_A_' ...
,'u_label_A_','n_label_A','label_A_each__' ...
,'u_label_A_enum_','index_label_A_' ...
,'n_X_GENE' ...
,'n_Z_index','label_B__' ...
,'u_label_B__','n_label_B_','label_B_each___','nlpv_each__' ...
,'markergene_auc_B___' ...
,'markergene_auz_B___' ...
,'markergene_nlp_B___' ...
,'lp0_A_each_vs_B__' ...
,'lpv_A_each_vs_B__' ...
);
end;%if (~exist(tmp_fname_mat,'file'));
load(tmp_fname_mat);
%%%%%%%%;

for nZ_index=0:n_Z_index-1;
%%%%%%%%;
tmp_dir_jpg = sprintf('%s/dir_jpg',tmp_dir_cluster);
if (~exist(tmp_dir_jpg,'dir')); sprintf(' %% mkdir %s',tmp_dir_jpg); mkdir(tmp_dir_jpg); end;
fname_fig = sprintf('%s/dir_jpg/%s_z%d_lpv_FIGA',tmp_dir_cluster,str_xfix,1+nZ_index);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf('%% %s not found, creating',fname_fig));
%%%%%%%%;
[tmp_lP0,tmp_cap_] = label_to_label_enrichment_lP0(label_A_,label_B__{1+nZ_index});
%%%%%%%%;
tmp_index_label_B_ = 0:n_label_B_(1+nZ_index)-1;
index_label_B_ = zeros(n_label_B_(1+nZ_index),1);
for nlabel_A=0:min(n_label_A,n_label_B_(1+nZ_index))-1;
[~,tmp_index] = max(nlpv_each__(1+index_label_A_(1+nlabel_A),1+tmp_index_label_B_)); tmp_index = tmp_index - 1;
index_label_B_(1+nlabel_A) = tmp_index_label_B_(1+tmp_index);
tmp_index_label_B_ = setdiff(tmp_index_label_B_,index_label_B_(1+nlabel_A));
end;%for nlabel_A=0:min(n_label_A,n_label_B_(1+nZ_index))-1;
%%%%%%%%;
figure(1);clf;
figbeach;
subplot(1,2,1);
imagesc(tmp_cap_(1+index_label_A_,1+index_label_B_),[0,100]); colorbar;
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gca,'XTick',1:n_label_B_(1+nZ_index),'XTickLabel',u_label_B__{1+nZ_index}(1+index_label_B_)); xtickangle(90);
axis image; title('cap');
subplot(1,2,2);
imagesc(nlpv_each__(1+index_label_A_,1+index_label_B_),[0,9]); colorbar;
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gca,'XTick',1:n_label_B_(1+nZ_index),'XTickLabel',u_label_B__{1+nZ_index}(1+index_label_B_)); xtickangle(90);
colorbar;
axis image; title('nlpv');
sgtitle(sprintf('%s_z%d',str_xfix,1+nZ_index),'Interpreter','none');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%for nZ_index=0:n_Z_index-1;

disp('returning'); return; 
