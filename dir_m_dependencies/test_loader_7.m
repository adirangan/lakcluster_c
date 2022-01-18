% testing a loader for table data ;
clear;
setup;
dir_trunk = sprintf('/data/rangan/dir_bcc/dir_jamison');
dir_data = sprintf('%s/data_summary_20190730',dir_trunk);
%%%%%%%%;
flag_load=0;
if flag_load;
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
end;%if flag_load;
%%%%%%%%;

%%%%%%%%;
load(sprintf('%s/EXON.Gx_50421.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'E_');
load(sprintf('%s/INTRON_INTERGENIC.Gx_25998.Cx_50.brain_table_merge_output.CPM.NotCtrl.NotEmpty.mat',dir_data),'I_');
load(sprintf('%s/20161026_covariate_table.format.mat',dir_data),'C_');
%%%%%%%%;

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
flag_plot=1;
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

%%%%%%%%;
E_ID_ = E_{:,1};
I_ID_ = E_{:,1};
C_ID_ = C_{:,1}; 
u_ID_ = intersectall({E_ID_,I_ID_,C_ID_});
n_u = length(u_ID_);
[~,~,E_to_u_ID_] = intersect(u_ID_,E_ID_,'stable');
[~,~,I_to_u_ID_] = intersect(u_ID_,I_ID_,'stable');
[~,~,C_to_u_ID_] = intersect(u_ID_,C_ID_,'stable');
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

%%%%%%%%;
% Extracting unique cluster labels. ;
%%%%%%%%;
u_CLabel_ = unique(C_.Cluster_ID_20161007(C_to_u_ID_));
n_CLabel = length(u_CLabel_);
n_CLabel_ = zeros(n_CLabel,1);
for nCLabel = 1:n_CLabel;
n_CLabel_(nCLabel) = length(find(strcmp(C_.Cluster_ID_20161007(C_to_u_ID_),u_CLabel_(nCLabel))));
end;%for nCLabel = 1:n_CLabel;
flag_plot=1;
if flag_plot;
bar(1:n_CLabel,n_CLabel_);
set(gca,'XTick',1:n_CLabel,'XTickLabel',u_CLabel_); xtickangle(90);
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
E_col_val_ = 52:size(E_,2)-1; %<-- gene list begins after 'Fail Confidence Sum' ;
E_val_ = decostand_total_0(E_{E_to_u_ID_,E_col_val_},'col');
n_E_GENE = size(E_val_,2);
I_col_val_ = 52:size(I_,2)-1;
I_val_ = decostand_total_0(I_{I_to_u_ID_,I_col_val_},'col');
n_I_GENE = size(I_val_,2);  %<-- gene list begins after 'Fail Confidence Sum' ;
%assert(size(I_val_,1)==size(E_val_,1)); %<-- fewer samples in E_. ;
E_rank_ = rank_normalize_0(E_val_,'row');;
I_rank_ = rank_normalize_0(I_val_,'row');;

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
C_VariableName_LC_ = C_VariableName_;
for nCCOV=1:n_CCOV;
if C_LC_(nCCOV) == -1; tmp_suffix = '[ - ]'; end;
if C_LC_(nCCOV) ==  0; tmp_suffix = '[    ]'; end;
if C_LC_(nCCOV) == +1; tmp_suffix = '[+]'; end;
C_VariableName_LC_{nCCOV} = sprintf('%s_%s',C_VariableName_{nCCOV},tmp_suffix);
end;%for nCCOV=1:n_CCOV;
clear tmp_suffix ;
C_val_ = [B_val_ , D_val_];
C_rank_ = [B_rank_ , D_rank_];
C_LC_ = [B_LC_ ; D_LC_];

%%%%%%%%;
% Now measure correlation across samples of E_ and C_. ;
%%%%%%%%;
tmp_rank_avg_ = mean(C_rank_,1); tmp_rank_std_ = max(1,std(C_rank_,1,1));
tmp_C_ = (C_rank_ - repmat(tmp_rank_avg_,n_u,1))./repmat(tmp_rank_std_,n_u,1);
tmp_rank_avg_ = mean(E_rank_,1); tmp_rank_std_ = max(1,std(E_rank_,1,1));
tmp_E_ = (E_rank_ - repmat(tmp_rank_avg_,n_u,1))./repmat(tmp_rank_std_,n_u,1);
clear tmp_rank_avg_ tmp_rank_std_;
CtEn_rank_ = transpose(tmp_C_)*tmp_E_ / n_u;
[tmp_U_,tmp_S_,tmp_V_] = svds(CtEn_rank_,1); [~,CtEn_rank_U_ij_] = sort(tmp_U_,'descend'); [~,CtEn_rank_V_ij_] = sort(tmp_V_,'descend');
CtEn_rank_ori_ = CtEn_rank_(CtEn_rank_U_ij_,CtEn_rank_V_ij_);
clear CtEn_rank_ tmp_U_ tmp_S_ tmp_V_ ;
flag_disp=1;
if flag_disp;
colormap(colormap_beach()); 
imagesc(transpose(CtEn_rank_ori_),[-1,+1]); 
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(CtEn_rank_U_ij_),'TickLabelInterpreter','none'); xtickangle(90);
set(gca,'FontSize',7);
xlabel('covariates'); ylabel('genes'); title('Correlation between Genes and Covariates');
set(gcf,'Position',1+[0,0,1024*2,1024]);
colorbar;
fname_base = sprintf('%s/dir_jpg/CtEn_ori',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_disp;
n_iteration = 32;
CtEn_rank_avg_ = zeros(n_CCOV,n_E_GENE);
CtEn_rank_std_ = zeros(n_CCOV,n_E_GENE);
for niteration=1:n_iteration;
if (mod(niteration,10)==0); disp(sprintf(' %% niteration %d/%d',niteration,n_iteration)); end;
[tmp_Q_,~] = qr(randn(n_u));
tmp_CtEn_rank_ = transpose(tmp_C_)*tmp_Q_*tmp_E_ / n_u; clear tmp_Q_;
[tmp_U_,tmp_S_,tmp_V_] = svds(tmp_CtEn_rank_,1); [~,tmp_CtEn_rank_U_ij_] = sort(tmp_U_,'descend'); [~,tmp_CtEn_rank_V_ij_] = sort(tmp_V_,'descend');
CtEn_rank_avg_ = CtEn_rank_avg_ + tmp_CtEn_rank_(tmp_CtEn_rank_U_ij_,tmp_CtEn_rank_V_ij_);
CtEn_rank_std_ = CtEn_rank_std_ + tmp_CtEn_rank_(tmp_CtEn_rank_U_ij_,tmp_CtEn_rank_V_ij_).^2;
clear tmp_CtEn_rank_ tmp_U_ tmp_S_ tmp_V_ tmp_CtEn_rank_U_ij_ tmp_CtEn_rank_V_ij_ ;
end;%for niteration=1:n_iteration;
CtEn_rank_avg_ = CtEn_rank_avg_/n_iteration;
CtEn_rank_std_ = sqrt(CtEn_rank_std_/n_iteration - CtEn_rank_avg_.^2);
CtEn_rank_Z_ = ( CtEn_rank_ori_ - CtEn_rank_avg_ ) ./ CtEn_rank_std_ ;
%CtEn_rank_p_ = 0.5 * erfc( CtEn_rank_Z_ / sqrt(2) ) ;
CtEn_rank_p_ = 1.0 * erfc( abs(CtEn_rank_Z_) / sqrt(2) ) ; %<-- two sided. ;
clear tmp_C_ tmp_E_;
flag_disp=1;
if flag_disp;
colormap(colormap_beach()); 
imagesc(transpose(-log10(CtEn_rank_p_)),[0,+15]); 
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(CtEn_rank_U_ij_),'TickLabelInterpreter','none'); xtickangle(90);
set(gca,'FontSize',7);
xlabel('covariates'); ylabel('genes'); title('log10(p-value) of Correlation between Genes and Covariates');
set(gcf,'Position',1+[0,0,1024*2,1024]);
colorbar;
fname_base = sprintf('%s/dir_jpg/CtEn_p',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_disp;

%%%%%%%%;
% Now measure correlation across samples of I_ and C_. ;
%%%%%%%%;
tmp_rank_avg_ = mean(C_rank_,1); tmp_rank_std_ = max(1,std(C_rank_,1,1));
tmp_C_ = (C_rank_ - repmat(tmp_rank_avg_,n_u,1))./repmat(tmp_rank_std_,n_u,1);
tmp_rank_avg_ = mean(I_rank_,1); tmp_rank_std_ = max(1,std(I_rank_,1,1));
tmp_I_ = (I_rank_ - repmat(tmp_rank_avg_,n_u,1))./repmat(tmp_rank_std_,n_u,1);
clear tmp_rank_avg_ tmp_rank_std_;
CtIn_rank_ = transpose(tmp_C_)*tmp_I_ / n_u;
[tmp_U_,tmp_S_,tmp_V_] = svds(CtIn_rank_,1); [~,CtIn_rank_U_ij_] = sort(tmp_U_,'descend'); [~,CtIn_rank_V_ij_] = sort(tmp_V_,'descend');
CtIn_rank_ori_ = CtIn_rank_(CtIn_rank_U_ij_,CtIn_rank_V_ij_);
clear CtIn_rank_ tmp_U_ tmp_S_ tmp_V_ ;
flag_disp=1;
if flag_disp;
colormap(colormap_beach()); 
imagesc(transpose(CtIn_rank_ori_),[-1,+1]); 
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(CtIn_rank_U_ij_),'TickLabelInterpreter','none'); xtickangle(90);
set(gca,'FontSize',7);
xlabel('covariates'); ylabel('genes'); title('Correlation between Genes and Covariates');
set(gcf,'Position',1+[0,0,1024*2,1024]);
colorbar;
fname_base = sprintf('%s/dir_jpg/CtIn_ori',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_disp;
n_iteration = 32;
CtIn_rank_avg_ = zeros(n_CCOV,n_I_GENE);
CtIn_rank_std_ = zeros(n_CCOV,n_I_GENE);
for niteration=1:n_iteration;
if (mod(niteration,10)==0); disp(sprintf(' %% niteration %d/%d',niteration,n_iteration)); end;
[tmp_Q_,~] = qr(randn(n_u));
tmp_CtIn_rank_ = transpose(tmp_C_)*tmp_Q_*tmp_I_ / n_u; clear tmp_Q_;
[tmp_U_,tmp_S_,tmp_V_] = svds(tmp_CtIn_rank_,1); [~,tmp_CtIn_rank_U_ij_] = sort(tmp_U_,'descend'); [~,tmp_CtIn_rank_V_ij_] = sort(tmp_V_,'descend');
CtIn_rank_avg_ = CtIn_rank_avg_ + tmp_CtIn_rank_(tmp_CtIn_rank_U_ij_,tmp_CtIn_rank_V_ij_);
CtIn_rank_std_ = CtIn_rank_std_ + tmp_CtIn_rank_(tmp_CtIn_rank_U_ij_,tmp_CtIn_rank_V_ij_).^2;
clear tmp_CtIn_rank_ tmp_U_ tmp_S_ tmp_V_ tmp_CtIn_rank_U_ij_ tmp_CtIn_rank_V_ij_ ;
end;%for niteration=1:n_iteration;
CtIn_rank_avg_ = CtIn_rank_avg_/n_iteration;
CtIn_rank_std_ = sqrt(CtIn_rank_std_/n_iteration - CtIn_rank_avg_.^2);
CtIn_rank_Z_ = ( CtIn_rank_ori_ - CtIn_rank_avg_ ) ./ CtIn_rank_std_ ;
%CtIn_rank_p_ = 0.5 * erfc( CtIn_rank_Z_ / sqrt(2) ) ;
CtIn_rank_p_ = 1.0 * erfc( abs(CtIn_rank_Z_) / sqrt(2) ) ; %<-- two sided. ;
clear tmp_C_ tmp_I_;
flag_disp=1;
if flag_disp;
colormap(colormap_beach()); 
imagesc(transpose(-log10(CtIn_rank_p_)),[0,+15]); 
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(CtIn_rank_U_ij_),'TickLabelInterpreter','none'); xtickangle(90);
set(gca,'FontSize',7);
xlabel('covariates'); ylabel('genes'); title('log10(p-value) of Correlation between Genes and Covariates');
set(gcf,'Position',1+[0,0,1024*2,1024]);
colorbar;
fname_base = sprintf('%s/dir_jpg/CtIn_p',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_disp;

flag_disp=1;
if flag_disp;
n_bin = 32;
%%%%%%%%;
subplot(2,4,1); hold on;
tmp_ng_ = 1:n_E_GENE;
%plot(tmp_ng_,log10(1+E_length_(E_col_val_(CtEn_rank_V_ij_))),'r.');
tmp_ij_ = find(E_length_(E_col_val_(CtEn_rank_V_ij_))>0);
plot(tmp_ng_(tmp_ij_),log10(1+E_length_(E_col_val_(CtEn_rank_V_ij_(tmp_ij_)))),'ro');
xlabel('gene-index (matching CtEn_p order)'); xlim([0,n_E_GENE]); tmp_yl_ = ylim();
ylabel('log10(gene length)');
title('E gene length vs gene ordering');
subplot(2,4,2); hold on;
tmp_ng_ = 1:n_E_GENE;
tmp_ij_ = find(E_length_(E_col_val_(CtEn_rank_V_ij_))>0);
colormap(colormap_beach());
imagesc(log10(1+hist2d_0(tmp_ng_(tmp_ij_),log10(1+E_length_(E_col_val_(CtEn_rank_V_ij_(tmp_ij_)))),n_bin,n_bin,[1,n_E_GENE],tmp_yl_)));
colorbar();
xlabel('gene-index (matching CtEn_p order)'); ylabel('log10(gene length)');
xlim([1,n_bin]); ylim([1,n_bin]);
set(gca,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
title('E gene length vs gene ordering');
%%%%%%%%;
subplot(2,4,3); hold on;
tmp_ng_ = 1:n_I_GENE;
%plot(tmp_ng_,log10(1+I_length_(I_col_val_(CtIn_rank_V_ij_))),'b.');
tmp_ij_ = find(I_length_(I_col_val_(CtIn_rank_V_ij_))>0);
plot(tmp_ng_(tmp_ij_),log10(1+I_length_(I_col_val_(CtIn_rank_V_ij_(tmp_ij_)))),'bo');
xlabel('gene-index (matching CtIn_p order)'); xlim([0,n_I_GENE]);
ylabel('log10(gene length)');
title('I gene length vs gene ordering');
subplot(2,4,4); hold on;
tmp_ng_ = 1:n_I_GENE;
tmp_ij_ = find(I_length_(I_col_val_(CtIn_rank_V_ij_))>0);
colormap(colormap_beach()); 
imagesc(log10(1+hist2d_0(tmp_ng_(tmp_ij_),log10(1+I_length_(I_col_val_(CtIn_rank_V_ij_(tmp_ij_)))),n_bin,n_bin,[1,n_I_GENE],tmp_yl_)));
colorbar();
xlabel('gene-index (matching CtIn_p order)'); ylabel('log10(gene length)');
xlim([1,n_bin]); ylim([1,n_bin]);
set(gca,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
title('I gene length vs gene ordering');
%%%%%%%%;
subplot(2,4,5); hold on;
tmp_ng_ = 1:n_E_GENE;
tmp_ij_ = find(E_length_(E_col_val_(CtEn_rank_V_ij_))>0);
tmp_X_ = corr(transpose(log10(sum(abs(CtEn_rank_ori_(:,tmp_ij_)).^2,1))),log10(E_length_(E_col_val_(CtEn_rank_V_ij_(tmp_ij_)))));
plot(log10(sum(abs(CtEn_rank_ori_(:,tmp_ij_)).^2,1)),log10(E_length_(E_col_val_(CtEn_rank_V_ij_(tmp_ij_)))),'ro');
xlabel('log10(sum squared correlation)'); ylabel('log10(gene length)');
tmp_yl_ = ylim(); tmp_xl_ = xlim();
title(sprintf('avg-corr^2 vs E gene-length (C %0.6f)',tmp_X_));
subplot(2,4,6); hold on;
colormap(colormap_beach());
imagesc(log10(1+hist2d_0(log10(sum(abs(CtEn_rank_ori_(:,tmp_ij_)).^2,1)),log10(E_length_(E_col_val_(CtEn_rank_V_ij_(tmp_ij_)))),n_bin,n_bin,tmp_xl_,tmp_yl_)));
colorbar();
colorbar();
xlabel('log10(sum squared correlation)'); ylabel('log10(gene length)');
xlim([1,n_bin]); ylim([1,n_bin]);
set(gca,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
title(sprintf('avg-corr^2 vs E gene-length (C %0.6f)',tmp_X_));
%%%%%%%%;
subplot(2,4,7); hold on;
tmp_ng_ = 1:n_I_GENE;
tmp_ij_ = find(I_length_(I_col_val_(CtIn_rank_V_ij_))>0);
tmp_X_ = corr(transpose(log10(sum(abs(CtIn_rank_ori_(:,tmp_ij_)).^2,1))),log10(I_length_(I_col_val_(CtIn_rank_V_ij_(tmp_ij_)))));
plot(log10(sum(abs(CtIn_rank_ori_(:,tmp_ij_)).^2,1)),log10(I_length_(I_col_val_(CtIn_rank_V_ij_(tmp_ij_)))),'bo');
xlabel('log10(sum squared correlation)'); ylabel('log10(gene length)');
tmp_yl_ = ylim(); tmp_xl_ = xlim();
title(sprintf('avg-corr^2 vs I gene-length (C %0.6f)',tmp_X_));
subplot(2,4,8); hold on;
colormap(colormap_beach());
imagesc(log10(1+hist2d_0(log10(sum(abs(CtIn_rank_ori_(:,tmp_ij_)).^2,1)),log10(I_length_(I_col_val_(CtIn_rank_V_ij_(tmp_ij_)))),n_bin,n_bin,tmp_xl_,tmp_yl_)));
colorbar();
colorbar();
xlabel('log10(sum squared correlation)'); ylabel('log10(gene length)');
xlim([1,n_bin]); ylim([1,n_bin]);
set(gca,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
title(sprintf('avg-corr^2 vs I gene-length (C %0.6f)',tmp_X_));
%%%%%%%%;
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_base = sprintf('%s/dir_jpg/CtXn_vs_length',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
clear tmp_xl_ tmp_yl_ tmp_X_ tmp_ij_ tmp_ng_ ;
end;%if flag_disp;
%%%%%%%%;
% Note that there is a modest correlation between length and covariate-correlation. ;
%%%%%%%%;

%%%%%%%%;
% Here we try and estimate the genes E_ or I_ using the covariates C_. ;
% For now we use all of C_ to estimate E_ and I_. ;
% Eventually, we might want to use only the C_ which are not associated with quality. ;
%%%%%%%%;
n_rank = 3+3;
%%%%%%%%;
% Can we estimate the link E_rank_ = C_rank_ * zeta_E_  using a low-rank zeta_E_ ?
%%%%%%%%;
zeta_E_vn_ = zeros(n_E_GENE,n_rank);
zeta_E_un_ = zeros(n_CCOV,n_rank);
tmp_Bn_ = pinv(C_rank_);
tmp_Yn_ = E_rank_;
for nrank=1:n_rank;
[tmp_un_,tmp_vn_] = rrr_0(C_rank_,tmp_Yn_,tmp_Bn_);
zeta_E_un_(:,nrank) = tmp_un_;
zeta_E_vn_(:,nrank) = tmp_vn_;
tmp_Yn_ = tmp_Yn_ - (C_rank_*zeta_E_un_(:,nrank))*transpose(zeta_E_vn_(:,nrank));
clear tmp_un_ tmp_vn_;
end;%for nrank=1:n_rank;
clear tmp_Yn_ tmp_Bn_;
%%%%%%%%
% Can we estimate the link I_rank_ = C_rank_ * zeta_I_ using a low-rank zeta_I_ ?
%%%%%%%%
zeta_I_vn_ = zeros(n_I_GENE,n_rank);
zeta_I_un_ = zeros(n_CCOV,n_rank);
tmp_Bn_ = pinv(C_rank_);
tmp_Yn_ = I_rank_;
for nrank=1:n_rank;
[tmp_un_,tmp_vn_] = rrr_0(C_rank_,tmp_Yn_,tmp_Bn_);
zeta_I_un_(:,nrank) = tmp_un_;
zeta_I_vn_(:,nrank) = tmp_vn_;
tmp_Yn_ = tmp_Yn_ - (C_rank_*zeta_I_un_(:,nrank))*transpose(zeta_I_vn_(:,nrank));
clear tmp_un_ tmp_vn_;
end;%for nrank=1:n_rank;
clear tmp_Yn_ tmp_Bn_;
%%%%%%%%;
% Appears to be rank 3 or so ;
%%%%%%%%;
flag_disp=1;
if flag_disp;
subplot(1,2,1);
plot(1:n_rank,sqrt(sum((C_rank_*zeta_E_un_).^2,1)),'ro');
xlabel('pc #'); ylabel('l2-norm');
title('numerical rank of zeta_E_','Interpreter','none');
subplot(1,2,2);
plot(1:n_rank,sqrt(sum((C_rank_*zeta_I_un_).^2,1)),'bo');
xlabel('pc #'); ylabel('l2-norm');
title('numerical rank of zeta_I_','Interpreter','none');
set(gcf,'Position',1+[0,0,512*2,512]);
fname_base = sprintf('%s/dir_jpg/zeta_X_',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_disp;
%%%%%%%%;
n_zeta_E = 3;
n_zeta_I = 3;

flag_disp=0;
if flag_disp;
clf;
subplot(1,2,1); tmp_ = (C_rank_*zeta_E_un_)*transpose(zeta_E_vn_); hist(tmp_(:),32);
tmp_n = length(find(tmp_(:)<1 | tmp_(:)>n_E_GENE));
disp(sprintf(' %% zeta_E_: %d/%d = %0.2f out of bounds',tmp_n,length(tmp_(:)),tmp_n/length(tmp_(:))))
subplot(1,2,2); tmp_ = (C_rank_*zeta_I_un_)*transpose(zeta_I_vn_); hist(tmp_(:),32);
tmp_n = length(find(tmp_(:)<1 | tmp_(:)>n_I_GENE));
disp(sprintf(' %% zeta_I_: %d/%d = %0.2f out of bounds',tmp_n,length(tmp_(:)),tmp_n/length(tmp_(:))))
clear tmp_ tmp_n;
end;%if flag_disp;
%%%%%%%%;
% only 0.00 out of bounds. ;
%%%%%%%%;

%%%%%%%%;
% Here we try and estimate the covariates C_ using the genes E_ or I_. ;
% For now we estimate all of C_. ;
% Eventually, we might want to only estimate the C_ which are not associated with a quality. ;
%%%%%%%%;
n_rank = 6+3;
%%%%%%%%;
% Can we estimate the link E_rank_ * beta_E_ = C_rank_ using a low-rank beta_E_ ?
%%%%%%%%;
beta_E_un_ = zeros(n_E_GENE,n_rank);
beta_E_vn_ = zeros(n_CCOV,n_rank);
tmp_Bn_ = pinv(E_rank_);
tmp_Yn_ = C_rank_;
for nrank=1:n_rank;
[tmp_un_,tmp_vn_] = rrr_0(E_rank_,tmp_Yn_,tmp_Bn_);
beta_E_un_(:,nrank) = tmp_un_;
beta_E_vn_(:,nrank) = tmp_vn_;
tmp_Yn_ = tmp_Yn_ - (E_rank_*beta_E_un_(:,nrank))*transpose(beta_E_vn_(:,nrank));
clear tmp_un_ tmp_vn_;
end;%for nrank=1:n_rank;
clear tmp_Yn_ tmp_Bn_;
%%%%%%%%
% Can we estimate the link I_rank_ * beta_I_ = C_rank_ using a low-rank beta_I_ ?
%%%%%%%%
beta_I_un_ = zeros(n_I_GENE,n_rank);
beta_I_vn_ = zeros(n_CCOV,n_rank);
tmp_Bn_ = pinv(I_rank_);
tmp_Yn_ = C_rank_;
for nrank=1:n_rank;
[tmp_un_,tmp_vn_] = rrr_0(I_rank_,tmp_Yn_,tmp_Bn_);
beta_I_un_(:,nrank) = tmp_un_;
beta_I_vn_(:,nrank) = tmp_vn_;
tmp_Yn_ = tmp_Yn_ - (I_rank_*beta_I_un_(:,nrank))*transpose(beta_I_vn_(:,nrank));
clear tmp_un_ tmp_vn_;
end;%for nrank=1:n_rank;
clear tmp_Yn_ tmp_Bn_;
%%%%%%%%;
% Appears to be rank 6 or so ;
%%%%%%%%;
flag_disp=1;
if flag_disp;
subplot(1,2,1);
plot(1:n_rank,sqrt(sum((E_rank_*beta_E_un_).^2,1)),'ro');
xlabel('pc #'); ylabel('l2-norm');
title('numerical rank of beta_E_','Interpreter','none');
subplot(1,2,2);
plot(1:n_rank,sqrt(sum((I_rank_*beta_I_un_).^2,1)),'bo');
xlabel('pc #'); ylabel('l2-norm');
title('numerical rank of beta_I_','Interpreter','none');
set(gcf,'Position',1+[0,0,512*2,512]);
fname_base = sprintf('%s/dir_jpg/beta_X_',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_disp;
%%%%%%%%;
n_beta_E = 6;
n_beta_I = 6;

flag_disp=0;
if flag_disp;
clf;
subplot(2,2,1); imagesc(C_rank_); title('Cn');
subplot(2,2,2); imagesc((I_rank_*beta_I_un_)*transpose(beta_I_vn_)); title('In*betan');
subplot(2,2,3); imagesc(transpose(I_rank_)*C_rank_); title('It*Cn');
subplot(2,2,4); imagesc(transpose(I_rank_)*((I_rank_*beta_I_un_)*transpose(beta_I_vn_))); title('It*In*betan');
end;%if flag_disp;
%%%%%%%%;
% is a logistic model necessary, or does linear model fit well? 
%%%%%%%%;
flag_disp=0;
if flag_disp;
clf;
subplot(1,2,1); tmp_ = (E_rank_*beta_E_un_)*transpose(beta_E_vn_); hist(tmp_(:),32);
tmp_n = length(find(tmp_(:)<1 | tmp_(:)>n_u));
disp(sprintf(' %% beta_E_: %d/%d = %0.2f out of bounds',tmp_n,length(tmp_(:)),tmp_n/length(tmp_(:))))
subplot(1,2,2); tmp_ = (I_rank_*beta_I_un_)*transpose(beta_I_vn_); hist(tmp_(:),32);
tmp_n = length(find(tmp_(:)<1 | tmp_(:)>n_u));
disp(sprintf(' %% beta_I_: %d/%d = %0.2f out of bounds',tmp_n,length(tmp_(:)),tmp_n/length(tmp_(:))))
clear tmp_ tmp_n;
end;%if flag_disp;
%%%%%%%%;
% only 0.03 out of bounds. ;
%%%%%%%%;

%%%%%%%%;
% Now measure rank-differences. ;
% Use label J_ to refer to E_rank_ * beta_E_un_ ;
% Use label K_ to refer to I_rank_ * beta_I_un_ ;
% Use lable F_ to refer to C_rank_ * zeta_E_un_ ;
% Use lable H_ to refer to C_rank_ * zeta_I_un_ ;
%%%%%%%%;
C_Label_ = C_.Cluster_ID_20161007(C_to_u_ID_);
u_CLabel_ = unique(C_Label_); n_CLabel = length(u_CLabel_);
E_aucX__ = zeros(n_CLabel,n_CLabel,n_E_GENE);
I_aucX__ = zeros(n_CLabel,n_CLabel,n_I_GENE);
C_aucX__ = zeros(n_CLabel,n_CLabel,n_CCOV);
J_aucX__ = zeros(n_CLabel,n_CLabel,n_CCOV);
K_aucX__ = zeros(n_CLabel,n_CLabel,n_CCOV);
E_rank_avg_ = zeros(n_CLabel,n_E_GENE);
I_rank_avg_ = zeros(n_CLabel,n_I_GENE);
C_rank_avg_ = zeros(n_CLabel,n_CCOV);
J_rank_avg_ = zeros(n_CLabel,n_CCOV);
K_rank_avg_ = zeros(n_CLabel,n_CCOV);
%%%%%%%%;
for nA = 1:n_CLabel;
disp(sprintf(' %% Rank average for Label %d (%s):',nA,u_CLabel_{nA}));
tmp_Label_A_ij_ = find(strcmp(C_Label_,u_CLabel_{nA})); tmp_n_A = length(tmp_Label_A_ij_);
E_rank_avg_(nA,:) = mean(E_rank_(tmp_Label_A_ij_,:),1);
I_rank_avg_(nA,:) = mean(I_rank_(tmp_Label_A_ij_,:),1);
C_rank_avg_(nA,:) = mean(C_rank_(tmp_Label_A_ij_,:),1);
J_rank_avg_(nA,:) = mean((E_rank_(tmp_Label_A_ij_,:)*beta_E_un_(:,1:n_beta_E))*transpose(beta_E_vn_(:,1:n_beta_E)),1);
K_rank_avg_(nA,:) = mean((I_rank_(tmp_Label_A_ij_,:)*beta_I_un_(:,1:n_beta_I))*transpose(beta_I_vn_(:,1:n_beta_I)),1);
end;%for nA = 1:n_CLabel;
%%%%%%%%;
for nA = 1:n_CLabel;
for nB = nA+1:n_CLabel;
disp(sprintf(' %% Comparing Label %d (%s) and %d (%s):',nA,u_CLabel_{nA},nB,u_CLabel_{nB}));
E_aucX__(nA,nB,:) = E_rank_avg_(nA,:) - E_rank_avg_(nB,:);
E_aucX__(nB,nA,:) = -E_aucX__(nA,nB,:);
I_aucX__(nA,nB,:) = I_rank_avg_(nA,:) - I_rank_avg_(nB,:);
I_aucX__(nB,nA,:) = -I_aucX__(nA,nB,:);
C_aucX__(nA,nB,:) = C_rank_avg_(nA,:) - C_rank_avg_(nB,:);
C_aucX__(nB,nA,:) = -C_aucX__(nA,nB,:);
J_aucX__(nA,nB,:) = [C_rank_avg_(nA,:) - C_rank_avg_(nB,:)] - [J_rank_avg_(nA,:) - J_rank_avg_(nB,:)];
J_aucX__(nB,nA,:) = -J_aucX__(nA,nB,:);
K_aucX__(nA,nB,:) = [C_rank_avg_(nA,:) - C_rank_avg_(nB,:)] - [K_rank_avg_(nA,:) - K_rank_avg_(nB,:)];
K_aucX__(nB,nA,:) = -K_aucX__(nA,nB,:);
end;%for nB = nA+1:n_CLabel;
end;%for nA = 1:n_CLabel;
clear tmp_Label_A_ij_ tmp_n_A ;

flag_test=0;
if flag_test;
%%%%%%%;
% Now we can test a pair of labels: ;
%%%%%%%%;
nA=find(strcmp(u_CLabel_,'28')); nB=find(strcmp(u_CLabel_,'33'));
disp(sprintf(' %% Comparing Label %d (%s) and %d (%s):',nA,u_CLabel_{nA},nB,u_CLabel_{nB}));
tmp_Label_A_ij_ = find(strcmp(C_Label_,u_CLabel_{nA})); tmp_n_A = length(tmp_Label_A_ij_);
tmp_Label_B_ij_ = find(strcmp(C_Label_,u_CLabel_{nB})); tmp_n_B = length(tmp_Label_B_ij_);
%%%%%%%%;
% Now sort genes by aucX and plot heatmap: ;
%%%%%%%%;
[~,tmp_E_ij_] = sort(E_aucX__(nA,nB,:));
tmp_AB_ = [E_rank_(tmp_Label_A_ij_,tmp_E_ij_) ; E_rank_(tmp_Label_B_ij_,tmp_E_ij_)];
tmp_AB_ = decostand_total_0(tmp_AB_,'row');
clim_ = mean(tmp_AB_(:)) + 1.5*std(tmp_AB_(:))*[-1,+1];
figure(1); clf;
subplot(5,1,[1:2]); imagesc(tmp_AB_(0*tmp_n_A + (1:tmp_n_A),:),clim_); ylabel(sprintf('label %s',u_CLabel_{nA}));
subplot(5,1,[3:4]); imagesc(tmp_AB_(1*tmp_n_A + (1:tmp_n_B),:),clim_); ylabel(sprintf('label %s',u_CLabel_{nB}));
subplot(5,1,5); plot(1:n_E_GENE,mean(E_rank_(tmp_Label_A_ij_,tmp_E_ij_),1) - mean(E_rank_(tmp_Label_B_ij_,tmp_E_ij_),1),'r-',1:n_E_GENE,0*(1:n_E_GENE),'k:');
xlim([1,n_E_GENE]); xlabel('gene'); ylabel('rank diff');
%%%%%%%%;
% Now sort genes by aucX and plot heatmap: ;
%%%%%%%%;
[~,tmp_I_ij_] = sort(I_aucX__(nA,nB,:));
tmp_AB_ = [I_rank_(tmp_Label_A_ij_,tmp_I_ij_) ; I_rank_(tmp_Label_B_ij_,tmp_I_ij_)];
tmp_AB_ = decostand_total_0(tmp_AB_,'row');
clim_ = mean(tmp_AB_(:)) + 1.5*std(tmp_AB_(:))*[-1,+1];
figure(2); clf;
subplot(5,1,[1:2]); imagesc(tmp_AB_(0*tmp_n_A + (1:tmp_n_A),:),clim_); ylabel(sprintf('label %s',u_CLabel_{nA}));
subplot(5,1,[3:4]); imagesc(tmp_AB_(1*tmp_n_A + (1:tmp_n_B),:),clim_); ylabel(sprintf('label %s',u_CLabel_{nB}));
subplot(5,1,5); plot(1:n_I_GENE,mean(I_rank_(tmp_Label_A_ij_,tmp_I_ij_),1) - mean(I_rank_(tmp_Label_B_ij_,tmp_I_ij_),1),'r-',1:n_I_GENE,0*(1:n_I_GENE),'k:');
xlim([1,n_I_GENE]); xlabel('gene'); ylabel('rank diff');
%%%%%%%%;
% Now sort covariates by p-value and plot heatmap: ;
%%%%%%%%;
[~,tmp_C_ij_] = sort(C_aucX__(nA,nB,:));
tmp_AB_ = [C_val_(tmp_Label_A_ij_,tmp_C_ij_) ; C_val_(tmp_Label_B_ij_,tmp_C_ij_)];
tmp_AB_ = decostand_total_0(tmp_AB_,'row');
clim_ = mean(tmp_AB_(:)) + 1.5*std(tmp_AB_(:))*[-1,+1];
figure(3); clf;
subplot(5,1,[1:2]); imagesc(tmp_AB_(0*tmp_n_A + (1:tmp_n_A),:),clim_); ylabel(sprintf('label %s',u_CLabel_{nA}));
subplot(5,1,[3:4]); imagesc(tmp_AB_(1*tmp_n_A + (1:tmp_n_B),:),clim_); ylabel(sprintf('label %s',u_CLabel_{nB}));
subplot(5,1,5); plot(1:n_CCOV,mean(C_rank_(tmp_Label_A_ij_,tmp_C_ij_),1) - mean(C_rank_(tmp_Label_B_ij_,tmp_C_ij_),1),'r-',1:n_CCOV,0*(1:n_CCOV),'k:');
xlim([1,n_CCOV]); xlabel('ccov'); ylabel('rank diff');
clear tmp_I_ij_ tmp_C_ij_ tmp_E_ij_ ;
clear tmp_Label_A_ij_ tmp_n_A ;
clear tmp_Label_B_ij_ tmp_n_B ;
end;%if flag_test;

%%%%%%%%;
% For each gene and covariate, we can record the sorted-list of Label-pairs from I_aucX__ and C_aucX__. ;
% We can compare this to the distribution of sorted-lists obtained in a (Label-Label) permutation-test. ;
% Similarly, for each Label-pair, we can record the sorted-list of genes and covariates from I_aucX__ and C_aucX__. ;
% Again, we can compare this to the distribution of sorted-lists obtained in a (Label-Label) permutation-test. ;
% Note that, for now, we perform this only for the I_ genes, leaving the E_ genes for a machine with more ram. ;
%%%%%%%%;
AB_I_sort_list_ori__ = zeros(n_CLabel,n_CLabel,n_I_GENE);
AB_C_sort_list_ori__ = zeros(n_CLabel,n_CLabel,n_CCOV);
AB_K_sort_list_ori__ = zeros(n_CLabel,n_CLabel,n_CCOV);
I_AB_sort_list_ori__ = zeros(n_CLabel.^2,n_I_GENE);
C_AB_sort_list_ori__ = zeros(n_CLabel.^2,n_CCOV);
K_AB_sort_list_ori__ = zeros(n_CLabel.^2,n_CCOV);
AB_I_sort_list_ij__ = zeros(n_CLabel,n_CLabel,n_I_GENE);
AB_C_sort_list_ij__ = zeros(n_CLabel,n_CLabel,n_CCOV);
AB_K_sort_list_ij__ = zeros(n_CLabel,n_CLabel,n_CCOV);
I_AB_sort_list_ij__ = zeros(n_CLabel.^2,n_I_GENE);
C_AB_sort_list_ij__ = zeros(n_CLabel.^2,n_CCOV);
K_AB_sort_list_ij__ = zeros(n_CLabel.^2,n_CCOV);
n_iteration = 64;
AB_I_sort_list_avg__ = zeros(n_CLabel,n_CLabel,n_I_GENE);
AB_C_sort_list_avg__ = zeros(n_CLabel,n_CLabel,n_CCOV);
AB_K_sort_list_avg__ = zeros(n_CLabel,n_CLabel,n_CCOV);
I_AB_sort_list_avg__ = zeros(n_CLabel.^2,n_I_GENE);
C_AB_sort_list_avg__ = zeros(n_CLabel.^2,n_CCOV);
K_AB_sort_list_avg__ = zeros(n_CLabel.^2,n_CCOV);
AB_I_sort_list_std__ = zeros(n_CLabel,n_CLabel,n_I_GENE);
AB_C_sort_list_std__ = zeros(n_CLabel,n_CLabel,n_CCOV);
AB_K_sort_list_std__ = zeros(n_CLabel,n_CLabel,n_CCOV);
I_AB_sort_list_std__ = zeros(n_CLabel.^2,n_I_GENE);
C_AB_sort_list_std__ = zeros(n_CLabel.^2,n_CCOV);
K_AB_sort_list_std__ = zeros(n_CLabel.^2,n_CCOV);
%%%%%%%%;
AB_I_sort_list_p__ = zeros(n_CLabel,n_CLabel,n_I_GENE);
AB_C_sort_list_p__ = zeros(n_CLabel,n_CLabel,n_CCOV);
AB_K_sort_list_p__ = zeros(n_CLabel,n_CLabel,n_CCOV);
I_AB_sort_list_p__ = zeros(n_CLabel.^2,n_I_GENE);
C_AB_sort_list_p__ = zeros(n_CLabel.^2,n_CCOV);
K_AB_sort_list_p__ = zeros(n_CLabel.^2,n_CCOV);
%%%%%%%%;
C_Label_ = C_.Cluster_ID_20161007(C_to_u_ID_);
u_CLabel_ = unique(C_Label_); n_CLabel = length(u_CLabel_);
rng(1);
for niteration=0:n_iteration;
if (mod(niteration,1)==0); disp(sprintf(' %% niteration %d/%d',niteration,n_iteration)); end;
if (niteration==0); tmp_p_ = 1:n_u; else; tmp_p_ = randperm(n_u); end;
tmp_I_aucX__ = zeros(n_CLabel,n_CLabel,n_I_GENE);
tmp_C_aucX__ = zeros(n_CLabel,n_CLabel,n_CCOV);
tmp_K_aucX__ = zeros(n_CLabel,n_CLabel,n_CCOV);
tmp_I_rank_avg_ = zeros(n_CLabel,n_I_GENE);
tmp_C_rank_avg_ = zeros(n_CLabel,n_CCOV);
tmp_K_rank_avg_ = zeros(n_CLabel,n_CCOV);
%%%%%%%%;
for nA = 1:n_CLabel;
tmp_Label_A_ij_ = find(strcmp(C_Label_(tmp_p_),u_CLabel_{nA})); tmp_n_A = length(tmp_Label_A_ij_);
tmp_I_rank_avg_(nA,:) = mean(I_rank_(tmp_Label_A_ij_,:),1);
tmp_C_rank_avg_(nA,:) = mean(C_rank_(tmp_Label_A_ij_,:),1);
tmp_K_rank_avg_(nA,:) = mean((I_rank_(tmp_Label_A_ij_,:)*beta_I_un_(:,1:n_beta_I))*transpose(beta_I_vn_(:,1:n_beta_I)),1);
end;%for nA = 1:n_CLabel;
%%%%%%%%;
for nA = 1:n_CLabel;
for nB = nA+1:n_CLabel;
tmp_I_aucX__(nA,nB,:) = tmp_I_rank_avg_(nA,:) - tmp_I_rank_avg_(nB,:);
tmp_I_aucX__(nB,nA,:) = -tmp_I_aucX__(nA,nB,:);
tmp_C_aucX__(nA,nB,:) = tmp_C_rank_avg_(nA,:) - tmp_C_rank_avg_(nB,:);
tmp_C_aucX__(nB,nA,:) = -tmp_C_aucX__(nA,nB,:);
tmp_K_aucX__(nA,nB,:) = [tmp_C_rank_avg_(nA,:) - tmp_C_rank_avg_(nB,:)] - [tmp_K_rank_avg_(nA,:) - tmp_K_rank_avg_(nB,:)];
tmp_K_aucX__(nB,nA,:) = -tmp_K_aucX__(nA,nB,:);
end;%for nB = nA+1:n_CLabel;
end;%for nA = 1:n_CLabel;
%%%%%%%%;
[tmp_AB_I_sort_list__,tmp_AB_I_sort_ij__] = sort(tmp_I_aucX__,3,'ascend');
[tmp_AB_C_sort_list__,tmp_AB_C_sort_ij__] = sort(tmp_C_aucX__,3,'ascend');
[tmp_AB_K_sort_list__,tmp_AB_K_sort_ij__] = sort(tmp_K_aucX__,3,'ascend');
[tmp_I_AB_sort_list__,tmp_I_AB_sort_ij__] = sort(reshape(tmp_I_aucX__,n_CLabel.^2,n_I_GENE),2,'ascend');
[tmp_C_AB_sort_list__,tmp_C_AB_sort_ij__] = sort(reshape(tmp_C_aucX__,n_CLabel.^2,n_CCOV),2,'ascend');
[tmp_K_AB_sort_list__,tmp_K_AB_sort_ij__] = sort(reshape(tmp_K_aucX__,n_CLabel.^2,n_CCOV),2,'ascend');
%%%%%%%%;
if (niteration==0); 
AB_I_sort_list_ori__ = tmp_AB_I_sort_list__;
AB_C_sort_list_ori__ = tmp_AB_C_sort_list__;
AB_K_sort_list_ori__ = tmp_AB_K_sort_list__;
I_AB_sort_list_ori__ = tmp_I_AB_sort_list__;
C_AB_sort_list_ori__ = tmp_C_AB_sort_list__;
K_AB_sort_list_ori__ = tmp_K_AB_sort_list__;
AB_I_sort_list_ij__ = tmp_AB_I_sort_ij__;
AB_C_sort_list_ij__ = tmp_AB_C_sort_ij__;
AB_K_sort_list_ij__ = tmp_AB_K_sort_ij__;
I_AB_sort_list_ij__ = tmp_I_AB_sort_ij__;
C_AB_sort_list_ij__ = tmp_C_AB_sort_ij__;
K_AB_sort_list_ij__ = tmp_K_AB_sort_ij__;
end; %if (niteration==0);
if (niteration>0);
AB_I_sort_list_avg__ = AB_I_sort_list_avg__ + tmp_AB_I_sort_list__;
AB_C_sort_list_avg__ = AB_C_sort_list_avg__ + tmp_AB_C_sort_list__;
AB_K_sort_list_avg__ = AB_K_sort_list_avg__ + tmp_AB_K_sort_list__;
I_AB_sort_list_avg__ = I_AB_sort_list_avg__ + tmp_I_AB_sort_list__;
C_AB_sort_list_avg__ = C_AB_sort_list_avg__ + tmp_C_AB_sort_list__;
K_AB_sort_list_avg__ = K_AB_sort_list_avg__ + tmp_K_AB_sort_list__;
AB_I_sort_list_std__ = AB_I_sort_list_std__ + tmp_AB_I_sort_list__.^2;
AB_C_sort_list_std__ = AB_C_sort_list_std__ + tmp_AB_C_sort_list__.^2;
AB_K_sort_list_std__ = AB_K_sort_list_std__ + tmp_AB_K_sort_list__.^2;
I_AB_sort_list_std__ = I_AB_sort_list_std__ + tmp_I_AB_sort_list__.^2;
C_AB_sort_list_std__ = C_AB_sort_list_std__ + tmp_C_AB_sort_list__.^2;
K_AB_sort_list_std__ = K_AB_sort_list_std__ + tmp_K_AB_sort_list__.^2;
end;%if (niteration>0);
%%%%%%%%;
clear tmp_I_aucX__ tmp_C_aucX__ tmp_K_aucX__ tmp_I_rank_avg_ tmp_C_rank_avg_ tmp_K_rank_avg_;
clear tmp_Label_A_ij_ tmp_n_A tmp_I_rank_avg_ tmp_C_rank_avg_ tmp_K_rank_avg_ ;
clear tmp_I_aucX__ tmp_C_aucX__ tmp_K_aucX__ ;
clear tmp_AB_I_sort_list__ tmp_AB_C_sort_list__ tmp_AB_K_sort_list__ tmp_I_AB_sort_list__ tmp_C_AB_sort_list__ tmp_K_AB_sort_list__ ;
clear tmp_AB_I_sort_ij__ tmp_AB_C_sort_ij__ tmp_AB_K_sort_ij__ tmp_I_AB_sort_ij__ tmp_C_AB_sort_ij__ tmp_K_AB_sort_ij__ ;
end;%for niteration=0:n_iteration;
AB_I_sort_list_avg__ = AB_I_sort_list_avg__/n_iteration;
AB_C_sort_list_avg__ = AB_C_sort_list_avg__/n_iteration;
AB_K_sort_list_avg__ = AB_K_sort_list_avg__/n_iteration;
I_AB_sort_list_avg__ = I_AB_sort_list_avg__/n_iteration;
C_AB_sort_list_avg__ = C_AB_sort_list_avg__/n_iteration;
K_AB_sort_list_avg__ = K_AB_sort_list_avg__/n_iteration;
AB_I_sort_list_std__ = sqrt(AB_I_sort_list_std__/n_iteration - AB_I_sort_list_avg__.^2);
AB_C_sort_list_std__ = sqrt(AB_C_sort_list_std__/n_iteration - AB_C_sort_list_avg__.^2);
AB_K_sort_list_std__ = sqrt(AB_K_sort_list_std__/n_iteration - AB_K_sort_list_avg__.^2);
I_AB_sort_list_std__ = sqrt(I_AB_sort_list_std__/n_iteration - I_AB_sort_list_avg__.^2);
C_AB_sort_list_std__ = sqrt(C_AB_sort_list_std__/n_iteration - C_AB_sort_list_avg__.^2);
K_AB_sort_list_std__ = sqrt(K_AB_sort_list_std__/n_iteration - K_AB_sort_list_avg__.^2);
%%%%%%%%;
AB_I_sort_list_p__ = erfc( abs(AB_I_sort_list_ori__ - AB_I_sort_list_avg__)./AB_I_sort_list_std__ / sqrt(2) );
AB_C_sort_list_p__ = erfc( abs(AB_C_sort_list_ori__ - AB_C_sort_list_avg__)./AB_C_sort_list_std__ / sqrt(2) );
AB_K_sort_list_p__ = erfc( abs(AB_K_sort_list_ori__ - AB_K_sort_list_avg__)./AB_K_sort_list_std__ / sqrt(2) );
I_AB_sort_list_p__ = erfc( abs(I_AB_sort_list_ori__ - I_AB_sort_list_avg__)./I_AB_sort_list_std__ / sqrt(2) );
C_AB_sort_list_p__ = erfc( abs(C_AB_sort_list_ori__ - C_AB_sort_list_avg__)./C_AB_sort_list_std__ / sqrt(2) );
K_AB_sort_list_p__ = erfc( abs(K_AB_sort_list_ori__ - K_AB_sort_list_avg__)./K_AB_sort_list_std__ / sqrt(2) );

%%%%%%%%;
% Now run through the cluster pairs, ;
% accumulating the covariates which contribute most to the differences. ;
%%%%%%%%;
AB_I_absZ_avg__ = zeros(n_I_GENE,1);
AB_C_absZ_avg__ = zeros(n_CCOV,1);
AB_K_absZ_avg__ = zeros(n_CCOV,1);
tmp_s = (n_CLabel - 1) * n_CLabel / 2 ;
for nA = 1:n_CLabel;
for nB = nA+1:n_CLabel;
tmp_ij_ = squeeze(AB_I_sort_list_ij__(nA,nB,:)); tmp_z_ = squeeze( abs(AB_I_sort_list_ori__(nA,nB,:) - AB_I_sort_list_avg__(nA,nB,:))./AB_I_sort_list_std__(nA,nB,:) / sqrt(2) ); AB_I_absZ_avg__(tmp_ij_) = AB_I_absZ_avg__(tmp_ij_) + tmp_z_/tmp_s;
tmp_ij_ = squeeze(AB_C_sort_list_ij__(nA,nB,:)); tmp_z_ = squeeze( abs(AB_C_sort_list_ori__(nA,nB,:) - AB_C_sort_list_avg__(nA,nB,:))./AB_C_sort_list_std__(nA,nB,:) / sqrt(2) ); AB_C_absZ_avg__(tmp_ij_) = AB_C_absZ_avg__(tmp_ij_) + tmp_z_/tmp_s;
tmp_ij_ = squeeze(AB_K_sort_list_ij__(nA,nB,:)); tmp_z_ = squeeze( abs(AB_K_sort_list_ori__(nA,nB,:) - AB_K_sort_list_avg__(nA,nB,:))./AB_K_sort_list_std__(nA,nB,:) / sqrt(2) ); AB_K_absZ_avg__(tmp_ij_) = AB_K_absZ_avg__(tmp_ij_) + tmp_z_/tmp_s;
end;%for nB = nA+1:n_CLabel;
end;%for nA = 1:n_CLabel;
clear tmp_ij_ tmp_z_ tmp_s;

flag_disp=1;
if flag_disp;
figure(1);clf;
subplot(2,1,1);
[~,tmp_ij_] = sort(AB_C_absZ_avg__,'ascend');
plot(1:n_CCOV,sort(AB_C_absZ_avg__,'ascend'),'ko');
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(tmp_ij_),'TickLabelInterpreter','none'); xtickangle(90);
xlim([1,n_CCOV]); ylabel('accumulated abs(z-score)'); title('Covariates sorted by raw relevance, accumulated across all label-pairs');
set(gca,'FontSize',7); %set(gca,'FontName','Courier');
subplot(2,1,2);
[~,tmp_ij_] = sort(AB_K_absZ_avg__,'ascend');
plot(1:n_CCOV,sort(AB_K_absZ_avg__,'ascend'),'ko');
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(tmp_ij_),'TickLabelInterpreter','none'); xtickangle(90);
xlim([1,n_CCOV]); ylabel('accumulated abs(z-score)'); title('Covariates sorted by relative relevance, accumulated across all label-pairs');
set(gca,'FontSize',7); %set(gca,'FontName','Courier');
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_base = sprintf('%s/dir_jpg/AB_X_absZ_avg',dir_trunk);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_disp;

flag_test=1;
if flag_test;
%%%%%%%;
% Now we can test a pair of labels: ;
%%%%%%%%;
for nl=1:3;
if nl==1; nA=find(strcmp(u_CLabel_,'28')); nB=find(strcmp(u_CLabel_,'33')); end;
if nl==2; nA=find(strcmp(u_CLabel_,'22')); nB=find(strcmp(u_CLabel_,'11')); end;
if nl==3; nA=find(strcmp(u_CLabel_,'19')); nB=find(strcmp(u_CLabel_,'15')); end;
disp(sprintf(' %% Comparing Label %d (%s) and %d (%s):',nA,u_CLabel_{nA},nB,u_CLabel_{nB}));
%%%%%%%%;
figure(1);clf; hold on;
plot(squeeze(AB_C_sort_list_avg__(nA,nB,:)+2.0*AB_C_sort_list_std__(nA,nB,:)),'k','LineWidth',0.5);
plot(squeeze(AB_C_sort_list_avg__(nA,nB,:)+0.0*AB_C_sort_list_std__(nA,nB,:)),'k','LineWidth',1.0);
plot(squeeze(AB_C_sort_list_avg__(nA,nB,:)-2.0*AB_C_sort_list_std__(nA,nB,:)),'k','LineWidth',0.5);
plot(squeeze(AB_C_sort_list_ori__(nA,nB,:)),'r','LineWidth',2.0);
xlim([1,n_CCOV]); xlabel('covariate'); ylabel('rank difference (auc)'); title(sprintf('sorted auc of covariates influencing labels %s-vs-%s',u_CLabel_{nA},u_CLabel_{nB}));
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(AB_C_sort_list_ij__(nA,nB,:)),'TickLabelInterpreter','none'); xtickangle(90);
set(gca,'FontSize',7);
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_base = sprintf('%s/dir_jpg/AB_C_sort_list_ori_%s_vs_%s',dir_trunk,u_CLabel_{nA},u_CLabel_{nB});
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
hold off;
%%%%%%%%;
figure(2);clf; hold on;
plot(-log10(squeeze(AB_C_sort_list_p__(nA,nB,:))),'r','LineWidth',2.0);
plot(1:n_CCOV,-log10(0.05)*ones(n_CCOV,1),'k--','LineWidth',2);
plot(1:n_CCOV,-log10(0.001)*ones(n_CCOV,1),'k:','LineWidth',2);
xlim([1,n_CCOV]); xlabel('covariate'); ylim([0,15]); ylabel('-log10(p)'); title(sprintf('-log10(p-value) for sorted auc of covariates influencing labels %s-vs-%s',u_CLabel_{nA},u_CLabel_{nB}));
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(AB_C_sort_list_ij__(nA,nB,:)),'TickLabelInterpreter','none'); xtickangle(90);
set(gca,'FontSize',7);
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_base = sprintf('%s/dir_jpg/AB_C_sort_list_p_%s_vs_%s',dir_trunk,u_CLabel_{nA},u_CLabel_{nB});
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
hold off;
%%%%%%%%;
figure(3);clf; hold on;
plot(squeeze(AB_K_sort_list_avg__(nA,nB,:)+2.0*AB_K_sort_list_std__(nA,nB,:)),'k','LineWidth',0.5);
plot(squeeze(AB_K_sort_list_avg__(nA,nB,:)+0.0*AB_K_sort_list_std__(nA,nB,:)),'k','LineWidth',1.0);
plot(squeeze(AB_K_sort_list_avg__(nA,nB,:)-2.0*AB_K_sort_list_std__(nA,nB,:)),'k','LineWidth',0.5);
plot(squeeze(AB_K_sort_list_ori__(nA,nB,:)),'r','LineWidth',2.0);
xlim([1,n_CCOV]); xlabel('covariate'); ylabel('rank difference (auc)'); title(sprintf('sorted auc of covariates influencing labels %s-vs-%s',u_CLabel_{nA},u_CLabel_{nB}));
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(AB_K_sort_list_ij__(nA,nB,:)),'TickLabelInterpreter','none'); xtickangle(90);
set(gca,'FontSize',7);
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_base = sprintf('%s/dir_jpg/AB_K_sort_list_ori_%s_vs_%s',dir_trunk,u_CLabel_{nA},u_CLabel_{nB});
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
hold off;
%%%%%%%%;
figure(4);clf; hold on;
plot(-log10(squeeze(AB_K_sort_list_p__(nA,nB,:))),'r','LineWidth',2.0);
plot(1:n_CCOV,-log10(0.05)*ones(n_CCOV,1),'k--','LineWidth',2);
plot(1:n_CCOV,-log10(0.001)*ones(n_CCOV,1),'k:','LineWidth',2);
xlim([1,n_CCOV]); xlabel('covariate'); ylim([0,15]); ylabel('-log10(p)'); title(sprintf('-log10(p-value) for sorted auc of covariates influencing labels %s-vs-%s',u_CLabel_{nA},u_CLabel_{nB}));
set(gca,'XTick',1:n_CCOV,'XTickLabel',C_VariableName_LC_(AB_K_sort_list_ij__(nA,nB,:)),'TickLabelInterpreter','none'); xtickangle(90);
set(gca,'FontSize',7);
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname_base = sprintf('%s/dir_jpg/AB_K_sort_list_p_%s_vs_%s',dir_trunk,u_CLabel_{nA},u_CLabel_{nB});
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
hold off;
%%%%%%%%;
end;%for nl=1:3;
end;%if flag_test;

save(sprintf('%s/test_loader_7.mat',dir_data),'-v7.3'...
,'AB_C_absZ_avg__'...
,'AB_C_sort_list_avg__'...
,'AB_C_sort_list_ij__'...
,'AB_C_sort_list_ori__'...
,'AB_C_sort_list_p__'...
,'AB_C_sort_list_std__'...
,'AB_I_absZ_avg__'...
,'AB_I_sort_list_avg__'...
,'AB_I_sort_list_ij__'...
,'AB_I_sort_list_ori__'...
,'AB_I_sort_list_p__'...
,'AB_I_sort_list_std__'...
,'AB_K_absZ_avg__'...
,'AB_K_sort_list_avg__'...
,'AB_K_sort_list_ij__'...
,'AB_K_sort_list_ori__'...
,'AB_K_sort_list_p__'...
,'AB_K_sort_list_std__'...
,'B_LC_'...
,'B_VariableName_'...
,'B_col_val_'...
,'B_rank_'...
,'B_val_'...
,'C_'...
,'C_AB_sort_list_avg__'...
,'C_AB_sort_list_ij__'...
,'C_AB_sort_list_ori__'...
,'C_AB_sort_list_p__'...
,'C_AB_sort_list_std__'...
,'C_ID_'...
,'C_LC_'...
,'C_Label_'...
,'C_VariableName_'...
,'C_VariableName_LC_'...
,'C_aucX__'...
,'C_rank_'...
,'C_rank_avg_'...
,'C_to_u_ID_'...
,'C_val_'...
,'CtEn_rank_U_ij_'...
,'CtEn_rank_V_ij_'...
,'CtEn_rank_Z_'...
,'CtEn_rank_avg_'...
,'CtEn_rank_ori_'...
,'CtEn_rank_p_'...
,'CtEn_rank_std_'...
,'CtIn_rank_U_ij_'...
,'CtIn_rank_V_ij_'...
,'CtIn_rank_Z_'...
,'CtIn_rank_avg_'...
,'CtIn_rank_ori_'...
,'CtIn_rank_p_'...
,'CtIn_rank_std_'...
,'D_LC_'...
,'D_VariableName_'...
,'D_col_val'...
,'D_col_val_'...
,'D_rank_'...
,'D_val_'...
,'G_'...
,'G_h_'...
,'G_h_mask_'...
,'I_ID_'...
,'I_VariableName_'...
,'I_VariableName_to_cap_'...
,'I_col_val_'...
,'I_length_'...
,'I_rank_avg_'...
,'I_to_u_ID_'...
,'K_AB_sort_list_avg__'...
,'K_AB_sort_list_ij__'...
,'K_AB_sort_list_ori__'...
,'K_AB_sort_list_p__'...
,'K_AB_sort_list_std__'...
,'K_aucX__'...
,'K_rank_avg_'...
,'Label_ID_'...
,'ans'...
,'beta_I_un_'...
,'beta_I_vn_'...
,'clim_'...
,'dir_data'...
,'dir_trunk'...
,'flag_disp'...
,'flag_load'...
,'flag_plot'...
,'flag_test'...
,'fname_base'...
,'gene_length_'...
,'gene_length_to_cap_'...
,'nA'...
,'nB'...
,'nCCOV'...
,'nCLabel'...
,'nC_ID'...
,'nE_ID'...
,'nG'...
,'nI_ID'...
,'nLabel_ID'...
,'n_BCOV'...
,'n_CCOV'...
,'n_CLabel'...
,'n_CLabel_'...
,'n_DCOV'...
,'n_D_col_val'...
,'n_E_GENE'...
,'n_FACTOR'...
,'n_G'...
,'n_I_GENE'...
,'n_Label_ID'...
,'n_Label_ID_'...
,'n_beta_E'...
,'n_beta_I'...
,'n_bin'...
,'n_i'...
,'n_iteration'...
,'n_rank'...
,'n_u'...
,'nc'...
,'ne'...
,'ng'...
,'ni'...
,'niteration'...
,'nl'...
,'nrank'...
,'nu'...
,'nx'...
,'str_exclude_'...
,'str_table_name'...
,'u_CLabel_'...
,'u_ID_'...
);


disp('returning');return;



