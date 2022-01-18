%%%%%%%%;
% data-set at: ;
% https://support.10xgenomics.com/single-cell-gene-expression/datasets/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_PBMCs_human_1 ;
%%%%%%%%;

clear;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

flag_recalc = 0;
flag_replot = 0;
tolerance_master = 1e-2;
nf=0;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root);
str_code = sprintf('%s/halfloop_dev',dir_code);
dir_data = sprintf('/%s/rangan/dir_bcc/dir_10xgenomics/dir_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_PBMCs_human_1',string_root);
dir_mat = sprintf('%s/dir_mat',dir_data); if ~exist(dir_mat); disp(sprintf(' %% mkdir %s',dir_mat)); mkdir(sprintf('%s',dir_mat)); end;
dir_jpg = sprintf('%s/dir_jpg',dir_data); if ~exist(dir_jpg); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(sprintf('%s',dir_jpg)); end;

fname_tsne = sprintf('%s/analysis/tsne/2_components/projection.csv',dir_data);
fp = fopen(fname_tsne,'r'); tmp_ = textscan(fp,'%s%s%s','headerlines',1,'Delimiter',','); fclose(fp);
loading_tsne_ = [ cellfun(@str2num,tmp_{2}) , cellfun(@str2num,tmp_{3}) ];
cell_id_tsne_ = tmp_{1};
clear tmp_ ;
fname_tsne = sprintf('%s/analysis/tsne/multiplexing_capture_2_components/projection.csv',dir_data);
fp = fopen(fname_tsne,'r'); tmp_ = textscan(fp,'%s%s%s','headerlines',1,'Delimiter',','); fclose(fp);
loading_mtsne_ = [ cellfun(@str2num,tmp_{2}) , cellfun(@str2num,tmp_{3}) ];
cell_id_mtsne_ = tmp_{1};
clear tmp_ ;

fname_azimuth_pred = sprintf('%s/azimuth_pred.tsv',dir_data);
fp = fopen(fname_azimuth_pred,'r'); tmp_ = textscan(fp,'%s%s%s%s','headerlines',1,'Delimiter','\t'); fclose(fp);
azimuth_pred_cell_id_ = tmp_{1};
azimuth_pred_celltype_ = tmp_{2};
azimuth_pred_celltype_score_ = cellfun(@str2num,tmp_{3});
azimuth_pred_mapping_score_ = cellfun(@str2num,tmp_{4});

%%%%%%%%;
% assert that the cell_ids are the same. ;
%%%%%%%%;
for nl=0:numel(cell_id_mtsne_)-1;
assert(strcmp(cell_id_mtsne_(1+nl),azimuth_pred_cell_id_(1+nl)));
assert(strcmp(cell_id_mtsne_(1+nl),cell_id_tsne_(1+nl)));
assert(strcmp(cell_id_tsne_(1+nl),azimuth_pred_cell_id_(1+nl)));
end;%for nl=0:numel(cell_id_mtsne_)-1;

%%%%%%%%;
% also load and check: ;
% azimuth_cell_label_.csv ;
% azimuth_umap__.mat ;
%%%%%%%%;
fname_azimuth_cell_label = sprintf('%s/azimuth_cell_label_.csv',dir_data);
fp = fopen(fname_azimuth_cell_label,'r'); tmp_ = textscan(fp,'%s','headerlines',1,'Delimiter','\n'); fclose(fp);
azimuth_cell_label_cell_id_ = tmp_{1};
for nl=0:numel(cell_id_mtsne_)-1;
assert(strfind(azimuth_cell_label_cell_id_{1+nl},azimuth_pred_cell_id_{1+nl})==2);
end;%for nl=0:numel(cell_id_mtsne_)-1;
fname_mat = sprintf('%s/azimuth_umap__.mat',dir_data);
azimuth_umap__ = load(fname_mat); loading_umap_ = azimuth_umap__.A;

[ ...
 celltype_enum_ ...
,n_u_celltype ...
,u_celltype_ ...
,index_nu_from_nall_ ...
,n_u_celltype_ ...
,index_nall_from_nu__ ...
] = ...
label_str_to_enum_1( ...
 azimuth_pred_celltype_ ...
);
[~,index_sort_u_celltype_] = sort(n_u_celltype_,'descend'); index_sort_u_celltype_ = index_sort_u_celltype_ - 1;
index_nu_srt_from_nu_ori_ = index_sort_u_celltype_;
n_all = sum(n_u_celltype_);
for nu_celltype=0:n_u_celltype-1;
tmp_n_u = n_u_celltype_(1+index_sort_u_celltype_(1+nu_celltype));
u_celltype = u_celltype_{1+index_sort_u_celltype_(1+nu_celltype)};
disp(sprintf(' %% nu %.2d/%.2d: %.4d/%.4d = %0.2f = %s',nu_celltype,n_u_celltype,tmp_n_u,n_all,tmp_n_u/n_all,u_celltype));
end;%for nu_celltype=0:n_u_celltype-1;

%%%%%%%%;
fname_fig = sprintf('%s/umap_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
markersize_final = 12;
markersize_start = 4;
str_symbol_ = {'o','^','s','p','h'};
c_lines__ = colormap('lines'); n_c_lines = size(c_lines__,1);
hold on;
for nl=0:numel(index_sort_u_celltype_)-1;
nc_lines = nl;
str_symbol = str_symbol_{1+mod(nl,5)};
markersize_use = round(markersize_start + (markersize_final-markersize_start)*nl/(numel(index_sort_u_celltype_)-1));
tmp_nu_celltype = index_sort_u_celltype_(1+nl);
tmp_index_ = index_nall_from_nu__{1+tmp_nu_celltype};
plot(loading_umap_(1+tmp_index_,1),loading_umap_(1+tmp_index_,2),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_lines__(1+nc_lines,:),'MarkerEdgeColor','k');
end;%for nl=0:numel(index_sort_u_celltype_)-1;
axis equal; axisnotick;
legend(u_celltype_(1+index_sort_u_celltype_));
title('umap');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/tsne_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
markersize_final = 12;
markersize_start = 4;
str_symbol_ = {'o','^','s','p','h'};
c_lines__ = colormap('lines'); n_c_lines = size(c_lines__,1);
hold on;
for nl=0:numel(index_sort_u_celltype_)-1;
nc_lines = nl;
str_symbol = str_symbol_{1+mod(nl,5)};
markersize_use = round(markersize_start + (markersize_final-markersize_start)*nl/(numel(index_sort_u_celltype_)-1));
tmp_nu_celltype = index_sort_u_celltype_(1+nl);
tmp_index_ = index_nall_from_nu__{1+tmp_nu_celltype};
plot(loading_tsne_(1+tmp_index_,1),loading_tsne_(1+tmp_index_,2),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_lines__(1+nc_lines,:),'MarkerEdgeColor','k');
end;%for nl=0:numel(index_sort_u_celltype_)-1;
axis equal; axisnotick;
legend(u_celltype_(1+index_sort_u_celltype_));
title('tsne');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/mtsne_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
markersize_final = 12;
markersize_start = 4;
str_symbol_ = {'o','^','s','p','h'};
c_lines__ = colormap('lines'); n_c_lines = size(c_lines__,1);
hold on;
for nl=0:numel(index_sort_u_celltype_)-1;
nc_lines = nl;
str_symbol = str_symbol_{1+mod(nl,5)};
markersize_use = round(markersize_start + (markersize_final-markersize_start)*nl/(numel(index_sort_u_celltype_)-1));
tmp_nu_celltype = index_sort_u_celltype_(1+nl);
tmp_index_ = index_nall_from_nu__{1+tmp_nu_celltype};
plot(loading_mtsne_(1+tmp_index_,1),loading_mtsne_(1+tmp_index_,2),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_lines__(1+nc_lines,:),'MarkerEdgeColor','k');
end;%for nl=0:numel(index_sort_u_celltype_)-1;
axis equal; axisnotick;
legend(u_celltype_(1+index_sort_u_celltype_));
title('mtsne');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

fname_mtx = sprintf('%s/sample_feature_bc_matrix/matrix.mtx',dir_data);
A_gs__ = read_mtx_0(fname_mtx);
[n_gene,n_smp] = size(A_gs__);
nnz_each_ = A_gs__(:);
nnz_gene_ = sum(A_gs__,2);
nnz_cell_ = sum(A_gs__,1);

%%%%%%%%;
% Now we could fit nnz_gene_ with an algebraic. ;
% assumed form is given by: ;
% p(n;alpha,j) = 1/(n+alpha)^(1+j). ;
%%%%%%%%;
% n_N = 64; plot(log(1+[0:n_N]),log(hist(nnz_gene_,[0:n_N])),'ro',log([0:n_N]),log(1./[0:n_N].^1.2) + 8.5,'k-');

%%%%%%%%;
fname_fig = sprintf('%s/nnz_each',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
n_N = 64;
subplot(1,2,1); semilogy(hist(nnz_each_,[0:n_N]),'ko-'); xlim([-1,1+n_N]); xlabel('n'); ylabel('hist');
subplot(1,2,2); loglog(hist(nnz_each_,[0:n_N]),'ko-');  xlabel('n'); ylabel('hist');
sgtitle('nnz_each','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/nnz_gene',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
n_N = 64;
subplot(1,2,1); semilogy(hist(nnz_gene_,[0:n_N]),'ko-'); xlim([-1,1+n_N]); xlabel('n'); ylabel('hist');
subplot(1,2,2); loglog(hist(nnz_gene_,[0:n_N]),'ko-');  xlabel('n'); ylabel('hist');
sgtitle('nnz_gene','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/nnz_cell',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
h_nnz_ = mean(nnz_cell_) + 4.5*std(nnz_cell_,1)*linspace(-1,1,128);
hist(nnz_cell_,h_nnz_);
xlabel('n'); ylabel('hist');
sgtitle('nnz_cell','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

nnz_gene_ = sum(A_gs__~=0,2);
index_gene_retain_ = efind(nnz_gene_> 1); %<-- at least 2 nonzero entries, since 1 seems to deviate from algebraic decay. ;
A_clr_sg__ = log(1+transpose(A_gs__(1+index_gene_retain_,:)));

[tmp_U__,tmp_S__,tmp_V__] = svds(A_clr_sg__,2); tmp_S_ = diag(tmp_S__);

fname_fig = sprintf('%s/pca_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
markersize_final = 12;
markersize_start = 4;
str_symbol_ = {'o','^','s','p','h'};
c_lines__ = colormap('lines'); n_c_lines = size(c_lines__,1);
hold on;
for nl=0:numel(index_sort_u_celltype_)-1;
nc_lines = nl;
str_symbol = str_symbol_{1+mod(nl,5)};
markersize_use = round(markersize_start + (markersize_final-markersize_start)*nl/(numel(index_sort_u_celltype_)-1));
tmp_nu_celltype = index_sort_u_celltype_(1+nl);
tmp_index_ = index_nall_from_nu__{1+tmp_nu_celltype};
plot(tmp_U__(1+tmp_index_,1),tmp_U__(1+tmp_index_,2),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_lines__(1+nc_lines,:),'MarkerEdgeColor','k');
end;%for nl=0:numel(index_sort_u_celltype_)-1;
axis equal; axisnotick;
legend(u_celltype_(1+index_sort_u_celltype_));
title('pca');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% estimate rank as 8 ;
%%%%%%%%;
% [tmp_U_A_clr__,tmp_S_A_clr__,tmp_V_A_clr__] = svds(mean_center_0(A_clr_sg__),24); tmp_S_A_clr_ = diag(tmp_S_A_clr__);
rank_estimate_A_clr = 8;
[tmp_U__,tmp_S__,tmp_V__] = svds(A_clr_sg__,2*rank_estimate_A_clr); tmp_S_ = diag(tmp_S__);
tmp_US__ = tmp_U__*tmp_S__;

%%%%%%%%;
% cluster full data-set. ;
%%%%%%%%;
flag_calculate = strcmp(platform,'access1');
if flag_calculate;
test_dexcluster_driver_5( ...
 str_code ...
,A_clr_sg__ ...
,index_nu_from_nall_ ...
,dir_mat ...
,rank_estimate_A_clr ...
);
end;%if flag_calculate;

flag_pair=strcmp(platform,'access1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_pair;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now isolate only two cell-types. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
snr_est_from_nu_ori__ = zeros(n_u_celltype,n_u_celltype);
snr_est_from_nu_srt__ = zeros(n_u_celltype,n_u_celltype);
for nl0=0:18; for nl1=nl0+1:18;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
nu_celltype_0 = index_sort_u_celltype_(1+nl0); %<-- most populous. ;
disp(sprintf(' %% nu_celltype_0 %d %s',nu_celltype_0,u_celltype_{1+nu_celltype_0}));
nu_celltype_1 = index_sort_u_celltype_(1+nl1); %<-- most populous. ;
disp(sprintf(' %% nu_celltype_1 %d %s',nu_celltype_1,u_celltype_{1+nu_celltype_1}));
index_nall_from_nu_0_ = index_nall_from_nu__{1+nu_celltype_0};
index_nall_from_nu_1_ = index_nall_from_nu__{1+nu_celltype_1};
str_vs = sprintf('%s_vs_%s',u_celltype_{1+nu_celltype_0},u_celltype_{1+nu_celltype_1}); str_vs(strfind(str_vs,' '))='_';
n_celltype_0 = n_u_celltype_(1+nu_celltype_0);
n_celltype_1 = n_u_celltype_(1+nu_celltype_1);
index_sub_nu_from_nall_ = [ zeros(n_celltype_0,1) ; ones(n_celltype_1,1) ];
%%%%%%%%;
A_sub_gs__ = A_gs__(:,1 + [index_nall_from_nu_0_;index_nall_from_nu_1_]);
nnz_sub_gene_ = sum(A_sub_gs__~=0,2);
index_sub_gene_retain_ = efind(nnz_sub_gene_> 1); %<-- at least 2 nonzero entries. ;
A_sub_clr_sg__ = log(1+transpose(A_sub_gs__(1+index_sub_gene_retain_,:)));
tmp_US_sub__ = tmp_US__(1 + [index_nall_from_nu_0_;index_nall_from_nu_1_],:);
[DDL_all_,USL_n_,VL_n_,SL_n_] = SVD_discat_0(tmp_US_sub__,index_sub_nu_from_nall_,2);
%%%%%%%%;
% Following random_matrix_planted_cluster_0, ;
% we define the signal-to-noise-ratio (snr) to be the ratio between: ;
% numerator: the dominant singular-value of the signal, and ;
% denominator: the dominant singular-value of the noise. ;
% With 2 cell-types we can set up an MA-by-NA sample-by-gene matrix: ;
% A_n = [ B_n - 1_n*C_avg_n ] = [ Bc_n ] ;
%       [ C_n - 1_n*C_avg_n ] = [ Cc_n ] ;
% where C_avg_n = mean(C_n,1), and 1_n is a ones-vector. ;
% Now we can examine the covariance matrix: ;
% A_n*A_t = [ Bc_n*Bc_t , Bc_n*Cc_t ] ;
%           [ Cc_n*Bc_t , Cc_n*Cc_t ] ;
% Now we can estimate the 'noise' as the dominant eigenvalue of Cc_n*Cc_t, ;
% and the 'signal' as the dominant eigenvalue of Bc_n*Bc_t. ;
% The square-root of the ratio of these two will be our estimate for the snr. ;
% Note that this is derived using the covariance-matrix A_n*A_t, ;
% rather than A_t*A_n, which is what we expect when calculating ;
% correlations between rows (i.e., searching for row-clusters). ;
%%%%%%%%;
tmp_Cc_n_ = mean(A_sub_clr_sg__(1:n_celltype_0,:),1);
tmp_A_n_ = A_sub_clr_sg__ - ones(n_celltype_0+n_celltype_1,1)*tmp_Cc_n_;
tmp_Cc_n_ = tmp_A_n_(1:n_celltype_0,:);
tmp_Bc_n_ = tmp_A_n_(n_celltype_0 + [1:n_celltype_1],:);
snr_est = svds(tmp_Bc_n_,1)/svds(tmp_Cc_n_,1);
disp(sprintf(' %% nu_celltype (%d,%d): %s: estimated snr: %0.4f',nu_celltype_0,nu_celltype_1,str_vs,snr_est));
snr_est_from_nu_ori__(1+nu_celltype_0,1+nu_celltype_1) = snr_est;
snr_est_from_nu_srt__(1+nl0,1+nl1) = snr_est;
%%%%%%%%;
% Now plot scatterplot restricted to only these cell-types. ;
%%%%%%%%;
fname_fig = sprintf('%s/pca_%s_FIGA',dir_jpg,str_vs);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
markersize_final = 12;
markersize_start = 4;
str_symbol_ = {'o','^','s','p','h'};
c_lines__ = colormap('lines'); n_c_lines = size(c_lines__,1);
hold on;
for nl=0:numel(index_sort_u_celltype_)-1;
nc_lines = nl;
str_symbol = str_symbol_{1+mod(nl,5)};
markersize_use = round(markersize_start + (markersize_final-markersize_start)*nl/(numel(index_sort_u_celltype_)-1));
tmp_nu_celltype = index_sort_u_celltype_(1+nl);
if (tmp_nu_celltype==nu_celltype_0);
plot(USL_n_(1:n_celltype_0,1),USL_n_(1:n_celltype_0,2),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_lines__(1+nc_lines,:),'MarkerEdgeColor','k');
end;%if (tmp_nu_celltype==nu_celltype_1);
if (tmp_nu_celltype==nu_celltype_1);
plot(USL_n_(n_celltype_0 + [1:n_celltype_1],1),USL_n_(n_celltype_0 + [1:n_celltype_1],2),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_lines__(1+nc_lines,:),'MarkerEdgeColor','k');
end;%if (tmp_nu_celltype==nu_celltype_1);
end;%for nl=0:numel(index_sort_u_celltype_)-1;
axis equal; axisnotick;
legend(u_celltype_(1+[nu_celltype_0 , nu_celltype_1]));
sgtitle(sprintf('%s: snr %0.4f',str_vs,snr_est),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Now call clustering. ;
%%%%%%%%;
flag_calculate = strcmp(platform,'access1');
if flag_calculate;
dir_mat_sub = sprintf('%s/dir_mat_%s',dir_data,str_vs);
if (~exist(dir_mat_sub,'dir')); disp(sprintf(' %% mkdir %s',dir_mat_sub)); mkdir(dir_mat_sub); end;
test_dexcluster_driver_5( ...
 str_code ...
,A_sub_clr_sg__ ...
,index_sub_nu_from_nall_ ...
,dir_mat_sub ...
,rank_estimate_A_clr ...
);
end;%if flag_calculate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;end;%for nl0=0:8; for nl1=nl0+1:8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_pair;

fname_mat = sprintf('%s/snr_est__.mat',dir_mat);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_mat ...
     ,'snr_est_from_nu_ori__' ...
     ,'snr_est_from_nu_srt__' ...
     ,'celltype_enum_' ...
     ,'u_celltype_' ...
     ,'n_u_celltype' ...
     ,'n_u_celltype_' ...
     ,'index_sort_u_celltype_' ...
     );
end;%if (~exist(fname_mat,file'));
load(fname_mat);

disp('returning'); return;

%%%%%%%%;
% collect statistics from full data-set. ;
%%%%%%%%;
[ ...
 flag_found_ ...
,lpv_quad_ ...
,n_cluster_found_ ...
,label_B__ ...
,prefix_method_ ...
,legend_method_ ...
,nrank_method_ ...
,symbol_method_ ...
] = ...
test_dexcluster_driver_collect_5( ...
 str_code ...
,index_nu_from_nall_ ...
,dir_mat ...
,rank_estimate_A_clr ...
);
%%%%%%%%;

%%%%%%%%;
% loading hnbr results from full data-set. ;
%%%%%%%%;
%infix_method_0 = 'hnbrtZRgumb_r8_r1'; infix_method_1 = sprintf('%s',infix_method_0);
infix_method_0 = 'hnbtZRgumb_r1'; infix_method_1 = sprintf('%s_g010_n32',infix_method_0);
dir_mat_method = sprintf('%s/dir_tmp_%s/dir_tmp_%s',dir_mat,infix_method_0,infix_method_1);
fname_output_label = sprintf('%s/output_label__.txt',dir_mat_method);
fp = fopen(fname_output_label); output_label__ = textscan(fp,'%s','Delimiter','\n'); output_label__ = output_label__{1}; fclose(fp);
fname_nlpbra_label = sprintf('%s/nlpbra_label__.txt',dir_mat_method);
fp = fopen(fname_nlpbra_label); nlpbra_label__ = textscan(fp,'%s','Delimiter','\n'); nlpbra_label__ = nlpbra_label__{1}; fclose(fp);
fname_nlpnex_label = sprintf('%s/nlpnex_label__.txt',dir_mat_method);
fp = fopen(fname_nlpnex_label); nlpnex_label__ = textscan(fp,'%s','Delimiter','\n'); nlpnex_label__ = nlpnex_label__{1}; fclose(fp);
[lpv,lP0,fla,cap_,cup_] = label_to_label_enrichment_quad_4(index_nu_from_nall_,label_str_to_enum_0(output_label__));
%%%%%%%%;
label_tree ...
= ...
label_to_tree_1( ...
 0 ...
,[] ...
,[] ...
,[] ...
,output_label__ ...
,nlpbra_label__ ...
,nlpnex_label__ ...
);
%%%%%%%%%%%%%%%%;
fname_fig = sprintf('%s/label_plot_recursive_nlpscale_%s_FIGA',dir_jpg,infix_method_1);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(1);figbig;;
c_nlpvt__ = colormap_nlpvt(); n_c_nlpvt = size(c_nlpvt__,1); nlp_lim_ = [0,27];
subplot(1,9,1);
cla;
tmp_x__ = zeros(4,n_nlpnex_threshold);
tmp_y__ = zeros(4,n_nlpnex_threshold);
tmp_c___ = zeros(1,n_nlpnex_threshold,3);
for nnlpnex_threshold=0:n_nlpnex_threshold-1;
nlpnex_threshold = nlpnex_threshold_(1+nnlpnex_threshold);
nlpv = -lpv_method_(1+nnlpnex_threshold);
tmp_x__(:,1+nnlpnex_threshold) = [0;nlpv;nlpv;0];
tmp_y__(:,1+nnlpnex_threshold) = nlpnex_threshold + 0.5*[-1;-1;+1;+1];
nc_nlpvt = max(0,min(n_c_nlpvt-1,floor(n_c_nlpvt*(nlpnex_threshold-min(nlp_lim_))/diff(nlp_lim_))));
tmp_c___(1,1+nnlpnex_threshold,:) = c_nlpvt__(1+nc_nlpvt,:);
end;%for nnlpnex_threshold=0:n_nlpnex_threshold-1;
hold on;
patch(tmp_x__,tmp_y__,tmp_c___);
hold off;
ylim([min(nlp_lim_),max(nlp_lim_)+1]);
xlim([0,max(-lpv_method_)*1.25]);
xlabel('entropy');
ylabel('nlp');
grid on;
%%%%%%%%;
subplot(1,9,[2:9]);
label_plot_recursive_nlpscale_2(output_label__,nlpbra_label__,nlpnex_label__,[],[],[]);
sgtitle(sprintf('%s',dir_mat_method),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%%%%%%%%%;
nlpnex_threshold_ = 12:1:27; n_nlpnex_threshold = numel(nlpnex_threshold_);
ouput_label_pid___ = cell(n_nlpnex_threshold,1);
lpv_method_ = zeros(n_nlpnex_threshold,1);
n_cluster_method_ = zeros(n_nlpnex_threshold,1);
for nnlpnex_threshold=0:n_nlpnex_threshold-1;
nlpnex_threshold = nlpnex_threshold_(1+nnlpnex_threshold);
tmp_output_label__ = label_tree_output_label_from_nlpnex_0(label_tree,nlpnex_threshold,[]);
lpv_method_(1+nnlpnex_threshold) = label_to_label_enrichment_quad_4(index_nu_from_nall_,label_str_to_enum_0(tmp_output_label__));
n_cluster_method_(1+nnlpnex_threshold) = numel(unique(tmp_output_label__));
output_label_pid___{1+nnlpnex_threshold} = tmp_output_label__;
end;%for nnlpnex_threshold=0:n_nlpnex_threshold-1;
%%%%%%%%%%%%%%%%;
for nnlpnex_threshold=0:n_nlpnex_threshold-1;
nlpnex_threshold = nlpnex_threshold_(1+nnlpnex_threshold);
%%%%;
fname_fig = sprintf('%s/umap_%s_nlp%.5d_FIGA',dir_jpg,infix_method_1,round(1000*nlpnex_threshold));
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
tmp_celltype_ = output_label_pid___{1+nnlpnex_threshold};
[ ...
 tmp_celltype_enum_ ...
,tmp_n_u_celltype ...
,tmp_u_celltype_ ...
,tmp_index_nu_from_nall_ ...
,tmp_n_u_celltype_ ...
,tmp_index_nall_from_nu__ ...
] = ...
label_str_to_enum_1( ...
 tmp_celltype_ ...
);
[~,tmp_index_sort_u_celltype_] = sort(tmp_n_u_celltype_,'descend'); tmp_index_sort_u_celltype_ = tmp_index_sort_u_celltype_ - 1;
tmp_index_nu_srt_from_nu_ori_ = tmp_index_sort_u_celltype_;
tmp_n_all = sum(tmp_n_u_celltype_);
for tmp_nu_celltype=0:tmp_n_u_celltype-1;
tmp_n_u = tmp_n_u_celltype_(1+tmp_index_sort_u_celltype_(1+tmp_nu_celltype));
tmp_u_celltype = tmp_u_celltype_{1+tmp_index_sort_u_celltype_(1+tmp_nu_celltype)};
disp(sprintf(' %% nu %.2d/%.2d: %.4d/%.4d = %0.2f = %s',tmp_nu_celltype,tmp_n_u_celltype,tmp_n_u,tmp_n_all,tmp_n_u/tmp_n_all,tmp_u_celltype));
end;%for tmp_nu_celltype=0:tmp_n_u_celltype-1;
%%%%;
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
markersize_final = 8;
markersize_start = 8;
str_symbol_ = {'o','^','s','p','h'};
c_lines__ = colormap('lines'); n_c_lines = size(c_lines__,1);
hold on;
for nl=0:numel(tmp_index_sort_u_celltype_)-1;
nc_lines = nl;
str_symbol = str_symbol_{1+mod(nl,5)};
markersize_use = round(markersize_start + (markersize_final-markersize_start)*nl/(numel(tmp_index_sort_u_celltype_)-1));
tmp_nu_celltype = tmp_index_sort_u_celltype_(1+nl);
tmp_index_ = tmp_index_nall_from_nu__{1+tmp_nu_celltype};
plot(loading_umap_(1+tmp_index_,1),loading_umap_(1+tmp_index_,2),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_lines__(1+nc_lines,:),'MarkerEdgeColor','k');
end;%for nl=0:numel(tmp_index_sort_u_celltype_)-1;
axis equal; axisnotick;
legend(tmp_u_celltype_(1+tmp_index_sort_u_celltype_));
title(sprintf('umap_%s_nlp%.5d <-- %d clusters',infix_method_1,round(1000*nlpnex_threshold),tmp_n_u_celltype),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%;
end;%for nnlpnex_threshold=0:n_nlpnex_threshold-1;
%%%%%%%%%%%%%%%%;
