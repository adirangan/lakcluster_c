function analyze_10xgenomics_5(dir_data_infix);
%%%%%%%%;
% 
%%%%%%%%;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

verbose=1;
flag_recalc = 0;
flag_replot = 0;
tolerance_master = 1e-2;
nf=0;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root);
str_code = sprintf('%s/halfloop_dev',dir_code);
dir_data = sprintf('/%s/rangan/dir_bcc/dir_10xgenomics/dir_%s',string_root,dir_data_infix);
dir_mat = sprintf('%s/dir_mat',dir_data); if ~exist(dir_mat); disp(sprintf(' %% mkdir %s',dir_mat)); mkdir(sprintf('%s',dir_mat)); end;
dir_jpg = sprintf('%s/dir_jpg',dir_data); if ~exist(dir_jpg); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(sprintf('%s',dir_jpg)); end;
dir_local = sprintf('%s/dir_local',dir_data); if (~exist(dir_local,'dir')); disp(sprintf(' %% mkdir %s',dir_local)); mkdir(dir_local); end;

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
[n_umi,n_smp] = size(A_gs__);
nnz_each_ = A_gs__(:);
nnz_umi_ = sum(A_gs__,2);
nnz_smp_ = sum(A_gs__,1);

%%%%%%%%;
% Now we could fit nnz_umi_ with an algebraic. ;
% assumed form is given by: ;
% p(n;alpha,j) = 1/(n+alpha)^(1+j). ;
%%%%%%%%;
% n_N = 64; plot(log(1+[0:n_N]),log(hist(nnz_umi_,[0:n_N])),'ro',log([0:n_N]),log(1./[0:n_N].^1.2) + 8.5,'k-');

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
fname_fig = sprintf('%s/nnz_umi',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
n_N = 64;
subplot(1,2,1); semilogy(hist(nnz_umi_,[0:n_N]),'ko-'); xlim([-1,1+n_N]); xlabel('n'); ylabel('hist');
subplot(1,2,2); loglog(hist(nnz_umi_,[0:n_N]),'ko-');  xlabel('n'); ylabel('hist');
sgtitle('nnz_umi','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/nnz_smp',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
h_nnz_ = mean(nnz_smp_) + 4.5*std(nnz_smp_,1)*linspace(-1,1,128);
hist(nnz_smp_,h_nnz_);
xlabel('n'); ylabel('hist');
sgtitle('nnz_smp','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

nnz_umi_ = sum(A_gs__~=0,2);
nnz_threshold_ = [1,2,round(n_smp*0.05),round(n_smp*0.10),round(n_smp*0.25),round(n_smp*0.50),round(n_smp*0.75)];
n_nnz_threshold = numel(nnz_threshold_);

%%%%%%%%;
% global variables for post-processing. ;
%%%%%%%%;
global_prefix_method_ = { ...
 'tsne00_isosplit5' ...
,'tsne50_isosplit5' ...
,'umap00_default' ...
,'umap00_isosplit5' ...
,'umap00_hdbscan' ...
,'louvain00_default' ...
,'hnbtZRgumb' ...
,'spectral_isosplit5' ...
,'tsne00pr_isosplit5' ...
,'tsne50pr_isosplit5' ...
,'umap00pr_default' ...
,'umap00pr_isosplit5' ...
,'umap00pr_hdbscan' ...
,'louvain00pr_default' ...
,'hnbrtZRgumb' ...
};
global_n_method = numel(global_prefix_method_);
global_index_method_proj_off_ = [0,1,2,3,4,5,6];
global_index_method_proj_0on_ = [7,8,9,10,11,12,13,14];
global_index_method_sub_ = [2,5,6,10,13,14];
%%%%%%%%;
global_legend_method_ = { ...
 'tsne00_isosplit5' ...
,'tsne50_isosplit5' ...
,'umap00_default' ...
,'umap00_isosplit5' ...
,'umap00_hdbscan' ...
,'louvain00_default' ...
,'halfloop00_c' ...
,'spectral_isosplit5' ...
,'tsne00pr_isosplit5' ...
,'tsne50pr_isosplit5' ...
,'umap00pr_default' ...
,'umap00pr_isosplit5' ...
,'umap00pr_hdbscan' ...
,'louvain00pr_default' ...
,'halfloop00pr_c' ...
};
for nl=0:numel(global_legend_method_)-1; global_legend_method_{1+nl}(strfind(global_legend_method_{1+nl},'_')) = ' '; end;
global_nrank_method_ = { ...
 [] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,1 ...
,1 ...
,1 ...
,1 ...
,1 ...
,1 ...
,1 ...
,1 ...
};
global_symbol_method_ = { ...
 'co-' ...
,'bo-' ...
,'kx-' ...
,'ko-' ...
,'k^-' ...
,'rs-' ...
,'mh-' ...
,'go-' ...
,'co-' ...
,'bo-' ...
,'kx-' ...
,'ko-' ...
,'k^-' ...
,'rs-' ...
,'mh-' ...
};
%%%%%%%%;

n_svd = 36;
date_diff_threshold = 0.5;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nnnz_threshold=n_nnz_threshold-1:-1:0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
nnz_threshold = nnz_threshold_(1+nnnz_threshold);
local_prefix = sprintf('nnzt%.4d',nnz_threshold);
local_index_umi_retain_ = efind(nnz_umi_> nnz_threshold);
local_n_smp = n_smp; local_n_umi = numel(local_index_umi_retain_);
local_B_logn_sg__ = log(1+transpose(A_gs__(1+local_index_umi_retain_,:)));
%%%%%%%%;
dir_local_mat = sprintf('%s/dir_%s_mat',dir_local,local_prefix);
if (~exist(dir_local_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_local_mat)); mkdir(dir_local_mat); end;
dir_local_jpg = sprintf('%s/dir_%s_jpg',dir_local,local_prefix);
if (~exist(dir_local_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_local_jpg)); mkdir(dir_local_jpg); end;
%%%%%%%%;
local_fname_pre = sprintf('%s/S_%s_B_logn_sg__',dir_local_mat,local_prefix);
[local_flag_skip,local_fname_mat] = open_fname_tmp(local_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~local_flag_skip;
[U_local_B_logn_sg__,S_local_B_logn_sg__,V_local_B_logn_sg__] = svds(local_B_logn_sg__,n_svd); S_local_B_logn_sg_ = diag(S_local_B_logn_sg__);
save(local_fname_mat ...
     ,'n_svd' ...
     ,'local_prefix','local_index_umi_retain_' ...
     ,'local_n_umi' ...
     ,'U_local_B_logn_sg__' ...
     ,'S_local_B_logn_sg_' ...
     ,'V_local_B_logn_sg__' ...
     );
close_fname_tmp(local_fname_pre)
end;%if ~local_flag_skip;
%%%%%%%%;
local_fname_shuffle_pre = sprintf('%s/S_%s_B_logn_sg_shuffle__',dir_local_mat,local_prefix);
[local_flag_shuffle_skip,local_fname_shuffle_mat] = open_fname_tmp(local_fname_shuffle_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~local_flag_shuffle_skip;
n_shuffle = 32;
S_local_B_logn_sg_shuffle__ = zeros(n_svd,n_shuffle);
for nshuffle=0:n_shuffle-1;
rng(nshuffle*1024);
tmp_local_B_logn_sg__ = local_B_logn_sg__;
for local_numi=0:local_n_umi-1;
if (mod(local_numi,1024)==0); disp(sprintf(' %% nshuffle %d/%d, local_numi %d/%d',nshuffle,n_shuffle,local_numi,local_n_umi)); end;
tmp_ij_ = randperm(local_n_smp);
tmp_local_B_logn_sg__(:,1+local_numi) = tmp_local_B_logn_sg__(tmp_ij_,1+local_numi);
end;%for local_numi=0:local_n_umi-1;
[tmp_U_local_B_logn_sg__,tmp_S_local_B_logn_sg__,tmp_V_local_B_logn_sg__] = svds(tmp_local_B_logn_sg__,n_svd); tmp_S_local_B_logn_sg_ = diag(tmp_S_local_B_logn_sg__);
S_local_B_logn_sg_shuffle__(:,1+nshuffle) = tmp_S_local_B_logn_sg_;
end;%for nshuffle=0:n_shuffle-1;
save(local_fname_shuffle_mat ...
     ,'n_svd' ...
     ,'n_shuffle' ...
     ,'S_local_B_logn_sg_shuffle__' ...
     );
close_fname_tmp(local_fname_shuffle_pre)
end;%if ~local_flag_shuffle_skip;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (  exist(local_fname_mat,'file') &  exist(local_fname_shuffle_mat,'file') );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
load(local_fname_mat); load(local_fname_shuffle_mat);
fname_fig = sprintf('%s/S_%s_B_logn_sg_FIGA',dir_local_jpg,local_prefix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;
semilogy(1:n_svd,S_local_B_logn_sg_,'kx',1:n_svd,S_local_B_logn_sg_shuffle__,'ro');
title(sprintf('S_%s_B_logn_',local_prefix),'Interpreter','none');
xlabel('rank'); ylabel('sigma');
xlim([0.5,n_svd+0.5]);
grid on;
set(gcf,'Position',1+[0,0,512,1024]);
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
local_rank_estimate_B_logn_sg = max(1+efind(S_local_B_logn_sg_ > prctile(S_local_B_logn_sg_shuffle__,95,2)));
if (verbose); disp(sprintf(' %% %s: local_rank_estimate_B_logn_sg %d/%d',local_prefix,local_rank_estimate_B_logn_sg,n_svd)); end;
%%%%%%%%;
label_xxxx_ = azimuth_pred_celltype_;
[ ...
 label_xxxx_enum_ ...
,n_u_label_xxxx ...
,u_label_xxxx_ ...
,index_nu_xxxx_from_nall_ ...
,n_u_label_xxxx_ ...
,index_nall_from_nu_xxxx__ ...
] = ...
label_str_to_enum_1( ...
 label_xxxx_ ...
);
%%%%%%%%;
[~,index_sort_u_label_xxxx_] = sort(n_u_label_xxxx_,'descend'); index_sort_u_label_xxxx_ = index_sort_u_label_xxxx_ - 1;
index_nu_srt_from_nu_ori_ = index_sort_u_label_xxxx_;
n_all = sum(n_u_label_xxxx_);
for nu_label_xxxx=0:n_u_label_xxxx-1;
tmp_n_u = n_u_label_xxxx_(1+index_sort_u_label_xxxx_(1+nu_label_xxxx));
u_label_xxxx = u_label_xxxx_{1+index_sort_u_label_xxxx_(1+nu_label_xxxx)};
disp(sprintf(' %% nu %.2d/%.2d: %.4d/%.4d = %0.2f = %s',nu_label_xxxx,n_u_label_xxxx,tmp_n_u,n_all,tmp_n_u/n_all,u_label_xxxx));
end;%for nu_label_xxxx=0:n_u_label_xxxx-1;
%%%%%%%%;
fname_fig = sprintf('%s/label_%s_B_logn_sg_pca_FIGA',dir_local_jpg,local_prefix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
try;
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
markersize_final = 12;
markersize_start = 4;
str_symbol_ = {'o','^','s','p','h'};
c_lines__ = colormap('lines'); n_c_lines = size(c_lines__,1);
hold on;
for nl=0:numel(index_sort_u_label_xxxx_)-1;
nc_lines = nl;
str_symbol = str_symbol_{1+mod(nl,5)};
markersize_use = round(markersize_start + (markersize_final-markersize_start)*nl/(numel(index_sort_u_label_xxxx_)-1));
tmp_nu_label_xxxx = index_sort_u_label_xxxx_(1+nl);
tmp_index_ = index_nall_from_nu_xxxx__{1+tmp_nu_label_xxxx};
plot(U_local_B_logn_sg__(1+tmp_index_,1),U_local_B_logn_sg__(1+tmp_index_,2),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_lines__(1+nc_lines,:),'MarkerEdgeColor','k');
end;%for nl=0:numel(index_sort_u_label_xxxx_)-1;
axis equal; axisnotick;
tmp_u_label_xxxx_ = u_label_xxxx_;
for nu_label_xxxx=0:n_u_label_xxxx-1;
tmp_u_label_xxxx_{1+nu_label_xxxx}(strfind(tmp_u_label_xxxx_{1+nu_label_xxxx},'_'))=' ';
end;%for nu_label_xxxx=0:n_u_label_xxxx-1;
legend(tmp_u_label_xxxx_(1+index_sort_u_label_xxxx_));
title('pca');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
catch; disp(sprintf(' %% Warning, could not create %s',fname_fig)); end;%try;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% cluster local data-set. ;
%%%%%%%%;
flag_calculate = strcmp(platform,'access1') | strcmp(platform,'eval1') | strcmp(platform,'OptiPlex');
parameter = struct('type','parameter');
parameter.verbose = verbose;
parameter.date_diff_threshold = date_diff_threshold;
parameter.flag_force_create_mat = flag_force_create_mat;
parameter.flag_force_create_tmp = flag_force_create_tmp;
parameter.str_code = str_code;
parameter.halfloop_recursion_limit = 32;
parameter.flag_skip_hnbtZRgumb = 1;
parameter.flag_skip_hnbrtZRgumb = 0;
if flag_calculate;
test_halfloop_driver_6( ...
 parameter ...
,local_B_logn_sg__ ...
,index_nu_xxxx_from_nall_ ...
,dir_local_mat ...
,local_rank_estimate_B_logn_sg ...
);
end;%if flag_calculate;
%%%%%%%%;
% Now postprocessing. ;
%%%%%%%%;

flag_replot=1;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now post-process local clustering results. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
%%%%%%%;
local_flag_found_ = zeros(global_n_method,1);
local_n_cluster_found_ = zeros(global_n_method,1);
local_label_B__ = cell(global_n_method,1);
local_lpv_xxxx_ = zeros(global_n_method,1);
%%%%%%%%;
for nmethod=0:global_n_method-1;
local_prefix_method = global_prefix_method_{1+nmethod};
global_nrank_method = global_nrank_method_{1+nmethod};
if (~isempty(global_nrank_method)); local_prefix_method = sprintf('%s_r%d',local_prefix_method,local_rank_estimate_B_logn_sg); end;
local_fname_mat = sprintf('%s/%s.mat',dir_local_mat,local_prefix_method);
if (~exist(local_fname_mat,'file'));
local_flag_found = 0;
local_n_cluster_found = 0;
local_lpv_xxxx = 0;
local_label_B_ = {};
if (verbose>0); disp(sprintf(' %% %s not found, skipping',local_fname_mat)); end;
end;%if (~exist(local_fname_mat,'file'));
if ( exist(local_fname_mat,'file'));
tmp_ = load(local_fname_mat,'lpv_','lP0_','fla_','label_B__');
local_flag_found = 1;
local_n_cluster_found = numel(unique(tmp_.label_B__{1}));
local_label_B_ = tmp_.label_B__{end};
local_lpv_xxxx = label_to_label_enrichment_quad_4(label_xxxx_enum_,local_label_B_);
end;%if ( exist(fname_mat,'file'));
local_flag_found_(1+nmethod) = local_flag_found;
local_n_cluster_found_(1+nmethod) = local_n_cluster_found;
local_lpv_xxxx_(1+nmethod) = local_lpv_xxxx;
local_label_B__{1+nmethod} = local_label_B_;
end;%for nmethod=0:global_n_method-1;
%%%%%%%%;
if (verbose);
disp(sprintf(' %% -------------------------------------------------------------------------------------------------------------------------------- '));
for nl=0:numel(global_index_method_sub_)-1;
nmethod = global_index_method_sub_(1+nl);
disp(sprintf(' %% %21s (%.3d): \t all %+12.2f ',global_prefix_method_{1+nmethod},local_n_cluster_found_(1+nmethod),local_lpv_xxxx_(1+nmethod)));
end;%for nl=0:numel(global_index_method_sub_)-1;
end;%if (verbose);
%%%%%%%%;
% loading local hnbrtZRgumb and hnbtZRgumb results. ;
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for str_local_method_ = {'hnbrtZRgumb'};%for str_local_method_ = {'hnbrtZRgumb','hnbtZRgumb'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
str_local_method = str_local_method_{1};
if (verbose); disp(sprintf(' %% postprocessing for %s',str_local_method)); end;
if (strcmp(str_local_method,'hnbrtZRgumb'));
tmp_infix_hnb_0 = sprintf('hnbrtZRgumb_r%d_r1',local_rank_estimate_B_logn_sg); 
tmp_infix_hnb_1 = sprintf('%s',tmp_infix_hnb_0);
tmp_fname_louvain_mat = sprintf('%s/louvain00pr_default_r%d.mat',dir_local_mat,local_rank_estimate_B_logn_sg);
tmp_fname_umap_mat = sprintf('%s/umap00pr_default_r%d.mat',dir_local_mat,local_rank_estimate_B_logn_sg);
end;%if (strcmp(str_local_method,'hnbrtZRgumb'));
if (strcmp(str_local_method,'hnbtZRgumb'));
tmp_infix_hnb_0 = sprintf('hnbtZRgumb_r1'); 
tmp_infix_hnb_1 = sprintf('%s_g010',tmp_infix_hnb_0);
tmp_fname_louvain_mat = sprintf('%s/louvain00_default_r%d.mat',dir_local_mat,local_rank_estimate_B_logn_sg);
tmp_fname_umap_mat = sprintf('%s/umap00_default_r%d.mat',dir_local_mat,local_rank_estimate_B_logn_sg);
end;%if (strcmp(str_local_method,'hnbtZRgumb'));
str_infix = tmp_infix_hnb_1;
dir_local_mat_hnb = sprintf('%s/dir_tmp_%s/dir_tmp_%s',dir_local_mat,tmp_infix_hnb_0,tmp_infix_hnb_1);
str_sgtitle = dir_local_mat_hnb;
tmp_local_fname_output_label = sprintf('%s/output_label__.txt',dir_local_mat_hnb);
tmp_local_fname_nlpbra_label = sprintf('%s/nlpbra_label__.txt',dir_local_mat_hnb);
tmp_local_fname_nlpnex_label = sprintf('%s/nlpnex_label__.txt',dir_local_mat_hnb);
%%%%;
parameter = struct('type','parameter');
parameter.flag_force_create_mat = 0;
parameter.flag_force_create_tmp = 0;
parameter.flag_replot = flag_replot;
label_tree_compilation_plot_0( ...
 parameter ...
,dir_local_mat ...
,str_infix ...
,str_sgtitle ...
,tmp_local_fname_output_label ...
,tmp_local_fname_nlpbra_label ...
,tmp_local_fname_nlpnex_label ...
,local_B_logn_sg__ ...
,local_rank_estimate_B_logn_sg ...
,label_xxxx_enum_ ...
,tmp_fname_louvain_mat ...
,tmp_fname_umap_mat ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for str_local_method_ = {'hnbrtZRgumb','hnbtZRgumb'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Finished post-processing. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
flag_replot=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (  exist(local_fname_mat,'file') &  exist(local_fname_shuffle_mat,'file') );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nnnz_threshold=0:n_nnz_threshold-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
