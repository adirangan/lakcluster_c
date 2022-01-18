function ...
[ ...
 parameter ...
] = ...
test_hnbtZRgumb_multi_collect_6( ...
 parameter ...
,M ...
,N ...
,snr ...
,n_cluster ...
,n_rank_ ...
,n_iteration ...
,n_shuffle ...
,flag_rerun ...
);
% collects results from test_hnbtZRgumb_multi_6.m ;
% test with: ;
% parameter=[]; M = 563; N=6325; snr = 0.80; n_cluster=6; n_rank_ = min(min(M,N),(1+n_cluster) + [0,1,2,4,8,16,32,64,128,256,512,1024]); n_iteration=32; n_shuffle=64; flag_rerun=0;
% parameter=[]; M = 563; N=6325; snr = 0.80; n_cluster=6; n_rank_ = min(min(M,N),(1+n_cluster) + [1,2,4,8,16,32,64,128,256,512]); n_iteration=32; n_shuffle=64; flag_rerun=0;

if (nargin<1);
%%%%%%%%;
platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
for nl=0:0;
na=0;   
if (nl==na); M =  563; N = 6325; n_cluster = 6; snr = 0.80; n_rank_ = min(min(M,N),(1+n_cluster) + [1,2,4,8,16,32,64,128,256,512]); n_iteration = 32; n_shuffle =  64; end; na=na+1;
flag_rerun = 0;
parameter = struct('type','parameter');
parameter = ...
test_hnbtZRgumb_multi_collect_6( ...
 parameter ...
,M ...
,N ...
,snr ...
,n_cluster ...
,n_rank_ ...
,n_iteration ...
,n_shuffle ...
,flag_rerun ...
);
end;%for nl=0:0;
disp('returning'); return;
end;%if (nargin<1);

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); string_root = 'home'; end;
if (strcmp(platform,'eval1')); string_root = 'home'; end;
if (strcmp(platform,'rusty')); string_root = 'mnt/home'; end;

flag_recalc = 0;
flag_replot = 1;
tolerance_master = 1e-2;
nf=0;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root);
str_code = sprintf('%s/halfloop_dev',dir_code);
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jamison',string_root);
dir_base = sprintf('%s/dir_hnbtZRgumb_multi_6',dir_trunk);
dir_jpg = sprintf('%s/dir_jpg',dir_base);
n_rank_sub = 1;

if (isempty(parameter)); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'verbose')); parameter.verbose = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'date_diff_threshold')); parameter.date_diff_threshold = 0.5; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_force_create_mat')); parameter.flag_force_create_mat = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_force_create_tmp')); parameter.flag_force_create_tmp = 0; end; %<-- parameter_bookmark. ;
verbose = parameter.verbose;
date_diff_threshold = parameter.date_diff_threshold;
flag_force_create_mat = parameter.flag_force_create_mat;
flag_force_create_tmp = parameter.flag_force_create_tmp;

if (verbose); disp(sprintf(' %% [entering test_hnbtZRgumb_multi_collect_6]')); end;

if (~exist(dir_base,'dir')); disp(sprintf(' %% mkdir %s',dir_base)); mkdir(dir_base); end;
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

%%%%%%%%;
% global variables for post-processing. ;
%%%%%%%%;
define_global_test_halfloop_driver_6;

flag_plot=0;
if flag_plot;
%%%%%%%%;
% Making specific tree figures for each iteration. ;
%%%%%%%%;
n_n_rank = numel(n_rank_);
for nn_rank=0:n_n_rank-1;
n_rank = n_rank_(1+nn_rank);
for niteration=1:n_iteration;
str_px = test_hnbtZRgumb_multi_xfix_6('px',M,N,snr,n_cluster,n_rank,niteration);
fname_px_mat = sprintf('%s/%s.mat',dir_base,str_px);
if (~exist(fname_px_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',fname_px_mat)); end;
end;%if (~exist(fname_px_mat,'file'));
if ( exist(fname_px_mat,'file'));
if (verbose); disp(sprintf(' %% %s found, not skipping',fname_px_mat)); end;
load(fname_px_mat);
str_dir = test_hnbtZRgumb_multi_xfix_6('dir',M,N,snr,n_cluster,n_rank,niteration);
dir_local_mat = sprintf('%s/%s',dir_base,str_dir);
if (~exist(dir_local_mat,'dir'));
if (verbose); disp(sprintf(' %% %s not found, skipping',dir_local_mat)); end;
end;%if (~exist(dir_local_mat,'dir')); 
if ( exist(dir_local_mat,'dir'));
if (verbose); disp(sprintf(' %% %s found, not skipping',dir_local_mat)); end;
rng(niteration);
[A_n_,label_A_] = random_matrix_planted_cluster_0(M,N,n_rank,n_cluster);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for str_local_method_ = {'hnbr0tZRgumb','hnbrtZRgumb','hnbtZRgumb'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
str_local_method = str_local_method_{1};
if (verbose); disp(sprintf(' %% postprocessing for %s',str_local_method)); end;
try;
%%%%;
if (strcmp(str_local_method,'hnbr0tZRgumb'));
tmp_infix_hnb_0 = sprintf('hnbr0tZRgumb_r%d_r1',n_rank); 
tmp_infix_hnb_1 = sprintf('%s',tmp_infix_hnb_0);
tmp_fname_louvain_mat = sprintf('%s/louvain00pr_default_r%d.mat',dir_local_mat,n_rank);
tmp_fname_umap_mat = sprintf('%s/umap00pr_default_r%d.mat',dir_local_mat,n_rank);
end;%if (strcmp(str_local_method,'hnbr0tZRgumb'));
%%%%;
if (strcmp(str_local_method,'hnbrtZRgumb'));
tmp_infix_hnb_0 = sprintf('hnbrtZRgumb_r%d_r1',n_rank); 
tmp_infix_hnb_1 = sprintf('%s',tmp_infix_hnb_0);
tmp_fname_louvain_mat = sprintf('%s/louvain00pr_default_r%d.mat',dir_local_mat,n_rank);
tmp_fname_umap_mat = sprintf('%s/umap00pr_default_r%d.mat',dir_local_mat,n_rank);
end;%if (strcmp(str_local_method,'hnbrtZRgumb'));
%%%%;
if (strcmp(str_local_method,'hnbtZRgumb'));
tmp_infix_hnb_0 = sprintf('hnbtZRgumb_r1'); 
tmp_infix_hnb_1 = sprintf('%s_g010',tmp_infix_hnb_0);
tmp_fname_louvain_mat = sprintf('%s/louvain00_default_r%d.mat',dir_local_mat,n_rank);
tmp_fname_umap_mat = sprintf('%s/umap00_default_r%d.mat',dir_local_mat,n_rank);
end;%if (strcmp(str_local_method,'hnbtZRgumb'));
%%%%;
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
,A_n_ ...
,n_rank ...
,label_A_ ...
,tmp_fname_louvain_mat ...
,tmp_fname_umap_mat ...
);
catch; disp(sprintf(' %% warning, could not process %s <-- %s',dir_local_mat,str_local_method));
end;%try;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for str_local_method_ = {'hnbr0tZRgumb','hnbrtZRgumb','hnbtZRgumb'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ( exist(dir_local_mat,'dir'));
end;%if ( exist(fname_px_mat,'file'));
end;%for niteration=1:n_iteration;
end;%for nn_rank=0:n_n_rank-1;
%%%%%%%%;
end;%if flag_plot;

%%%%%%%%;
% Collecting statistics across runs. ;
%%%%%%%%;
n_n_rank = numel(n_rank_);
lpv_rmi___ = zeros(n_n_rank,global_n_method,n_iteration);
n_cluster_rmi___ = zeros(n_n_rank,global_n_method,n_iteration);
flag_found_rmi___ = zeros(n_n_rank,global_n_method,n_iteration);
for nn_rank=0:n_n_rank-1;
n_rank = n_rank_(1+nn_rank);
if (verbose); disp(sprintf(' %% nn_rank %d/%d n_rank %0.2f',nn_rank,n_n_rank,n_rank)); end;
for niteration=1:n_iteration;
str_px = test_hnbtZRgumb_multi_xfix_6('px',M,N,snr,n_cluster,0*n_rank,niteration);
fname_px_mat = sprintf('%s/%s.mat',dir_base,str_px);
if (~exist(fname_px_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',fname_px_mat)); end;
end;%if (~exist(fname_px_mat,'file'));
if ( exist(fname_px_mat,'file'));
if (verbose); disp(sprintf(' %% %s found, not skipping',fname_px_mat)); end;
load(fname_px_mat);
str_dir = test_hnbtZRgumb_multi_xfix_6('dir',M,N,snr,n_cluster,0*n_rank,niteration);
dir_local_mat = sprintf('%s/%s',dir_base,str_dir);
if (~exist(dir_local_mat,'dir'));
if (verbose); disp(sprintf(' %% %s not found, skipping',dir_local_mat)); end;
end;%if (~exist(dir_local_mat,'dir')); 
if ( exist(dir_local_mat,'dir'));
if (verbose); disp(sprintf(' %% %s found, not skipping',dir_local_mat)); end;
for nmethod=0:global_n_method-1;
tmp_str_method = sprintf('%s',global_prefix_method_{1+nmethod});
if global_nrank_method_{1+nmethod};
tmp_str_method = sprintf('%s_r%d',global_prefix_method_{1+nmethod},n_rank);
end;%if global_nrank_method_{1+nmethod};
tmp_fname_mat = sprintf('%s/%s.mat',dir_local_mat,tmp_str_method);
if (~exist(tmp_fname_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, skipping',tmp_fname_mat)); end;
end;%if (~exist(tmp_fname_mat,'file'));
if ( exist(tmp_fname_mat,'file'));
if (verbose); disp(sprintf(' %% %s found, not skipping',tmp_fname_mat)); end;
tmp_ = load(tmp_fname_mat);
local_label_B_ = tmp_.label_B__{end};
local_lpv = label_to_label_enrichment_quad_4(label_A_,local_label_B_);
flag_found_rmi___(1+nn_rank,1+nmethod,niteration) = 1;
lpv_rmi___(1+nn_rank,1+nmethod,niteration) = local_lpv;
n_cluster_rmi___(1+nn_rank,1+nmethod,niteration) = numel(unique(local_label_B_));
clear tmp_ local_lpv;
end;%if ( exist(tmp_fname_mat,'file'));
end;%for nmethod=0:global_n_method-1;
end;%if ( exist(dir_local_mat,'dir'));
end;%if ( exist(fname_px_mat,'file'));
end;%for niteration=1:n_iteration;
end;%for nn_rank=0:n_n_rank-1;
%%%%%%%%;

%%%%%%%%;
lpv_avg_rm__ = -ones(n_n_rank,global_n_method);
n_cluster_avg_rm__ = -ones(n_n_rank,global_n_method);
for nn_rank=0:n_n_rank-1;
for nmethod=0:global_n_method-1;
tmp_index_ = efind(flag_found_rmi___(1+nn_rank,1+nmethod,:));
if numel(tmp_index_)>0;
lpv_avg_rm__(1+nn_rank,1+nmethod) = mean(lpv_rmi___(1+nn_rank,1+nmethod,1+tmp_index_));
n_cluster_avg_rm__(1+nn_rank,1+nmethod) = mean(n_cluster_rmi___(1+nn_rank,1+nmethod,1+tmp_index_));
end;%if numel(tmp_index_)>0;
end;%for nmethod=0:global_n_method-1;
end;%for nn_rank=0:n_n_rank-1;
%%%%%%%%;
prctile_ = [ 5 , 15 , 50 , 85 , 95 ]; n_prctile = numel(prctile_);
lpv_prctile_prm___ = -ones(n_prctile,n_n_rank,global_n_method);
n_cluster_prctile_prm___ = -ones(n_prctile,n_n_rank,global_n_method);
for nn_rank=0:n_n_rank-1;
for nmethod=0:global_n_method-1;
tmp_index_ = efind(flag_found_rmi___(1+nn_rank,1+nmethod,:));
if numel(tmp_index_)>0;
lpv_prctile_prm___(:,1+nn_rank,1+nmethod) = prctile(lpv_rmi___(1+nn_rank,1+nmethod,1+tmp_index_),prctile_);
n_cluster_prctile_prm___(:,1+nn_rank,1+nmethod) = prctile(n_cluster_rmi___(1+nn_rank,1+nmethod,1+tmp_index_),prctile_);
end;%if numel(tmp_index_)>0;
end;%for nmethod=0:global_n_method-1;
end;%for nn_rank=0:n_n_rank-1;
%%%%%%%%;

%%%%%%%%;
% Making general figures. ;
%%%%%%%%;
n_rank=0;
log2_n_rank_ = log2(n_rank_ - (1+n_cluster)); %log2_n_rank_(1) = -1;
str_fig = test_hnbtZRgumb_multi_xfix_6('fig_trialmean',M,N,snr,n_cluster,n_rank,niteration);
fname_fig = sprintf('%s/%s',dir_jpg,str_fig);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;
markersize_use = 16;
%%%%;
subplot(1,2,1);
hold on;
for nl=0:numel(global_index_method_proj_0on_)-1;
index_method = global_index_method_proj_0on_(1+nl);
symbol_method = global_symbol_method_{1+index_method};
plot(log2_n_rank_,-lpv_avg_rm__(:,1+index_method),symbol_method,'LineWidth',2,'MarkerSize',markersize_use);
end;%for nl=0:numel(global_index_method_proj_0on_)-1;
hold off;
xlim([min(log2_n_rank_),max(log2_n_rank_)]);
xlabel('log2(n_rank - (1+n_cluster))','Interpreter','none');
ylabel('-log(P)');
legend(global_legend_method_(1+global_index_method_proj_0on_),'Location','NorthEast');
title('proj_0on','Interpreter','none');
%%%%;
subplot(1,2,2);
hold on;
for nl=0:numel(global_index_method_proj_0on_)-1;
index_method = global_index_method_proj_0on_(1+nl);
symbol_method = global_symbol_method_{1+index_method};
plot(log2_n_rank_,n_cluster_avg_rm__(:,1+index_method),symbol_method,'LineWidth',2,'MarkerSize',markersize_use);
end;%for nl=0:numel(global_index_method_proj_0on_)-1;
hold off;
xlim([min(log2_n_rank_),max(log2_n_rank_)]);
xlabel('log2(n_rank - (1+n_cluster))','Interpreter','none');
ylabel('n_cluster','Interpreter','none');
legend(global_legend_method_(1+global_index_method_proj_0on_),'Location','NorthEast');
title('proj_0on','Interpreter','none');
%%%%;
sgtitle(str_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

flag_plot=0;
if flag_plot;
n_rank=0;
log2_n_rank_ = log2(n_rank_ - (1+n_cluster)); %log2_n_rank_(1) = -1;
str_fig = test_hnbtZRgumb_multi_xfix_6('fig_trialprctl',M,N,snr,n_cluster,n_rank,niteration);
fname_fig = sprintf('%s/%s',dir_jpg,str_fig);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;
markersize_use = 16;
dr = mean(diff(log2_n_rank_));
%%%%;
tmp_index_ = intersect(global_index_method_sub_,global_index_method_proj_0on_);
tmp_drnl_ = linspace(-0.25*dr,+0.25*dr,numel(tmp_index_));
%%%%;
subplot(1,2,1);
hold on;
for nl=0:numel(tmp_index_)-1;
tmp_drnl = tmp_drnl_(1+nl);
index_method = tmp_index_(1+nl);
symbol_method = global_symbol_method_{1+index_method};
tmp_ymid_ = -lpv_prctile_prm___(3,:,1+index_method);
tmp_yneg_ = -lpv_prctile_prm___(4,:,1+index_method); tmp_yneg_ = tmp_ymid_ - tmp_yneg_;
tmp_ypos_ = -lpv_prctile_prm___(2,:,1+index_method); tmp_ypos_ = tmp_ypos_ - tmp_ymid_;
errorbar(log2_n_rank_ + tmp_drnl,tmp_ymid_,tmp_yneg_,tmp_ypos_,symbol_method(1:2),'LineWidth',4,'MarkerSize',markersize_use);
tmp_ymid_ = -lpv_prctile_prm___(3,:,1+index_method);
tmp_yneg_ = -lpv_prctile_prm___(5,:,1+index_method); tmp_yneg_ = tmp_ymid_ - tmp_yneg_;
tmp_ypos_ = -lpv_prctile_prm___(1,:,1+index_method); tmp_ypos_ = tmp_ypos_ - tmp_ymid_;
errorbar(log2_n_rank_ + tmp_drnl,tmp_ymid_,tmp_yneg_,tmp_ypos_,symbol_method(1:2),'LineWidth',2,'MarkerSize',markersize_use);
end;%for nl=0:numel(tmp_index_)-1;
hold off;
xlim([min(log2_n_rank_),max(log2_n_rank_)]);
xlabel('log2(n_rank - (1+n_cluster))','Interpreter','none');
ylabel('-log(P)');
legend(global_legend_method_(1+tmp_index_),'Location','NorthWest');
title('proj_0on','Interpreter','none');
%%%%;
subplot(1,2,2);
hold on;
for nl=0:numel(tmp_index_)-1;
tmp_drnl = tmp_drnl_(1+nl);
index_method = tmp_index_(1+nl);
symbol_method = global_symbol_method_{1+index_method};
tmp_ymid_ = n_cluster_prctile_prm___(3,:,1+index_method);
tmp_yneg_ = n_cluster_prctile_prm___(2,:,1+index_method); tmp_yneg_ = tmp_ymid_ - tmp_yneg_;
tmp_ypos_ = n_cluster_prctile_prm___(4,:,1+index_method); tmp_ypos_ = tmp_ypos_ - tmp_ymid_;
errorbar(log2_n_rank_ + tmp_drnl,tmp_ymid_,tmp_yneg_,tmp_ypos_,symbol_method(1:2),'LineWidth',4,'MarkerSize',markersize_use);
tmp_ymid_ = n_cluster_prctile_prm___(3,:,1+index_method);
tmp_yneg_ = n_cluster_prctile_prm___(1,:,1+index_method); tmp_yneg_ = tmp_ymid_ - tmp_yneg_;
tmp_ypos_ = n_cluster_prctile_prm___(5,:,1+index_method); tmp_ypos_ = tmp_ypos_ - tmp_ymid_;
errorbar(log2_n_rank_ + tmp_drnl,tmp_ymid_,tmp_yneg_,tmp_ypos_,symbol_method(1:2),'LineWidth',2,'MarkerSize',markersize_use);
end;%for nl=0:numel(tmp_index_)-1;
hold off;
xlim([min(log2_n_rank_),max(log2_n_rank_)]);
xlabel('log2(n_rank - (1+n_cluster))','Interpreter','none');
ylabel('n_cluster','Interpreter','none');
legend(global_legend_method_(1+tmp_index_),'Location','NorthWest');
title('proj_0on','Interpreter','none');
%%%%;
sgtitle(str_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

if (verbose); disp(sprintf(' %% [finished test_hnbtZRgumb_multi_collect_6]')); end;
