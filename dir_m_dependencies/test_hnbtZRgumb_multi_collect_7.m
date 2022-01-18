function ...
[ ...
 parameter ...
] = ...
test_hnbtZRgumb_multi_collect_7( ...
 parameter ...
,M ...
,N ...
,snr_ ...
,n_cluster ...
,n_rank_ ...
,n_iteration ...
,n_shuffle ...
,flag_rerun ...
);
% collects results from test_hnbtZRgumb_multi_7.m ;
% test with: ;
% parameter=[]; M = 563; N=6325; n_step = 21; snr_ = linspace(0.5,1.0,n_step); n_cluster=6; n_rank_ = min(min(M,N),(1+n_cluster) + [1,2,4,8,16,32,64,128,256,512]); n_iteration = 16; n_shuffle=64; flag_rerun=0;

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
if (nl==na); M =  563; N = 6325; n_step = 21; snr_ = linspace(0.5,1.0,n_step); n_cluster = 6; n_rank_ = min(min(M,N),(1+n_cluster) + [1,2,4,8,16,32,64,128,256,512]); n_iteration = 16; n_shuffle =  64; end; na=na+1;
flag_rerun = 0;
parameter = struct('type','parameter');
parameter = ...
test_hnbtZRgumb_multi_collect_7( ...
 parameter ...
,M ...
,N ...
,snr_ ...
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
dir_base = sprintf('%s/dir_hnbtZRgumb_multi_7',dir_trunk);
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

if (verbose); disp(sprintf(' %% [entering test_hnbtZRgumb_multi_collect_7]')); end;

if (~exist(dir_base,'dir')); disp(sprintf(' %% mkdir %s',dir_base)); mkdir(dir_base); end;
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

%%%%%%%%;
% global variables for post-processing. ;
%%%%%%%%;
define_global_test_halfloop_driver_6;

%%%%%%%%;
% Collecting statistics across runs. ;
%%%%%%%%;
n_snr = numel(snr_);
n_n_rank = numel(n_rank_);
lpv_srmi____ = zeros(n_snr,n_n_rank,global_n_method,n_iteration);
n_cluster_srmi____ = zeros(n_snr,n_n_rank,global_n_method,n_iteration);
flag_found_srmi____ = zeros(n_snr,n_n_rank,global_n_method,n_iteration);
for nsnr=0:n_snr-1;
snr = snr_(1+nsnr);
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
flag_found_srmi____(1+nsnr,1+nn_rank,1+nmethod,niteration) = 1;
lpv_srmi____(1+nsnr,1+nn_rank,1+nmethod,niteration) = local_lpv;
n_cluster_srmi____(1+nsnr,1+nn_rank,1+nmethod,niteration) = numel(unique(local_label_B_));
clear tmp_ local_lpv;
end;%if ( exist(tmp_fname_mat,'file'));
end;%for nmethod=0:global_n_method-1;
end;%if ( exist(dir_local_mat,'dir'));
end;%if ( exist(fname_px_mat,'file'));
end;%for niteration=1:n_iteration;
end;%for nn_rank=0:n_n_rank-1;
end;%for nsnr=0:n_snr-1;
%%%%%%%%;

%%%%%%%%;
lpv_avg_srm___ = -ones(n_snr,n_n_rank,global_n_method);
n_cluster_avg_srm___ = -ones(n_snr,n_n_rank,global_n_method);
for nsnr=0:n_snr-1;
for nn_rank=0:n_n_rank-1;
for nmethod=0:global_n_method-1;
tmp_index_ = efind(flag_found_srmi____(1+nsnr,1+nn_rank,1+nmethod,:));
if numel(tmp_index_)>0;
lpv_avg_srm___(1+nsnr,1+nn_rank,1+nmethod) = mean(lpv_srmi____(1+nsnr,1+nn_rank,1+nmethod,1+tmp_index_));
n_cluster_avg_srm___(1+nsnr,1+nn_rank,1+nmethod) = mean(n_cluster_srmi____(1+nsnr,1+nn_rank,1+nmethod,1+tmp_index_));
end;%if numel(tmp_index_)>0;
end;%for nmethod=0:global_n_method-1;
end;%for nn_rank=0:n_n_rank-1;
end;%for nsnr=0:n_snr-1;
%%%%%%%%;

%%%%%%%%;
% Making general figures. ;
%%%%%%%%;
snr=0;
n_rank=0;
log2_n_rank_ = log2(n_rank_ - (1+n_cluster)); %log2_n_rank_(1) = -1;
str_fig = test_hnbtZRgumb_multi_xfix_6('fig_trialmean',M,N,snr,n_cluster,n_rank,niteration);
fname_fig = sprintf('%s/%s',dir_jpg,str_fig);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;figbeach();
%%%%;
local_index_method_sub_ = ...
[ ...
,efind(cellfun(@(x)strcmp(x,'umap00pr_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'louvain00pr_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'hnbrtZRgumb'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'hnbr0tZRgumb'),global_prefix_method_)) ...
];
%%%%;
for nl=0:numel(local_index_method_sub_)-1;
subplot(2,numel(local_index_method_sub_),1+nl+0*numel(local_index_method_sub_));
index_method = local_index_method_sub_(1+nl);
prefix_method = global_prefix_method_{1+index_method};
imagesc(squeeze(-lpv_avg_srm___(:,:,1+index_method)),[0,550]);
set(gca,'YTick',1:n_snr,'YTickLabel',snr_); set(gca,'Ydir','normal');
set(gca,'XTick',1:n_n_rank,'XTickLabel',log2_n_rank_);
ylabel('snr'); xlabel('log2(n_rank - (1+n_cluster))','Interpreter','none');
title(sprintf('%s: lpv',prefix_method),'Interpreter','none');
colorbar;
subplot(2,numel(local_index_method_sub_),1+nl+1*numel(local_index_method_sub_));
index_method = local_index_method_sub_(1+nl);
prefix_method = global_prefix_method_{1+index_method};
imagesc(squeeze(n_cluster_avg_srm___(:,:,1+index_method)),[0,2*(1+n_cluster)]);
set(gca,'YTick',1:n_snr,'YTickLabel',snr_); set(gca,'Ydir','normal');
set(gca,'XTick',1:n_n_rank,'XTickLabel',log2_n_rank_);
ylabel('snr'); xlabel('log2(n_rank - (1+n_cluster))','Interpreter','none');
title(sprintf('%s: n_cluster',prefix_method),'Interpreter','none');
colorbar;
end;%for nl=0:numel(local_index_method_sub_)-1;
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
if (verbose); disp(sprintf(' %% [finished test_hnbtZRgumb_multi_collect_7]')); end;
