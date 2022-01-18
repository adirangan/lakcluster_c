function ...
[ ...
 parameter ...
] = ...
test_hnbtZRgumb_multi_7( ...
 parameter ...
,M ...
,N ...
,snr ...
,n_cluster ...
,n_rank ...
,n_iteration ...
,n_shuffle ...
,flag_rerun ...
);
% tries to find up to 8 clusters. ;
% single rank. ;
% see test_hnbtZRgumb_multi_7_collect.m ;

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
if (nl==na); M =  563; N = 6325; n_cluster = 6; n_iteration = 16; n_shuffle =  64; end; na=na+1;
flag_rerun = 0;
n_rank_ = min(min(M,N),(1+n_cluster) + [1,2,4,8,16,32,64,128,256,512]);
n_n_rank = numel(n_rank_);
n_step = 21;
snr_ = linspace(0.5,1.0,n_step);
for nn_rank=0:n_n_rank-1;
n_rank = n_rank_(1+nn_rank);
for nstep=1:n_step;
snr = snr_(nstep);
parameter = struct('type','parameter');
parameter = ...
test_hnbtZRgumb_multi_7( ...
 parameter ...
,M ...
,N ...
,snr ...
,n_cluster ...
,n_rank ...
,n_iteration ...
,n_shuffle ...
,flag_rerun ...
);
end;%for nn_rank=0:n_n_rank-1;
end;%for nstep=1:n_step;
end;%for nl=0:0;
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); string_root = 'home'; end;
if (strcmp(platform,'eval1')); string_root = 'home'; end;
if (strcmp(platform,'rusty')); string_root = 'mnt/home'; end;

flag_recalc = 0;
flag_replot = 0;
tolerance_master = 1e-2;
nf=0;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root);
str_code = sprintf('%s/halfloop_dev',dir_code);
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_jamison',string_root);
dir_base = sprintf('%s/dir_hnbtZRgumb_multi_7',dir_trunk);
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

if (verbose); disp(sprintf(' %% [entering test_hnbtZRgumb_multi_7]')); end;

if (~exist(dir_base,'dir')); disp(sprintf(' %% mkdir %s',dir_base)); mkdir(dir_base); end;
%%%%%%%%;
for niteration=1:n_iteration;
disp(sprintf(' %% niteration %d/%d',niteration,n_iteration));
%%%%%%%%;
A_n_ = [];
str_px = test_hnbtZRgumb_multi_xfix_6('px',M,N,snr,n_cluster,0*n_rank,niteration);
fname_px_mat = sprintf('%s/%s.mat',dir_base,str_px);
if (~exist(fname_px_mat,'file'));
if (isempty(A_n_));
rng(niteration); [A_n_,label_A_,n_label_A_,pf_,pi_,snr_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);
end;%if (isempty(A_n_));
disp(sprintf(' %% %s not found, creating',fname_px_mat));
save(fname_px_mat,'M','N','snr','label_A_','n_label_A_','pf_','pi_','n_cluster');
end;%if (~exist(fname_px_mat,'file'));
if ( exist(fname_px_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_px_mat));
load(fname_px_mat);
end;%if ( exist(fname_px_mat,'file'));

str_dir = test_hnbtZRgumb_multi_xfix_6('dir',M,N,snr,n_cluster,0*n_rank,niteration);
dir_local_mat = sprintf('%s/%s',dir_base,str_dir);
if (~exist(dir_local_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_local_mat)); mkdir(dir_local_mat); end;

flag_calculate = strcmp(platform,'access1') | strcmp(platform,'eval1') | strcmp(platform,'OptiPlex');
if flag_calculate;
if (isempty(A_n_));
rng(niteration); [A_n_,label_A_,n_label_A_,pf_,pi_,snr_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);
end;%if (isempty(A_n_));
parameter = struct('type','parameter');
parameter.verbose = verbose;
parameter.date_diff_threshold = date_diff_threshold;
parameter.flag_force_create_mat = flag_force_create_mat;
parameter.flag_force_create_tmp = flag_force_create_tmp;
parameter.str_code = str_code;
parameter.halfloop_recursion_limit = 32;
parameter.flag_skip_hnbtZRgumb = 1;
parameter.flag_skip_hnbrtZRgumb = 0;
parameter.flag_skip_tsne = 1;
parameter.flag_skip_spectral = 1;
parameter.flag_skip_nondefault = 1;
parameter.n_shuffle = n_shuffle;
parameter.flag_orth_brute = 0;
if (M< 1024);
parameter.flag_orth_brute = 1; %<-- sample from full space of rotations. ;
end;%if (M< 1024);
test_halfloop_driver_6( ...
 parameter ...
,A_n_ ...
,label_A_ ...
,dir_local_mat ...
,n_rank ...
);
A_n_ = [];
end;%if flag_calculate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for niteration=1:n_iteration;

if (verbose); disp(sprintf(' %% [finished test_hnbtZRgumb_multi_7]')); end;
