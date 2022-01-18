function ...
[ ...
 parameter ...
] = ...
test_hnbtZRgumb_multi_5( ...
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
% see test_hnbtZRgumb_multi_5_collect.m ;

if (nargin<1);
%%%%%%%%;
platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
for nl=3;%for nl=0:5;
na=0;   
if (nl==na); M =  178; N =  2e3; n_cluster = 1; n_rank =  2; n_iteration =  8; n_shuffle = 256; end; na=na+1;
if (nl==na); M =  178; N =  2e3; n_cluster = 3; n_rank =  4; n_iteration =  8; n_shuffle = 256; end; na=na+1;
if (nl==na); M =  563; N = 6325; n_cluster = 1; n_rank =  3; n_iteration =  4; n_shuffle =  64; end; na=na+1;
if (nl==na); M =  563; N = 6325; n_cluster = 6; n_rank =  9; n_iteration = 32; n_shuffle =  64; end; na=na+1;
if (nl==na); M = 1781; N =  2e4; n_cluster = 1; n_rank =  2; n_iteration =  4; n_shuffle =  64; end; na=na+1;
if (nl==na); M = 1781; N =  2e4; n_cluster = 1; n_rank =  4; n_iteration =  4; n_shuffle =  64; end; na=na+1;
if (nl==na); M = 1781; N =  2e4; n_cluster = 8; n_rank =  9; n_iteration =  4; n_shuffle =  64; end; na=na+1;
if (nl==na); M = 1781; N =  2e4; n_cluster = 8; n_rank = 12; n_iteration =  4; n_shuffle =  64; end; na=na+1;
flag_rerun = 0;
n_step = 21;
snr_ = linspace(0.5,1.0,n_step);
for nstep=1:n_step;
snr = snr_(nstep);
parameter = struct('type','parameter');
parameter = ...
test_hnbtZRgumb_multi_5( ...
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
end;%for nstep=1:n_step;
end;%for nl=0:5;
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
dir_base = sprintf('%s/dir_hnbtZRgumb_multi_5',dir_trunk);
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

if (verbose); disp(sprintf(' %% [entering test_hnbtZRgumb_multi_5]')); end;

if (~exist(dir_base,'dir')); disp(sprintf(' %% mkdir %s',dir_base)); mkdir(dir_base); end;
%%%%%%%%;
for niteration=1:n_iteration;
disp(sprintf(' %% niteration %d/%d',niteration,n_iteration));
%%%%%%%%;
rng(niteration);
[A_n_,label_A_,n_label_A_,pf_,pi_,snr_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);
[tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank);
UVS_n_ = tmp_U_*tmp_S_*transpose(tmp_V_);
US_n_ = tmp_U_*tmp_S_;
clear tmp_U_ tmp_S_ tmp_V_;
str_px = test_hnbtZRgumb_multi_xfix_5('px',M,N,snr,n_cluster,n_rank,niteration);
fname_px_mat = sprintf('%s/%s.mat',dir_base,str_px);
if (~exist(fname_px_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_px_mat));
save(fname_px_mat,'n_rank','M','N','snr','label_A_','n_label_A_','pf_','pi_','n_cluster');
end;%if (~exist(fname_px_mat,'file'));
if ( exist(fname_px_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_px_mat));
load(fname_px_mat);
end;%if ( exist(fname_px_mat,'file'));

str_dir = test_hnbtZRgumb_multi_xfix_5('dir',M,N,snr,n_cluster,n_rank,niteration);
dir_local_mat = sprintf('%s/%s',dir_base,str_dir);
if (~exist(dir_local_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_local_mat)); mkdir(dir_local_mat); end;

flag_calculate = strcmp(platform,'access1') | strcmp(platform,'eval1') | strcmp(platform,'OptiPlex');
parameter = struct('type','parameter');
parameter.verbose = verbose;
parameter.date_diff_threshold = date_diff_threshold;
parameter.flag_force_create_mat = flag_force_create_mat;
parameter.flag_force_create_tmp = flag_force_create_tmp;
parameter.str_code = str_code;
parameter.halfloop_recursion_limit = 32;
parameter.flag_skip_hnbtZRgumb = 0;
parameter.flag_skip_hnbrtZRgumb = 0;
parameter.n_shuffle = n_shuffle;
parameter.flag_orth_brute = 0;
if (M< 1024);
parameter.flag_orth_brute = 1; %<-- sample from full space of rotations. ;
end;%if (M< 1024);
if flag_calculate;
test_halfloop_driver_6( ...
 parameter ...
,A_n_ ...
,label_A_ ...
,dir_local_mat ...
,n_rank ...
);
end;%if flag_calculate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for niteration=1:n_iteration;

if (verbose); disp(sprintf(' %% [finished test_hnbtZRgumb_multi_5]')); end;
