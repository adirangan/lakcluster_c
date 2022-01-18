clear;

platform = 'rusty';%platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

verbose=1;
flag_recalc = 0;
flag_replot = 0;
flag_center = 0;
tolerance_master = 1e-2;
nf=0;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root);
str_code = sprintf('%s/halfloop_dev',dir_code);
dir_weilin = sprintf('/%s/rangan/dir_bcc/dir_weilin',string_root);
dir_mat = sprintf('%s/dir_mat',dir_weilin); if (~exist(dir_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_mat)); mkdir(dir_mat); end;
dir_jpg = sprintf('%s/dir_jpg',dir_weilin); if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
dir_local = sprintf('%s/dir_local',dir_weilin); if (~exist(dir_local,'dir')); disp(sprintf(' %% mkdir %s',dir_local)); mkdir(dir_local); end;
fname_count_sparse = sprintf('%s/count_sparse_.mat',dir_mat);
if (~exist('B_','var')); load(fname_count_sparse); end;
fname_count_logn_sparse = sprintf('%s/count_logn_sparse_.mat',dir_mat);
if (~exist('B_logn_','var')); load(fname_count_logn_sparse); end;

n_smp = size(B_,1);
n_umi = size(B_,2);
B_sum_1_ = sum(B_,1);
B_sum_2_ = sum(B_,2);

fname_count_pasca = sprintf('%s/count_pasca_.csv',dir_weilin);
fp = fopen(fname_count_pasca);
label_pasca_ = textscan(fp,'%s','Delimiter','\n','headerlines',1);
fclose(fp);
label_pasca_ = label_pasca_{1};
label_pasca_ = label_pasca_(1:n_smp);
for nl=0:numel(label_pasca_)-1; label_pasca_{1+nl}(strfind(label_pasca_{1+nl},'"'))=''; end;
for nl=0:numel(label_pasca_)-1; label_pasca_{1+nl}(strfind(label_pasca_{1+nl},'_'))=''; end;
[ ...
 label_pasca_enum_ ...
,n_u_label_pasca ...
,u_label_pasca_ ...
,index_nu_pasca_from_nall_ ...
,n_u_label_pasca_ ...
,index_nall_from_nu_pasca__ ...
] = ...
label_str_to_enum_1( ...
 label_pasca_ ...
);

fname_count_stage = sprintf('%s/count_stage_.csv',dir_weilin);
fp = fopen(fname_count_stage);
label_stage_ = textscan(fp,'%s','Delimiter','\n','headerlines',1);
fclose(fp);
label_stage_ = label_stage_{1};
label_stage_ = label_stage_(1:n_smp);
for nl=0:numel(label_stage_)-1; label_stage_{1+nl}(strfind(label_stage_{1+nl},'"'))=''; end;
for nl=0:numel(label_stage_)-1; label_stage_{1+nl}(strfind(label_stage_{1+nl},'_'))=''; end;
[ ...
 label_stage_enum_ ...
,n_u_label_stage ...
,u_label_stage_ ...
,index_nu_stage_from_nall_ ...
,n_u_label_stage_ ...
,index_nall_from_nu_stage__ ...
] = ...
label_str_to_enum_1( ...
 label_stage_ ...
);

fname_count_isoge = sprintf('%s/count_Isogenic_type_.csv',dir_weilin);
fp = fopen(fname_count_isoge);
label_isoge_ = textscan(fp,'%s','Delimiter','\n','headerlines',1);
fclose(fp);
label_isoge_ = label_isoge_{1};
label_isoge_ = label_isoge_(1:n_smp);
for nl=0:numel(label_isoge_)-1; label_isoge_{1+nl}(strfind(label_isoge_{1+nl},'"'))=''; end;
for nl=0:numel(label_isoge_)-1; label_isoge_{1+nl}(strfind(label_isoge_{1+nl},'_'))=''; end;
[ ...
 label_isoge_enum_ ...
,n_u_label_isoge ...
,u_label_isoge_ ...
,index_nu_isoge_from_nall_ ...
,n_u_label_isoge_ ...
,index_nall_from_nu_isoge__ ...
] = ...
label_str_to_enum_1( ...
 label_isoge_ ...
);

label_pasca_isoge_stage_ = cell(n_smp,1);
for nsmp=0:n_smp-1;
label_pasca_isoge_stage_{1+nsmp} = strcat(label_pasca_{1+nsmp},'_',label_isoge_{1+nsmp},'_',label_stage_{1+nsmp});
end;%for nsmp=0:n_smp-1;
[ ...
 label_pasca_isoge_stage_enum_ ...
,n_u_label_pasca_isoge_stage ...
,u_label_pasca_isoge_stage_ ...
,index_nu_pasca_isoge_stage_from_nall_ ...
,n_u_label_pasca_isoge_stage_ ...
,index_nall_from_nu_pasca_isoge_stage__ ...
] = ...
label_str_to_enum_1( ...
 label_pasca_isoge_stage_ ...
);

%%%%%%%%;
% set index 0 to represent the 'empty' restriction (i.e., full set). ;
%%%%%%%%;
fname_mat = sprintf('%s/index_pis___.mat',dir_mat);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
index_pis___ = cell(1+n_u_label_pasca,1+n_u_label_isoge,1+n_u_label_stage);
prefix_pis___ = cell(1+n_u_label_pasca,1+n_u_label_isoge,1+n_u_label_stage);
numel_pis___ = zeros(1+n_u_label_pasca,1+n_u_label_isoge,1+n_u_label_stage);
for nu_label_pasca=0:1+n_u_label_pasca-1;
tmp_index_pasca_ = transpose(0:n_smp-1); tmp_prefix_pasca = 'all';
if (nu_label_pasca> 0);
tmp_u_label_pasca = u_label_pasca_{1+nu_label_pasca-1};
tmp_prefix_pasca = tmp_u_label_pasca;
tmp_index_pasca_ = efind(cellfun(@(x)strcmp(x,tmp_u_label_pasca),label_pasca_));
end;%if (nu_label_pasca> 0);
for nu_label_isoge=0:1+n_u_label_isoge-1;
tmp_index_isoge_ = transpose(0:n_smp-1); tmp_prefix_isoge = 'all';
if (nu_label_isoge> 0);
tmp_u_label_isoge = u_label_isoge_{1+nu_label_isoge-1};
tmp_prefix_isoge = tmp_u_label_isoge;
tmp_index_isoge_ = efind(cellfun(@(x)strcmp(x,tmp_u_label_isoge),label_isoge_));
end;%if (nu_label_isoge> 0);
for nu_label_stage=0:1+n_u_label_stage-1;
tmp_index_stage_ = transpose(0:n_smp-1); tmp_prefix_stage = 'all';
if (nu_label_stage> 0);
tmp_u_label_stage = u_label_stage_{1+nu_label_stage-1};
tmp_prefix_stage = tmp_u_label_stage;
tmp_index_stage_ = efind(cellfun(@(x)strcmp(x,tmp_u_label_stage),label_stage_));
end;%if (nu_label_stage> 0);
tmp_prefix = sprintf('%s_%s_%s',tmp_prefix_pasca,tmp_prefix_isoge,tmp_prefix_stage);
tmp_index_ = intersectall({tmp_index_pasca_,tmp_index_isoge_,tmp_index_stage_});
if (verbose); disp(sprintf(' %% (%d,%d,%d) %s <-- %d',nu_label_pasca,nu_label_isoge,nu_label_stage,tmp_prefix,numel(tmp_index_))); end;
index_pis___{1+nu_label_pasca,1+nu_label_isoge,1+nu_label_stage} = tmp_index_;
prefix_pis___{1+nu_label_pasca,1+nu_label_isoge,1+nu_label_stage} = tmp_prefix;
numel_pis___(1+nu_label_pasca,1+nu_label_isoge,1+nu_label_stage) = numel(tmp_index_);
end;%for nu_label_stage=0:1+n_u_label_stage-1;
end;%for nu_label_isoge=0:1+n_u_label_isoge-1;
end;%for nu_label_pasca=0:1+n_u_label_pasca-1;
save(fname_mat ...
     ,'label_pasca_','n_u_label_pasca','n_u_label_pasca_','u_label_pasca_','label_pasca_enum_' ...
     ,'label_isoge_','n_u_label_isoge','n_u_label_isoge_','u_label_isoge_','label_isoge_enum_' ...
     ,'label_stage_','n_u_label_stage','n_u_label_stage_','u_label_stage_','label_stage_enum_' ...
     ,'index_pis___' ...
     ,'prefix_pis___' ...
     ,'numel_pis___' ...
     );
end;%if (~exist(fname_mat,file'));
load(fname_mat);

%%%%%%%%;
% global variables for post-processing. ;
%%%%%%%%;
define_global_test_halfloop_driver_6;

n_svd = 36;
date_diff_threshold = 0.5;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;
%%%%%%%%;
% Now step through each set of restrictions. ;
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nu_label_pasca=0:1+n_u_label_pasca-1;
for nu_label_isoge=0:1+n_u_label_isoge-1;
for nu_label_stage=0:1+n_u_label_stage-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% test with: ;
% nu_label_pasca = 0; nu_label_isoge = 1; nu_label_stage = 1;
%%%%%%%%;
flag_at_least_one_nu_zero = (nu_label_pasca==0) | (nu_label_isoge==0) | (nu_label_stage==0) ;
if flag_at_least_one_nu_zero;
try;
%%%%;
tmp_index_pasca_ = transpose(0:n_smp-1); tmp_prefix_pasca = 'all';
if (nu_label_pasca> 0);
tmp_u_label_pasca = u_label_pasca_{1+nu_label_pasca-1};
tmp_prefix_pasca = tmp_u_label_pasca;
tmp_index_pasca_ = efind(cellfun(@(x)strcmp(x,tmp_u_label_pasca),label_pasca_));
end;%if (nu_label_pasca> 0);
%%%%;
tmp_index_isoge_ = transpose(0:n_smp-1); tmp_prefix_isoge = 'all';
if (nu_label_isoge> 0);
tmp_u_label_isoge = u_label_isoge_{1+nu_label_isoge-1};
tmp_prefix_isoge = tmp_u_label_isoge;
tmp_index_isoge_ = efind(cellfun(@(x)strcmp(x,tmp_u_label_isoge),label_isoge_));
end;%if (nu_label_isoge> 0);
%%%%;
tmp_index_stage_ = transpose(0:n_smp-1); tmp_prefix_stage = 'all';
if (nu_label_stage> 0);
tmp_u_label_stage = u_label_stage_{1+nu_label_stage-1};
tmp_prefix_stage = tmp_u_label_stage;
tmp_index_stage_ = efind(cellfun(@(x)strcmp(x,tmp_u_label_stage),label_stage_));
end;%if (nu_label_stage> 0);
%%%%;
tmp_prefix = sprintf('%s_%s_%s',tmp_prefix_pasca,tmp_prefix_isoge,tmp_prefix_stage);
tmp_index_ = intersectall({tmp_index_pasca_,tmp_index_isoge_,tmp_index_stage_});
tmp_numel = numel(tmp_index_);
%%%%;
assert(strcmp(tmp_prefix,prefix_pis___{1+nu_label_pasca,1+nu_label_isoge,1+nu_label_stage}));
%%%%;
local_prefix = prefix_pis___{1+nu_label_pasca,1+nu_label_isoge,1+nu_label_stage};
local_index_ = index_pis___{1+nu_label_pasca,1+nu_label_isoge,1+nu_label_stage};
local_n_smp = numel_pis___(1+nu_label_pasca,1+nu_label_isoge,1+nu_label_stage);
assert(local_n_smp==numel(local_index_));
%%%%%%%%;
flag_local_n_smp_lt_24904 = (local_n_smp< 24904);
if flag_local_n_smp_lt_24904;
%%%%%%%%;
dir_local_mat = sprintf('%s/dir_%s_mat',dir_local,local_prefix);
if (~exist(dir_local_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_local_mat)); mkdir(dir_local_mat); end;
dir_local_jpg = sprintf('%s/dir_%s_jpg',dir_local,local_prefix);
if (~exist(dir_local_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_local_jpg)); mkdir(dir_local_jpg); end;
%%%%%%%%;
local_B_ = B_(1+local_index_,:);
local_nnz_umi_ = sum(local_B_~=0,1);
local_index_umi_retain_ = efind(local_nnz_umi_>=1); %<-- at least 1 nonzero entries. ;
local_n_umi = numel(local_index_umi_retain_);
local_B_logn_sg__ = B_logn_(1+local_index_,1+local_index_umi_retain_);
%%%%%%%%;
local_fname_pre = sprintf('%s/S_%s_B_logn_sg__',dir_local_mat,local_prefix);
[local_flag_skip,local_fname_mat] = open_fname_tmp(local_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~local_flag_skip;
[U_local_B_logn_sg__,S_local_B_logn_sg__,V_local_B_logn_sg__] = svds(local_B_logn_sg__,n_svd); S_local_B_logn_sg_ = diag(S_local_B_logn_sg__);
save(local_fname_mat ...
     ,'n_svd' ...
     ,'local_prefix','local_index_' ...
     ,'local_nnz_umi_','local_index_umi_retain_' ...
     ,'local_n_smp','local_n_umi' ...
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
label_xxxx_ = label_pasca_isoge_stage_(1+local_index_);
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
flag_calculate = strcmp(platform,'access1');
parameter = struct('type','parameter');
parameter.verbose = verbose;
parameter.date_diff_threshold = date_diff_threshold;
parameter.flag_force_create_mat = flag_force_create_mat;
parameter.flag_force_create_tmp = flag_force_create_tmp;
parameter.str_code = str_code;
parameter.halfloop_recursion_limit = 32;
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
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now post-process local clustering results. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
label_xxxx_pasca_ = label_pasca_(1+local_index_);
[ ...
 label_xxxx_pasca_enum_ ...
,n_u_label_xxxx_pasca ...
,u_label_xxxx_pasca_ ...
,index_nu_xxxx_pasca_from_nall_ ...
,n_u_label_xxxx_pasca_ ...
,index_nall_from_nu_xxxx_pasca__ ...
] = ...
label_str_to_enum_1( ...
 label_xxxx_pasca_ ...
);
%%%%%%%%;
label_xxxx_isoge_ = label_isoge_(1+local_index_);
[ ...
 label_xxxx_isoge_enum_ ...
,n_u_label_xxxx_isoge ...
,u_label_xxxx_isoge_ ...
,index_nu_xxxx_isoge_from_nall_ ...
,n_u_label_xxxx_isoge_ ...
,index_nall_from_nu_xxxx_isoge__ ...
] = ...
label_str_to_enum_1( ...
 label_xxxx_isoge_ ...
);
%%%%%%%%;
label_xxxx_stage_ = label_stage_(1+local_index_);
[ ...
 label_xxxx_stage_enum_ ...
,n_u_label_xxxx_stage ...
,u_label_xxxx_stage_ ...
,index_nu_xxxx_stage_from_nall_ ...
,n_u_label_xxxx_stage_ ...
,index_nall_from_nu_xxxx_stage__ ...
] = ...
label_str_to_enum_1( ...
 label_xxxx_stage_ ...
);
%%%%%%%;
local_flag_found_ = zeros(global_n_method,1);
local_n_cluster_found_ = zeros(global_n_method,1);
local_label_B__ = cell(global_n_method,1);
local_lpv_xxxx_ = zeros(global_n_method,1);
local_lpv_xxxx_pasca_ = zeros(global_n_method,1);
local_lpv_xxxx_isoge_ = zeros(global_n_method,1);
local_lpv_xxxx_stage_ = zeros(global_n_method,1);
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
local_lpv_xxxx_pasca = 0;
local_lpv_xxxx_isoge = 0;
local_lpv_xxxx_stage = 0;
local_label_B_ = {};
if (verbose>0); disp(sprintf(' %% %s not found, skipping',local_fname_mat)); end;
end;%if (~exist(local_fname_mat,'file'));
if ( exist(local_fname_mat,'file'));
tmp_ = load(local_fname_mat,'lpv_','lP0_','fla_','label_B__');
local_flag_found = 1;
local_n_cluster_found = numel(unique(tmp_.label_B__{1}));
local_label_B_ = tmp_.label_B__{end};
local_lpv_xxxx = label_to_label_enrichment_quad_4(label_xxxx_enum_,local_label_B_);
local_lpv_xxxx_pasca = label_to_label_enrichment_quad_4(label_xxxx_pasca_enum_,local_label_B_);
local_lpv_xxxx_isoge = label_to_label_enrichment_quad_4(label_xxxx_isoge_enum_,local_label_B_);
local_lpv_xxxx_stage = label_to_label_enrichment_quad_4(label_xxxx_stage_enum_,local_label_B_);
end;%if ( exist(fname_mat,'file'));
local_flag_found_(1+nmethod) = local_flag_found;
local_n_cluster_found_(1+nmethod) = local_n_cluster_found;
local_lpv_xxxx_(1+nmethod) = local_lpv_xxxx;
local_lpv_xxxx_pasca_(1+nmethod) = local_lpv_xxxx_pasca;
local_lpv_xxxx_isoge_(1+nmethod) = local_lpv_xxxx_isoge;
local_lpv_xxxx_stage_(1+nmethod) = local_lpv_xxxx_stage;
local_label_B__{1+nmethod} = local_label_B_;
end;%for nmethod=0:global_n_method-1;
%%%%%%%%;
if (verbose);
disp(sprintf(' %% -------------------------------------------------------------------------------------------------------------------------------- '));
for nl=0:numel(global_index_method_sub_)-1;
nmethod = global_index_method_sub_(1+nl);
disp(sprintf(' %% %21s (%.3d): \t all %+12.2f \t pasca %+12.2f \t isoge %+12.2f \t stage %+12.2f',global_prefix_method_{1+nmethod},local_n_cluster_found_(1+nmethod),local_lpv_xxxx_(1+nmethod),local_lpv_xxxx_pasca_(1+nmethod),local_lpv_xxxx_isoge_(1+nmethod),local_lpv_xxxx_stage_(1+nmethod)));
end;%for nl=0:numel(global_index_method_sub_)-1;
end;%if (verbose);
for str_method_ = {'hnbr0tZRgumb','hnbrtZRgumb','hnbtZRgumb'};
str_method = str_method_{1};
%%%%;
if strcmp(str_method,'hnbr0tZRgumb');
infix_hnbr0tZRgumb_0 = sprintf('hnbr0tZRgumb_r%d_r1',local_rank_estimate_B_logn_sg); 
infix_hnbr0tZRgumb_1 = sprintf('%s',infix_hnbr0tZRgumb_0);
str_infix = infix_hnbr0tZRgumb_1;
dir_local_mat_hnbr0tZRgumb = sprintf('%s/dir_tmp_%s/dir_tmp_%s',dir_local_mat,infix_hnbr0tZRgumb_0,infix_hnbr0tZRgumb_1);
str_sgtitle = dir_local_mat_hnbr0tZRgumb;
local_fname_output_label = sprintf('%s/output_label__.txt',dir_local_mat_hnbr0tZRgumb);
local_fname_nlpbra_label = sprintf('%s/nlpbra_label__.txt',dir_local_mat_hnbr0tZRgumb);
local_fname_nlpnex_label = sprintf('%s/nlpnex_label__.txt',dir_local_mat_hnbr0tZRgumb);
local_fname_louvain_mat = sprintf('%s/%s_r%d.mat',dir_local_mat,'louvain00pr_default',local_rank_estimate_B_logn_sg);
local_fname_umap_mat = sprintf('%s/%s_r%d.mat',dir_local_mat,'umap00pr_default',local_rank_estimate_B_logn_sg);
end;%if strcmp(str_method,'hnbr0tZRgumb');
%%%%;
if strcmp(str_method,'hnbrtZRgumb');
infix_hnbrtZRgumb_0 = sprintf('hnbrtZRgumb_r%d_r1',local_rank_estimate_B_logn_sg); 
infix_hnbrtZRgumb_1 = sprintf('%s',infix_hnbrtZRgumb_0);
str_infix = infix_hnbrtZRgumb_1;
dir_local_mat_hnbrtZRgumb = sprintf('%s/dir_tmp_%s/dir_tmp_%s',dir_local_mat,infix_hnbrtZRgumb_0,infix_hnbrtZRgumb_1);
str_sgtitle = dir_local_mat_hnbrtZRgumb;
local_fname_output_label = sprintf('%s/output_label__.txt',dir_local_mat_hnbrtZRgumb);
local_fname_nlpbra_label = sprintf('%s/nlpbra_label__.txt',dir_local_mat_hnbrtZRgumb);
local_fname_nlpnex_label = sprintf('%s/nlpnex_label__.txt',dir_local_mat_hnbrtZRgumb);
local_fname_louvain_mat = sprintf('%s/%s_r%d.mat',dir_local_mat,'louvain00pr_default',local_rank_estimate_B_logn_sg);
local_fname_umap_mat = sprintf('%s/%s_r%d.mat',dir_local_mat,'umap00pr_default',local_rank_estimate_B_logn_sg);
end;%if strcmp(str_method,'hnbrtZRgumb');
%%%%;
if strcmp(str_method,'hnbtZRgumb');
infix_hnbtZRgumb_0 = sprintf('hnbtZRgumb_r1'); 
infix_hnbtZRgumb_1 = sprintf('%s_g010',infix_hnbtZRgumb_0);
str_infix = infix_hnbtZRgumb_1;
dir_local_mat_hnbtZRgumb = sprintf('%s/dir_tmp_%s/dir_tmp_%s',dir_local_mat,infix_hnbtZRgumb_0,infix_hnbtZRgumb_1);
str_sgtitle = dir_local_mat_hnbtZRgumb;
local_fname_output_label = sprintf('%s/output_label__.txt',dir_local_mat_hnbtZRgumb);
local_fname_nlpbra_label = sprintf('%s/nlpbra_label__.txt',dir_local_mat_hnbtZRgumb);
local_fname_nlpnex_label = sprintf('%s/nlpnex_label__.txt',dir_local_mat_hnbtZRgumb);
local_fname_louvain_mat = sprintf('%s/%s.mat',dir_local_mat,'louvain00_default');
local_fname_umap_mat = sprintf('%s/%s.mat',dir_local_mat,'umap00_default');
end;%if strcmp(str_method,'hnbtZRgumb');
%%%%;
if numel(unique(label_xxxx_pasca_enum_))> 1;
parameter = struct('type','parameter');
parameter.flag_force_create_mat = 0;
parameter.flag_force_create_tmp = 0;
parameter.flag_replot = flag_replot;
label_tree_compilation_plot_0( ...
 parameter ...
,dir_local_mat ...
,sprintf('%s_pasca',str_infix) ...
,sprintf('%s (pasca)',str_sgtitle) ...
,local_fname_output_label ...
,local_fname_nlpbra_label ...
,local_fname_nlpnex_label ...
,local_B_logn_sg__ ...
,local_rank_estimate_B_logn_sg ...
,label_xxxx_pasca_enum_ ...
,local_fname_louvain_mat ...
,local_fname_umap_mat ...
);
end;%if numel(unique(label_xxxx_pasca_enum_))> 1;
%%%%;
if numel(unique(label_xxxx_isoge_enum_))> 1;
parameter = struct('type','parameter');
parameter.flag_force_create_mat = 0;
parameter.flag_force_create_tmp = 0;
parameter.flag_replot = flag_replot;
label_tree_compilation_plot_0( ...
 parameter ...
,dir_local_mat ...
,sprintf('%s_isoge',str_infix) ...
,sprintf('%s (isoge)',str_sgtitle) ...
,local_fname_output_label ...
,local_fname_nlpbra_label ...
,local_fname_nlpnex_label ...
,local_B_logn_sg__ ...
,local_rank_estimate_B_logn_sg ...
,label_xxxx_isoge_enum_ ...
,local_fname_louvain_mat ...
,local_fname_umap_mat ...
);
end;%if numel(unique(label_xxxx_isoge_enum_))> 1;
%%%%;
if numel(unique(label_xxxx_stage_enum_))> 1;
parameter = struct('type','parameter');
parameter.flag_force_create_mat = 0;
parameter.flag_force_create_tmp = 0;
parameter.flag_replot = flag_replot;
label_tree_compilation_plot_0( ...
 parameter ...
,dir_local_mat ...
,sprintf('%s_stage',str_infix) ...
,sprintf('%s (stage)',str_sgtitle) ...
,local_fname_output_label ...
,local_fname_nlpbra_label ...
,local_fname_nlpnex_label ...
,local_B_logn_sg__ ...
,local_rank_estimate_B_logn_sg ...
,label_xxxx_stage_enum_ ...
,local_fname_louvain_mat ...
,local_fname_umap_mat ...
);
end;%if numel(unique(label_xxxx_stage_enum_))> 1;
%%%%;
end;%for str_method_ = {'hnbr0tZRgumb','hnbrtZRgumb','hnbtZRgumb'};
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Finished post-processing. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (  exist(local_fname_mat,'file') &  exist(local_fname_shuffle_mat,'file') );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_local_n_smp_lt_24904;
catch; disp(sprintf(' %% Warning, could not run nu_label_pasca %d nu_label_isoge %d nu_label_stage %d ',nu_label_pasca,nu_label_isoge,nu_label_stage)); end;%try;
end;%if flag_at_least_one_nu_zero;
end;%for nu_label_stage=0:1+n_u_label_stage-1;
end;%for nu_label_isoge=0:1+n_u_label_isoge-1;
end;%for nu_label_pasca=0:1+n_u_label_pasca-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

disp('returning');return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;



