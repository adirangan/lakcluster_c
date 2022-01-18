function ...
[] = ...
test_loader_cluster_wrap_8( ...
 str_code ...
,dir_trunk ...
,str_label_A_ ...
,E_ ...
,rank_estimate_E ...
,str_X ...
,str_infix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);
% single clustering, no covariate correction. ;

%%%%%%%%;
na=0;
if (nargin<1+na); dir_trunk=pwd; end; na=na+1;
if (nargin<1+na); str_label_A_=[]; end; na=na+1;
if (nargin<1+na); E_=[]; end; na=na+1;
if (nargin<1+na); rank_estimate_E=[]; end; na=na+1;
if (nargin<1+na); str_X=[]; end; na=na+1;
if (nargin<1+na); str_infix=[]; end; na=na+1;
if (nargin<1+na); date_diff_threshold=[]; end; na=na+1;
if (nargin<1+na); flag_force_create_mat=[]; end; na=na+1;
if (nargin<1+na); flag_force_create_tmp=[]; end; na=na+1;
%%%%;
if ( isempty(rank_estimate_E)); rank_estimate_E = []; end;
if ( isempty(date_diff_threshold)); date_diff_threshold = 1.5; end;
if ( isempty(flag_force_create_mat)); flag_force_create_mat = 0; end;
if ( isempty(flag_force_create_tmp)); flag_force_create_tmp = 0; end;
%%%%%%%%;

dir_base = sprintf('%s/dir_loader_cluster',dir_trunk);
if (~exist(dir_base,'dir')); disp(sprintf(' %% mkdir %s',dir_base)); mkdir(dir_base); end;
str_prefix = sprintf('%s_%s',str_X,str_infix);
dir_cluster = sprintf('%s/dir_%s_cluster',dir_base,str_prefix);
if (~exist(dir_cluster,'dir')); disp(sprintf(' %% mkdir %s',dir_cluster)); mkdir(dir_cluster); end;
test_loader_cluster_5( ...
 str_code ...
,dir_cluster ...
,str_label_A_ ...
,E_ ...
,rank_estimate_E ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);