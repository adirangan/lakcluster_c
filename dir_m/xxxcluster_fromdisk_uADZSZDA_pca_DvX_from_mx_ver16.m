function ...
[ ...
 parameter ...
,AZnV_ ...
,AnV_ ...
,ZnV_ ...
,V_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_DvX_from_mx_ver16( ...
 parameter ...
,pca_mr_A_ ...
,pca_mr_Z_ ...
,pca_mc_A ...
,pca_str_infix ...
,mx__ ...
);

str_thisfunction = 'xxxcluster_fromdisk_uADZSZDA_pca_DvX_from_mx_ver16';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); pca_mr_A_=[]; end; na=na+1;
if (nargin<1+na); pca_mr_Z_=[]; end; na=na+1;
if (nargin<1+na); pca_mc_A=[]; end; na=na+1;
if (nargin<1+na); mx__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter');
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;

if isempty(mx__); mx__ = load_mx__from_parameter_ver0(parameter); end;
%%%%;
if isempty(parameter,'dir_0in'); disp(sprintf(' %% Warning, parameter.dir_0in undefined in %s',str_thisfunction)); end;
dir_0in = parameter.dir_0in;
str_suffix = sprintf('%s','analyze');
dir_tmp = sprintf('%s_%s',dir_0in,str_suffix);
if ~exist(dir_tmp,'dir'); disp(sprintf(' %% mkdir %s',dir_tmp)); mkdir(dir_tmp); end;
%%%%;
if isempty(parameter,'str_name_s0000');
parameter.str_name_s0000 = 'default';
end;%if isempty(parameter,'str_name_s0000');
str_name_s0000 = parameter.str_name_s0000;
%%%%;
if isempty(parameter,'dir_out_s0000');
str_suffix = sprintf('%s','analyze');
dir_out_s0000 = sprintf('%s_%s/dir_%s',dir_0in,str_suffix,str_name_s0000);
if ~exist(dir_out_s0000,'dir'); disp(sprintf(' %% mkdir %s',dir_out_s0000)); mkdir(dir_out_s0000); end;
parameter.dir_out_s0000 = dir_out_s0000;
end;%if isempty(parameter,'dir_out_s0000');
%%%%;
str_name_s0000 = parameter.str_name_s0000;
dir_out_s0000 = parameter.dir_out_s0000;

%%%%;
n_study = parameter.n_study;
flag_reverse = parameter.flag_reverse;

%%%%%%%%;
% pca. ;
%%%%%%%%;
parameter_pca = struct('type','parameter_pca');
parameter_pca.str_driver = 'pca_driver';
%parameter_pca.flag_force_create = 1;
parameter_pca.pca_mc_T = zeros(size(mx__.mc_T_)); parameter_pca.pca_mc_T(1)=1;
parameter_pca.pca_mc_A = mx__.mc_A_;
if ~isempty(pca_mc_A); parameter_pca.pca_mc_A = pca_mc_A; end;
parameter_pca.pca_mr_A_ = cell(n_study,1);
parameter_pca.pca_mr_Z_ = cell(n_study,1);
for nstudy=0:n_study-1;
parameter_pca.pca_mr_A_{1+nstudy} = mx__.mr_A__{1+nstudy};
parameter_pca.pca_mr_Z_{1+nstudy} = mx__.mr_Z__{1+nstudy};
end;%for nstudy=0:n_study-1;
if ~isempty(pca_mr_A_); parameter_pca.pca_mr_A_ = pca_mr_A_; end;
if ~isempty(pca_mr_Z_); parameter_pca.pca_mr_Z_ = pca_mr_Z_; end;
flag_continue=1;
str_tmp = sprintf('%.16d',floor(1e16*rem(now,1)));
parameter_pca.str_infix = sprintf('pca_DvX_%d',ni);
[ ...
 parameter ...
 parameter_pca ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_ver16( ...
 parameter ...
,parameter_pca ...
);
V_ = mda_read_r8(parameter_pca.str_V);
%%%%%%%%;

%%%%%%%%;
% projection. ;
%%%%%%%%;
parameter_pca_proj = struct('type','parameter_pca_proj');
parameter_pca_proj.str_driver = 'pca_proj_driver';
parameter_pca_proj.str_V = parameter_pca.str_V;
%parameter_pca_proj.flag_force_create = 1;
parameter_pca_proj.pca_mc_T = zeros(size(mx__.mc_T_)); parameter_pca_proj.pca_mc_T(1)=1;
parameter_pca_proj.pca_mc_A = mx__.mc_A_;
if ~isempty(pca_mc_A); parameter_pca_proj.pca_mc_A = pca_mc_A; end;
parameter_pca_proj.pca_mr_A_ = cell(n_study,1);
parameter_pca_proj.pca_mr_Z_ = cell(n_study,1);
for nstudy=0:n_study-1;
parameter_pca_proj.pca_mr_A_{1+nstudy} = mx__.mr_A__{1+nstudy};
parameter_pca_proj.pca_mr_Z_{1+nstudy} = mx__.mr_Z__{1+nstudy};
end;%for nstudy=0:n_study-1;
%if ~isempty(pca_mr_A_); parameter_pca_proj.pca_mr_A_ = pca_mr_A_; end;
%if ~isempty(pca_mr_Z_); parameter_pca_proj.pca_mr_Z_ = pca_mr_Z_; end;
parameter_pca_proj.str_infix = sprintf('pca_proj_DvX_ni%d',ni);
[ ...
 parameter ...
 parameter_pca_proj ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_ver16( ...
 parameter ...
,parameter_pca_proj ...
);
%%%%%%%%;
AnV_ = mda_read_r8(parameter_pca_proj.str_AnV);
ZnV_ = mda_read_r8(parameter_pca_proj.str_ZnV);
AZnV_ = AnV_ + ZnV_;
