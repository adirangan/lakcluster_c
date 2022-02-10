function ...
[ ...
 AZnV_ ...
,AnV_ ...
,ZnV_ ...
,V_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_DvX_from_ni_ver16( ...
 parameter ...
,ni ...
,xdrop_ ...
,trace__ ...
,mx__ ...
);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); ni=[]; end; na=na+1;
if (nargin<1+na); xdrop_=[]; end; na=na+1;
if (nargin<1+na); trace__=[]; end; na=na+1;
if (nargin<1+na); mx__=[]; end; na=na+1;

if isempty(ni); ni=0; end;
if isempty(xdrop_); xdrop_ = load_out_xdrop_from_str_ver0(parameter.str_out_xdrop_a_s0000); end;
if isempty(trace__); trace__ = load_trace__from_dir_ver0(parameter.dir_out_trace); end;
if isempty(mx__); mx__ = load_mx__from_parameter_ver0(parameter); end;

%%%%;
n_study = parameter.n_study;
flag_reverse = parameter.flag_reverse;
n_r_tmp = trace__.r_rem_s0000_(1+ni);
n_c_tmp = trace__.c_rem_s0000_(1+ni);

%%%%%%%%;
% pca. ;
%%%%%%%%;
parameter_pca = struct('type','parameter_pca');
parameter_pca.str_driver = 'pca_driver';
%parameter_pca.flag_force_create = 1;
parameter_pca.pca_mc_T = zeros(size(mx__.mc_T_)); parameter_pca.pca_mc_T(1)=1;
parameter_pca.pca_mc_A = zeros(size(mx__.mc_A_)); parameter_pca.pca_mc_A(xdrop_.ij_ckeep_(1:n_c_tmp))=1;
parameter_pca.pca_mr_A_ = cell(n_study,1);
parameter_pca.pca_mr_Z_ = cell(n_study,1);
for nstudy=0:n_study-1;
parameter_pca.pca_mr_A_{1+nstudy} = zeros(size(mx__.mr_A__{1+nstudy}));
parameter_pca.pca_mr_Z_{1+nstudy} = zeros(size(mx__.mr_Z__{1+nstudy}));
if (flag_reverse==0);
parameter_pca.pca_mr_A_{1+nstudy}(xdrop_.ij_rkeep_(1:n_r_tmp))=1;
parameter_pca.pca_mr_Z_{1+nstudy} = mx__.mr_Z__{1+nstudy};
end;%if (flag_reverse==0);
if (flag_reverse==1);
parameter_pca.pca_mr_Z_{1+nstudy}(xdrop_.ij_rkeep_(1:n_r_tmp))=1;
parameter_pca.pca_mr_A_{1+nstudy} = mx__.mr_A__{1+nstudy};
end;%if (flag_reverse==1);
end;%for nstudy=0:n_study-1;
parameter_pca.str_infix = sprintf('pca_DvX_ni%d',ni);
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
parameter_pca_proj.pca_mc_A = zeros(size(mx__.mc_A_)); parameter_pca_proj.pca_mc_A(xdrop_.ij_ckeep_(1:n_c_tmp))=1;
parameter_pca_proj.pca_mr_A_ = cell(n_study,1);
parameter_pca_proj.pca_mr_Z_ = cell(n_study,1);
for nstudy=0:n_study-1;
parameter_pca_proj.pca_mr_A_{1+nstudy} = mx__.mr_A__{1+nstudy};
parameter_pca_proj.pca_mr_Z_{1+nstudy} = mx__.mr_Z__{1+nstudy};
end;%for nstudy=0:n_study-1;
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
