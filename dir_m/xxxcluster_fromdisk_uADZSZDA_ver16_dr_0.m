function ...
[ ...
 parameter ...
] = ...
  xxxcluster_fromdisk_uADZSZDA_ver16_dr_0( ...
 parameter ...
);

if nargin<1;
disp(sprintf(' %% Testing xxxcluster_fromdisk_uADZSZDA_ver16_dr_0. ;'));
disp(sprintf(' %% First running test_stripped_xxxcluster_fromdisk_uADZSZDA_ver16. ;'));
disp(sprintf(' %% Now running xxxcluster_fromdisk_uADZSZDA_ver16_dr_0 on same files. ;'));
dir_trunk = sprintf('%s/dir_test_xxxcluster_fromdisk_ver16',pwd);
test_stripped_xxxcluster_fromdisk_ver16(struct('flag_verbose',0));
parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.dir_trunk = dir_trunk;
parameter.study_name_of_branch_s_ = {'dir_study00'};
parameter.study_name_without_extension_s_ = {'study00'};
parameter.mds_name_plus_extension_s_ = {'study00_mds_tsv.txt'};
parameter.ent_cutoff = 0.0045022835561375; %<-- test_stripped_xxxcluster_fromdisk_uADZSZDA_ver16;
xxxcluster_fromdisk_uADZSZDA_ver16_dr_0(parameter);
disp(sprintf(' %% Now loading results from test_stripped_xxxcluster_fromdisk_uADZSZDA_ver16. ;'));

%%%%%%%%;
dir_out_0_ = sprintf('%s/dir_test_maf01_analyze/dir_test_maf01_dex_p10_D_m2r1_g050',dir_trunk);
fname_bimext_0 = sprintf('%s/dir_test_maf01/test_maf01_bim.ext',dir_trunk);
[ ...
 n_snp_bimext_0 ...
,bimext_khr_0_ ...
,bimext_vid_0_ ...
,bimext_gdi_0_ ...
,bimext_pdi_0_ ...
,bimext_al1_0_ ...
,bimext_al2_0_ ...
,bimext_alt_0_ ...
,bimext_ent_0_ ...
,bimext_frq_0_ ...
,bimext_mss_0_ ...
,bimext_maf_0_ ...
,bimext_name_0_ ...
,bimext_0_ ...
] = ...
load_bimext_ver1( ...
fname_bimext_0 ...
);
fname_famext_0 = sprintf('%s/dir_test_maf01/test_maf01_fam.ext',dir_trunk);
[ ...
 n_patient_famext_0 ...
,famext_fid_0_ ...
,famext_iid_0_ ...
,famext_yid_0_ ...
,famext_xid_0_ ...
,famext_sex_0_ ...
,famext_dvx_0_ ...
,famext_dir_0_ ...
,famext_fidandiid_0_ ...
,famext_0_ ...
  ] = ...
load_famext_ver1( ...
 fname_famext_0 ...
);
%%%%%%%%;
dir_out_1_ = sprintf('%s/dir_test_dr_maf01_analyze/dir_test_dr_maf01_dex_p10_D_m2r1_g050',dir_trunk);
fname_bimext_1 = sprintf('%s/dir_test_dr_maf01/test_dr_maf01_bim.ext',dir_trunk);
[ ...
 n_snp_bimext_1 ...
,bimext_khr_1_ ...
,bimext_vid_1_ ...
,bimext_gdi_1_ ...
,bimext_pdi_1_ ...
,bimext_al1_1_ ...
,bimext_al2_1_ ...
,bimext_alt_1_ ...
,bimext_ent_1_ ...
,bimext_frq_1_ ...
,bimext_mss_1_ ...
,bimext_maf_1_ ...
,bimext_name_1_ ...
,bimext_1_ ...
] = ...
load_bimext_ver1( ...
fname_bimext_1 ...
);
fname_famext_1 = sprintf('%s/dir_test_dr_maf01/test_dr_maf01_fam.ext',dir_trunk);
[ ...
 n_patient_famext_1 ...
,famext_fid_1_ ...
,famext_iid_1_ ...
,famext_yid_1_ ...
,famext_xid_1_ ...
,famext_sex_1_ ...
,famext_dvx_1_ ...
,famext_dir_1_ ...
,famext_fidandiid_1_ ...
,famext_1_ ...
  ] = ...
load_famext_ver1( ...
 fname_famext_1 ...
);
%%%%%%%%;
[bimext_name_cap_,ij_ns0_from_ncap_,ij_ns1_from_ncap_] = intersect(bimext_name_0_,bimext_name_1_,'stable');
disp(sprintf(' %% bimext_0 error %0.16f',fnorm(ij_ns0_from_ncap_ - transpose([1:n_snp_bimext_0]))));
disp(sprintf(' %% bimext_1 error %0.16f',fnorm(ij_ns1_from_ncap_ - transpose([1:n_snp_bimext_1]))));
%%%%%%%%;
[famext_fidandiid_cap_,ij_ns0_from_ncap_,ij_ns1_from_ncap_] = intersect(famext_fidandiid_0_,famext_fidandiid_1_,'stable');
disp(sprintf(' %% famext_0 error %0.16f',fnorm(ij_ns0_from_ncap_ - transpose([1:n_patient_famext_0]))));
disp(sprintf(' %% famext_1 error %0.16f',fnorm(ij_ns1_from_ncap_ - transpose([1:n_patient_famext_1]))));
%%%%%%%%;
fname_xdrop_0 = sprintf('%s/dir_test_maf01_analyze/dir_test_maf01_dex_p10_D_m2r1_g050/out_xdrop_a.txt',dir_trunk);
xdrop_0_ = load_out_xdrop_from_str_ver0(fname_xdrop_0);
fname_xdrop_1 = sprintf('%s/dir_test_dr_maf01_analyze/dir_test_dr_maf01_dex_p10_D_m2r1_g050/out_xdrop_a.txt',dir_trunk);
xdrop_1_ = load_out_xdrop_from_str_ver0(fname_xdrop_1);
disp(sprintf(' %% xdrop_.index_rdrop_ error %0.16f',fnorm(xdrop_0_.index_rdrop_ - xdrop_1_.index_rdrop_)));
disp(sprintf(' %% xdrop_.index_cdrop_ error %0.16f',fnorm(xdrop_0_.index_cdrop_ - xdrop_1_.index_cdrop_)));
disp('returning'); return;
end;%if nargin<1;

na=0;  
if (nargin<1+na); parameter=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
%%%%%%%%;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=1; end;
flag_verbose=parameter.flag_verbose;
flag_disp = flag_verbose;
if (flag_verbose); disp(sprintf(' %% [entering xxxcluster_fromdisk_ver16_dr_0]')); end;
%%%%%%%%;
if ~isfield(parameter,'flag_force_create'); parameter.flag_force_create=0; end;
flag_force_create=parameter.flag_force_create;
%%%%%%%%;
if ~isfield(parameter,'rseed'); parameter.rseed=0; end;
rseed=parameter.rseed;
%%%%%%%%;
if ~isfield(parameter,'str_lak_vs_dex'); parameter.str_lak_vs_dex='dex'; end;
str_lak_vs_dex=parameter.str_lak_vs_dex;
%%%%%%%%;
if ~isfield(parameter,'dir_trunk'); parameter.dir_trunk=pwd; end;
dir_trunk = parameter.dir_trunk;
%%%%%%%%;
if ~isfield(parameter,'dir_code'); parameter.dir_code=sprintf('%s/..',pwd); end;
dir_code = parameter.dir_code;
%%%%%%%%;
if ~isfield(parameter,'study_name_of_branch_s_'); parameter.study_name_of_branch_s_=[]; end;
study_name_of_branch_s_ = parameter.study_name_of_branch_s_;
%%%%%%%%;
if ~isfield(parameter,'study_name_without_extension_s_'); parameter.study_name_without_extension_s_=[]; end;
study_name_without_extension_s_ = parameter.study_name_without_extension_s_;
%%%%%%%%;
if ~isfield(parameter,'mds_name_plus_extension_s_'); parameter.mds_name_plus_extension_s_=[]; end;
mds_name_plus_extension_s_ = parameter.mds_name_plus_extension_s_;
%%%%%%%%;
if ~isfield(parameter,'mds_delimiter'); parameter.mds_delimiter='\t'; end;
mds_delimiter = parameter.mds_delimiter;
%%%%%%%%;
if ~isfield(parameter,'mds_headerlines'); parameter.mds_headerlines=0; end;
mds_headerlines = parameter.mds_headerlines;
%%%%%%%%;
if ~isfield(parameter,'str_output_prefix'); parameter.str_output_prefix='test_dr'; end;
str_output_prefix = parameter.str_output_prefix;
%%%%%%%%;
if ~isfield(parameter,'ent_cutoff'); parameter.ent_cutoff=0.001; end;
ent_cutoff = parameter.ent_cutoff;
%%%%%%%%;
if ~isfield(parameter,'maf_cutoff'); parameter.maf_cutoff=0.010; end;
maf_cutoff = parameter.maf_cutoff;
%%%%%%%%;
if ~isfield(parameter,'gamma'); parameter.gamma=0.05; end;
gamma = parameter.gamma;
%%%%%%%%;
if ~isfield(parameter,'n_mds_repl'); parameter.n_mds_repl=1; end;
n_mds_repl = parameter.n_mds_repl;
%%%%%%%%;
if ~isfield(parameter,'n_shuffle'); parameter.n_shuffle=32; end;
n_shuffle = parameter.n_shuffle;
%%%%%%%%;

setup_0;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% dir_trunk: %s;',dir_trunk)); end;

if ischar(study_name_without_extension_s_);
n_study = 1;
tmp_char = study_name_without_extension_s_;
study_name_without_extension_s_ = cell(n_study,1); for nstudy=0:n_study-1; study_name_without_extension_s_{1+nstudy} = tmp_char; end;
clear tmp_char;
end;%if ischar(study_name_without_extension_s_);
n_study = numel(study_name_without_extension_s_);
if (n_study<1); disp(sprintf(' %% Warning, no studies specified in xxxcluster_fromdisk_uADZSZDA_ver16_dr_0, returning')); return; end;
if isempty(study_name_of_branch_s_);
study_name_of_branch_s_ = cell(n_study,1); for nstudy=0:n_study-1; study_name_of_branch_s_{1+nstudy} = pwd; end;
end;%if isempty(study_name_of_branch_s_);
if ischar(study_name_of_branch_s_);
tmp_char = study_name_of_branch_s_;
study_name_of_branch_s_ = cell(n_study,1); for nstudy=0:n_study-1; study_name_of_branch_s_{1+nstudy} = tmp_char; end;
clear tmp_char;
end;%if ischar(study_name_of_branch_s_);
if ischar(mds_name_plus_extension_s_);
tmp_char = mds_name_plus_extension_s_;
mds_name_plus_extension_s_ = cell(n_study,1); for nstudy=0:n_study-1; mds_name_plus_extension_s_{1+nstudy} = tmp_char; end;
clear tmp_char;
end;%if ischar(mds_name_plus_extension_s_);

flag_found = 1;
n_patient_s_ = zeros(n_study,1);
n_snp_s_ = zeros(n_study,1);
for nstudy=0:n_study-1;
%%%%;
tmp_fname_bed = sprintf('%s/%s/%s.bed',dir_trunk,study_name_of_branch_s_{1+nstudy},study_name_without_extension_s_{1+nstudy});
if ~exist(tmp_fname_bed,'file'); disp(sprintf(' %% Warning, %s not found in xxxcluster_fromdisk_uADZSZDA_ver16_dr_0, returning',tmp_fname_bed)); return; end;
flag_found = flag_found*exist(tmp_fname_bed,'file');
if (flag_verbose & flag_found); disp(sprintf(' %% %s found',tmp_fname_bed)); end;
%%%%;
tmp_fname_bim = sprintf('%s/%s/%s.bim',dir_trunk,study_name_of_branch_s_{1+nstudy},study_name_without_extension_s_{1+nstudy});
if ~exist(tmp_fname_bim,'file'); disp(sprintf(' %% Warning, %s not found in xxxcluster_fromdisk_uADZSZDA_ver16_dr_0, returning',tmp_fname_bim)); return; end;
flag_found = flag_found*exist(tmp_fname_bim,'file');
if (flag_verbose & flag_found); disp(sprintf(' %% %s found',tmp_fname_bim)); end;
n_snp_s_(1+nstudy) = wc_0(tmp_fname_bim);
%%%%;
tmp_fname_fam = sprintf('%s/%s/%s.fam',dir_trunk,study_name_of_branch_s_{1+nstudy},study_name_without_extension_s_{1+nstudy});
if ~exist(tmp_fname_fam,'file'); disp(sprintf(' %% Warning, %s not found in xxxcluster_fromdisk_uADZSZDA_ver16_dr_0, returning',tmp_fname_fam)); return; end;
flag_found = flag_found*exist(tmp_fname_fam,'file');
if (flag_verbose & flag_found); disp(sprintf(' %% %s found',tmp_fname_fam)); end;
n_patient_s_(1+nstudy) = wc_0(tmp_fname_fam);
%%%%;
if (flag_verbose); disp(sprintf(' %% n_patient_s_(1+%d) = %d',nstudy,n_patient_s_(1+nstudy))); end;
end;%for nstudy=0:n_study-1;

n_patient_cup = sum(n_patient_s_);
if isempty(mds_name_plus_extension_s_);
if (flag_verbose); disp(sprintf(' %% no mds file given')); end;
mds_pv__ = zeros(n_patient_cup,0);
end;%if isempty(mds_name_plus_extension_s_);
%%%%%%%%;
if numel(mds_name_plus_extension_s_)==1;
flag_found=0;
%%%%;
tmp_fname_mds = sprintf('%s/%s/%s',dir_trunk,study_name_of_branch_s_{1},mds_name_plus_extension_s_{1});
if ~flag_found & exist(tmp_fname_mds,'file');
if (flag_verbose); disp(sprintf(' %% %s found',tmp_fname_mds)); end;
mds_pv__ = textread(tmp_fname_mds,'','delimiter',mds_delimiter,'headerlines',mds_headerlines);
flag_found = 1;
end;%if ~flag_found & exist(tmp_fname_mds,'file');
%%%%;
tmp_fname_mds = sprintf('%s',mds_name_plus_extension_s_{1});
if ~flag_found & exist(tmp_fname_mds,'file');
if (flag_verbose); disp(sprintf(' %% %s found',tmp_fname_mds)); end;
mds_pv__ = textread(tmp_fname_mds,'','delimiter',mds_delimiter,'headerlines',mds_headerlines);
flag_found = 1;
end;%if ~flag_found & exist(tmp_fname_mds,'file');
%%%%;
if (~flag_found); disp(sprintf(' %% Warning, %s not found in xxxcluster_fromdisk_uADZSZDA_ver16_dr_0, returning',tmp_fname_mds)); return; end;
end;%if numel(mds_name_plus_extension_s_)==1;
%%%%%%%%;
if numel(mds_name_plus_extension_s_)> 1;
for nstudy=0:n_study-1;
flag_found=0;
%%%%;
tmp_fname_mds = sprintf('%s/%s/%s',dir_trunk,study_name_of_branch_s_{1+nstudy},mds_name_plus_extension_s_{1+nstudy});
if ~flag_found & exist(tmp_fname_mds,'file');
if (flag_verbose); disp(sprintf(' %% %s found',tmp_fname_mds)); end;
tmp_mds_pv__ = textread(tmp_fname_mds,'','delimiter',mds_delimiter,'headerlines',mds_headerlines);
flag_found = 1;
end;%if ~flag_found & exist(tmp_fname_mds,'file');
%%%%;
tmp_fname_mds = sprintf('%s',mds_name_plus_extension_s_{1+nstudy});
if ~flag_found & exist(tmp_fname_mds,'file');
if (flag_verbose); disp(sprintf(' %% %s found',tmp_fname_mds)); end;
tmp_mds_pv__ = textread(tmp_fname_mds,'','delimiter',mds_delimiter,'headerlines',mds_headerlines);
flag_found = 1;
end;%if ~flag_found & exist(tmp_fname_mds,'file');
%%%%;
if (~flag_found); disp(sprintf(' %% Warning, %s not found in xxxcluster_fromdisk_uADZSZDA_ver16_dr_0, returning',tmp_fname_mds)); return; end;
mds_pv__ = [mds_pv__ ; tmp_mds_pv__];
end;%for nstudy=0:n_study-1;
end;%if numel(mds_name_plus_extension_s_)> 1;
%%%%%%%%;
[n_mds_p,n_mds_v] = size(mds_pv__);
if (flag_verbose); disp(sprintf(' %% mds_pv__ of size (%d,%d)',n_mds_p,n_mds_v)); end;
mds_fidandiid_p_ = [];
if (flag_verbose); disp(sprintf(' %% setting random seed to %d',rseed)); end;
rng(rseed); nf=0;

%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% We first convert the dosage-data to b16-data. ;')); end;
if (flag_verbose); disp(sprintf(' %% This encodes each dosage value as a 3-bit binary array. ;')); end;
if (flag_verbose); disp(sprintf(' %% dosage ''0'' is encoded as [ 1 0 0 ]. ;')); end;
if (flag_verbose); disp(sprintf(' %% dosage ''1'' is encoded as [ 0 1 0 ]. ;')); end;
if (flag_verbose); disp(sprintf(' %% dosage ''2'' is encoded as [ 0 0 1 ]. ;')); end;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% You can also set certain people to be excluded (similar to a famex file). ;')); end;
if (flag_verbose); disp(sprintf(' %% For now we skip this step. ;')); end;
n_famex = 0;
famex_fidandiid_p_ = cell(n_famex,1);
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% When converting the bed, bim and fam data into a b16 array, ;')); end;
if (flag_verbose); disp(sprintf(' %% we must specify a few parameters. ;')); end;
if (flag_verbose); disp(sprintf(' %% Some of the more important parameter are: ;')); end;
if (flag_verbose); disp(sprintf(' %% maf_cutoff: scalar double: this is the minor-allele-frequency (maf) cutoff for snps. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t Only snps with a strictly higher maf will be considered. ;')); end;
if (flag_verbose); disp(sprintf(' %% ent_cutoff: scalar double: this is the entropy-cutoff for snps. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t We typically expect the dosage-data for any particular snp ;')); end;
if (flag_verbose); disp(sprintf(' %%\t to have a frequency of 0, 1 and 2 values ;')); end;
if (flag_verbose); disp(sprintf(' %%\t (i.e., ''aa'', ''Aa'' and ''AA'' values) ;')); end;
if (flag_verbose); disp(sprintf(' %%\t which is close to pp, 2pq and qq, respectively, ;')); end;
if (flag_verbose); disp(sprintf(' %%\t where p is the associated minor-allele-frequency (maf). ;')); end;
if (flag_verbose); disp(sprintf(' %%\t The entropy (I_opt) measures the deviation from this expectation, ;')); end;
if (flag_verbose); disp(sprintf(' %%\t with higher values indicating a greater deviation. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t Only snps with a strictly lower I_opt will be considered.  ;')); end;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% This is also where you need to specify the mds-components you are correcting for, ;')); end;
if (flag_verbose); disp(sprintf(' %% As well as the patients to be excluded (i.e., famex). ;')); end;
parameter.flag_verbose = max(0,flag_verbose-1);
[ ...
 parameter ...
] = ...
bed_to_b16_flip_ver7( ...
 parameter ...
,n_study ...
,study_name_of_branch_s_ ...
,study_name_without_extension_s_ ...
,n_mds_p ...
,n_mds_v ...
,mds_pv__ ...
,mds_fidandiid_p_ ...
,n_famex ...
,famex_fidandiid_p_ ...
);
%%%%%%%%;
str_maf = sprintf('maf%.2d',floor(100*parameter.maf_cutoff));
str_output_prefix_plus_maf = sprintf('%s_%s',str_output_prefix,str_maf);
str_output_prefix_local = sprintf('%s_',str_output_prefix_plus_maf);
dir_out = sprintf('%s/dir_%s',parameter.dir_trunk,str_output_prefix_plus_maf);
if (flag_verbose>0); disp(sprintf(' %% str_output_prefix_local %s ;\n %% dir_out %s ;\n',str_output_prefix_local,dir_out)); end;
%%%%%%%%;

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we do not uncompress the b16 file. ;')); end;
if (flag_verbose); disp(sprintf(' %% This step is too slow for larger data-sets. ;')); end;
if (flag_verbose); disp(sprintf(' %% We can still load the fam.ext and bim.ext files. ;')); end;
fname_b16 = sprintf('%s/%s_A_full_n.b16',dir_out,str_output_prefix_plus_maf);
if (flag_verbose); disp(sprintf(' %% fname_b16: %s',fname_b16)); end;
if ~exist(fname_b16,'file'); disp(sprintf(' %% Warning, %s not found in xxxcluster_fromdisk_uADZSZDA_ver16_dr_0, returning',fname_b16)); return; end;
%A_full_n_ = binary_uncompress(fname_b16);
%%%%%%%%;
fname_famext = sprintf('%s/%s_fam.ext',dir_out,str_output_prefix_plus_maf);
if (flag_verbose); disp(sprintf(' %% fname_famext: %s',fname_famext)); end;
if ~exist(fname_famext,'file'); disp(sprintf(' %% Warning, %s not found in xxxcluster_fromdisk_uADZSZDA_ver16_dr_0, returning',fname_famext)); return; end;
[ ...
 n_patient_famext ...
,famext_fid_ ...
,famext_iid_ ...
,famext_yid_ ...
,famext_xid_ ...
,famext_sex_ ...
,famext_dvx_ ...
,famext_dir_ ...
,famext_fidandiid_ ...
,famext_ ...
  ] = ...
load_famext_ver1( ...
 fname_famext ...
);
mr_dvx_ = famext_dvx_;
%%%%%%%%;
fname_bimext = sprintf('%s/%s_bim.ext',dir_out,str_output_prefix_plus_maf);
if (flag_verbose); disp(sprintf(' %% fname_bimext: %s',fname_bimext)); end;
if ~exist(fname_bimext,'file'); disp(sprintf(' %% Warning, %s not found in xxxcluster_fromdisk_uADZSZDA_ver16_dr_0, returning',fname_bimext)); return; end;
[ ...
 n_snp_bimext ...
,bimext_khr_ ...
,bimext_vid_ ...
,bimext_gdi_ ...
,bimext_pdi_ ...
,bimext_al1_ ...
,bimext_al2_ ...
,bimext_alt_ ...
,bimext_ent_ ...
,bimext_frq_ ...
,bimext_mss_ ...
,bimext_maf_ ...
,bimext_name_ ...
,bimext_ ...
] = ...
load_bimext_ver1( ...
fname_bimext ...
);
[u_bimext_vid_,index_bim_ext_from_u_,index_bim_u_from_ext_] = unique(bimext_vid_,'stable');
index_bim_ext_from_u_ = index_bim_ext_from_u_ - 1;
index_bim_u_from_ext_ = index_bim_u_from_ext_ - 1;
%%%%;

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we run the basic biclustering method on the b16-data. ;')); end;
if (flag_verbose); disp(sprintf(' %% Some of the important parameters are: ;')); end;
if (flag_verbose); disp(sprintf(' %% n_mds: scalar integer: number of mds-components to correct for. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t Note that the mds-components will be binarized into orthants. ;')); end;
if (flag_verbose); disp(sprintf(' %% n_mds_repl: scalar integer: number of times to replicate the mds-components. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t This specifically refers to rotation+replication across orthants. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t Larger numbers are better, but 1 or 2 is usually sufficient. ;')); end;
if (flag_verbose); disp(sprintf(' %% ij_mds_use_: integer array of size n_mds: 1-based indices of mds-components to correct for. ;')); end;
if (flag_verbose); disp(sprintf(' %% gamma: scalar double: fraction of rows or columns eliminated at each iteration of the bicluster algorithm. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t Smaller numbers are more accurate, but larger numbers require less computation time. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t Empirically, gamma = 0.01 or 0.05 is usually sufficient. ;')); end;
if (flag_verbose); disp(sprintf(' %% flag_force_create: scalar integer: set to 1 to force recreation of all temporary files. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t if set to 0 will load temporary files from disk. ;')); end;
n_mds_0in = n_mds_v; ij_mds_use_ = [1:n_mds_0in];
str_prefix = str_output_prefix_plus_maf;
parameter.str_prefix = str_prefix;
parameter.n_mds = n_mds_0in;
parameter.ij_mds_use_ = ij_mds_use_;
[ ...
 parameter ...
] = ...
xxxcluster_fromdisk_uADZSZDA_ver16( ...
 parameter ...
);

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we load the results of this single biclustering run. ;')); end;
parameter.str_out_xdrop_a_s0000 = sprintf('%s/out_xdrop_a.txt',parameter.dir_out_s0000);
[ ...
 xdrop_ ...
] = ...
load_out_xdrop_from_str_ver0( ...
parameter.str_out_xdrop_a_s0000 ...
);
index_rdrop_ = xdrop_.index_rdrop_;
index_cdrop_ = xdrop_.index_cdrop_;

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we run a permutation-test. ;')); end;
if (flag_verbose); disp(sprintf(' %% This involves rerunning the biclustering-algorithm on ;')); end;
if (flag_verbose); disp(sprintf(' %% several label-permuted data-sets. ;')); end;
if (flag_verbose); disp(sprintf(' %% The parameter nshuffle indicates the (0-based) shuffle-index. ')); end;
if (flag_verbose); disp(sprintf(' %% When nshuffle==0 we do not shuffle, and process the original data-set. ;')); end;
if (flag_verbose); disp(sprintf(' %% This was already done above (as nshuffle==0 by default). ;')); end;
if (flag_verbose); disp(sprintf(' %% Each run should take as long as the original. ;')); end;
for nshuffle=1:n_shuffle;
parameter.nshuffle = nshuffle;
[ ...
 parameter ...
] = ...
xxxcluster_fromdisk_uADZSZDA_ver16( ...
 parameter ...
);
end;%for nshuffle=1:n_shuffle;

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we load the traces from each run. ;')); end;
if (flag_verbose); disp(sprintf(' %% Each trace records the (reweighted) mean-squared-correlation of the remaining array (across iterations).')); end;
trace__ = load_trace__from_dir_ver0(parameter.dir_out_trace);
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we display the original trace (red) against the label-shuffled traces (black). ;')); end;
if (flag_verbose); disp(sprintf(' %% The horizontal indicates iteration, and the vertical indicates z-score. ;')); end;
if (flag_verbose); disp(sprintf(' %% The maximum negative-log-p-value for the original trace is circled (red).')); end;
[tmp_nlpR,ij_nlpR] = max(trace__.nlpR_s0000_);
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
linewidth_sml = 0.5;
linewidth_big = 3;
markersize_big = 16;
hold on;
plot(trace__.niter_s0000_,trace__.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(trace__.niter_s0000_,trace__.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(trace__.niter_s0000_(ij_nlpR),trace__.ZR_s0000_(ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([min(trace__.niter_s0000_),max(trace__.niter_s0000_)]); xlabel('iteration');
ylabel('negative-log-p');
end;%if flag_disp;

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% We use this maximum negative-log-p-value to define a ''dominant bicluster''. ;')); end;
if (flag_verbose); disp(sprintf(' %% We extract this dominant bicluster (i.e., using the rows and cols at that iteration), ;')); end;
if (flag_verbose); disp(sprintf(' %% and then measure the principal-components of that bicluster. ;')); end;
if (flag_verbose); disp(sprintf(' %% By projecting each of the patients onto those principal-components, ;')); end;
if (flag_verbose); disp(sprintf(' %% we can visualize the effect of the biclustering. ;')); end;
if (flag_verbose); disp(sprintf(' %% A scatterplot of these patient-projections is shown below (right). ;')); end;
if (flag_verbose); disp(sprintf(' %% Each patient is colored to indicate its case-ctrl status. ;')); end;
if (flag_verbose); disp(sprintf(' %% We use our ''pca_driver'' code to calculate the principal-components (using power iteration). ;')); end;
if flag_disp;
figure(1+nf);nf=nf+1;clf;
figmed;fig80s;
end;%if flag_disp;
p_row = 1; p_col = 2; np=0;
markersize_sml = 4;
markersize_med = 8;
markersize_big = 16;
xdrop_ij_rkeep_ = xdrop_.ij_rkeep_;
xdrop_ij_ckeep_ = xdrop_.ij_ckeep_;
%%%%;
[A_p_c_,A_p_0_,AZ_rsum_] = load_A_p_c_from_dir_0(parameter.dir_out_s0000);
alpha_c_ = A_p_c_ - (1-A_p_c_);
D_c_ = sqrt(1./max(0.01,4.0*A_p_c_.*(1-A_p_c_)));
mx__ = load_mx__from_parameter_ver0(parameter);
%%%%%%%%;
ni=0;
[ ...
 AZnV_ni0_driver_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_ni_ver16( ...
 parameter ...
,ni ...
,xdrop_ ...
,trace__ ...
,mx__ ...
);
if flag_disp;
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_ni0_driver_(:,1),AZnV_ni0_driver_(:,2),16,mr_dvx_,'filled','MarkerEdgeColor','k');
axisnotick; title('full data-set (driver)'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;
%%%%%%%%;
ni=ij_nlpR-1;
[ ...
 AZnV_nix_driver_ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_ni_ver16( ...
 parameter ...
,ni ...
,xdrop_ ...
,trace__ ...
,mx__ ...
);
if flag_disp;
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_nix_driver_(:,1),AZnV_nix_driver_(:,2),16,mr_dvx_,'filled','MarkerEdgeColor','k');
axisnotick; title('bicluster-informed (driver)'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;

if (flag_verbose); disp(sprintf(' %% [finished xxxcluster_fromdisk_ver16_dr_0]')); end;






