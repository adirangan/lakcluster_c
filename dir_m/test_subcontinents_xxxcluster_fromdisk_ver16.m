
verbose=1;
if (verbose); disp(sprintf(' %% [entering test_subcontinents_xxxcluster_fromdisk_ver16]')); end;
setup_0;
rng(0); nf=0;

if (verbose); disp(sprintf(' %% ;')); end;
if (verbose); disp(sprintf(' %% First we generate a very simple example of dosage-data ;')); end;
if (verbose); disp(sprintf(' %% which exhibits two large clumps (i.e., subcontinents), ;')); end;
if (verbose); disp(sprintf(' %% each with a different trend separating the cases from controls. ;')); end;
str_lak_vs_dex = 'dex';
r_fraction_scA = 0.45;
r_fraction_scB = 1 - r_fraction_scA ;
c_fraction_scA = 0.1250;
c_fraction_scB = 0.1250;
maf_lob = 0.10;
maf_upb = 1.00 - maf_lob;
flag_force_create = 0;

%%%%%%%%;
n_snp_cup = 1024*2;
bim_cup_khr_ = zeros(n_snp_cup,1);
bim_cup_vid_ = cell(n_snp_cup,1);
bim_cup_pos_ = zeros(n_snp_cup,1);
bim_cup_bpc_ = zeros(n_snp_cup,1);
bim_cup_al1_ = cell(n_snp_cup,1);
bim_cup_al2_ = cell(n_snp_cup,1);
tmp_al1_ = {'A','T','G','C'};
tmp_al2_ = {'G','C','A','T'};
for nsnp_cup=0:n_snp_cup-1;
bim_cup_khr_(1+nsnp_cup) = 1+mod(nsnp_cup,22);
bim_cup_vid_{1+nsnp_cup} = sprintf('vid%.4d',1+nsnp_cup);
bim_cup_pos_(1+nsnp_cup) = 1+nsnp_cup;
bim_cup_bpc_(1+nsnp_cup) = 1+nsnp_cup;
bim_cup_al1_{1+nsnp_cup} = sprintf('%s',tmp_al1_{1+mod(nsnp_cup,4)});
bim_cup_al2_{1+nsnp_cup} = sprintf('%s',tmp_al2_{1+mod(nsnp_cup,4)});
end;%for nsnp_cup=0:n_snp_cup-1;
cij_ = transpose(1:n_snp_cup);
cij_tmp_ = cij_;
cij_scA_ = cij_tmp_(randperm(numel(cij_tmp_),floor(n_snp_cup*c_fraction_scA)));
n_snp_scA = numel(cij_scA_);
cij_tmp_ = setdiff(cij_,cij_scA_);
cij_scB_ = cij_tmp_(randperm(numel(cij_tmp_),floor(n_snp_cup*c_fraction_scB)));
n_snp_scB = numel(cij_scB_);
%%%%%%%%;
n_patient_cup = 1024*3;
%%%%;
n_mds_p = n_patient_cup;
n_mds_v = 2;
mds_pv__ = randn(n_mds_p,n_mds_v);
mds_r_p_ = sqrt( mds_pv__(:,1+0).^2 + mds_pv__(:,1+1).^2 );
mds_w_p_ = atan2( mds_pv__(:,1+1) , mds_pv__(:,1+0) );
%%%%;
fam_cup_fid_ = cell(n_patient_cup,1);
fam_cup_iid_ = cell(n_patient_cup,1);
fam_cup_yid_ = cell(n_patient_cup,1);
fam_cup_xid_ = cell(n_patient_cup,1);
fam_cup_sex_ = zeros(n_patient_cup,1);
fam_cup_dvx_ = zeros(n_patient_cup,1);
fam_cup_fidandiid_ = cell(n_patient_cup,1);
tmp_sex_ = [1,2,1,2]; %<-- 1==male, 2==fema, ignored. ;
tmp_dvx_ = [1,1,1,2,2]; %<-- 1==ctrl, 2==case, not ignored. ;
for npatient_cup=0:n_patient_cup-1;
fam_cup_fid_{1+npatient_cup} = sprintf('fid%.4d',1+npatient_cup);
fam_cup_iid_{1+npatient_cup} = sprintf('iid%.4d',1+npatient_cup);
fam_cup_yid_{1+npatient_cup} = sprintf('yid%.4d',1+npatient_cup);
fam_cup_xid_{1+npatient_cup} = sprintf('xid%.4d',1+npatient_cup);
fam_cup_sex_(1+npatient_cup) = tmp_sex_(1+mod(npatient_cup,4));
fam_cup_dvx_(1+npatient_cup) = tmp_dvx_(1+floor(5*npatient_cup/n_patient_cup)); %<-- first 3/5 of patients are ctrls. ;
end;%for npatient_cup=0:n_patient_cup-1;
for npatient_cup=0:n_patient_cup-1;
fam_cup_fidandiid_{1+npatient_cup} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for npatient_cup=0:n_patient_cup-1;
rij_case_ = find(fam_cup_dvx_==2); rij_ctrl_ = find(fam_cup_dvx_==1); 
n_patient_case = numel(rij_case_); n_patient_ctrl = numel(rij_ctrl_);
rij_case_tmp_ = rij_case_;
rij_case_scA_ = rij_case_tmp_(randperm(numel(rij_case_tmp_),floor(r_fraction_scA*n_patient_case)));
rij_case_scB_ = setdiff(rij_case_tmp_,rij_case_scA_);
n_patient_case_scA = numel(rij_case_scA_);
n_patient_case_scB = numel(rij_case_scB_);
rij_ctrl_tmp_ = rij_ctrl_;
rij_ctrl_scA_ = rij_ctrl_tmp_(randperm(numel(rij_ctrl_tmp_),floor(r_fraction_scA*n_patient_ctrl)));
rij_ctrl_scB_ = setdiff(rij_ctrl_tmp_,rij_ctrl_scA_);
n_patient_ctrl_scA = numel(rij_ctrl_scA_);
n_patient_ctrl_scB = numel(rij_ctrl_scB_);
rij_case_scA_cmp_ = setdiff(rij_case_,rij_case_scA_); n_patient_case_scA_cmp = numel(rij_case_scA_cmp_);
rij_case_scB_cmp_ = setdiff(rij_case_,rij_case_scB_); n_patient_case_scB_cmp = numel(rij_case_scB_cmp_);
rij_ctrl_scB_cmp_ = setdiff(rij_ctrl_,rij_ctrl_scB_); n_patient_ctrl_scB_cmp = numel(rij_ctrl_scB_cmp_);
%%%%%%%%;
frq_s_ = linspace(maf_lob,maf_upb,n_snp_cup);
ds1_cup_ps__ = zeros(n_patient_cup,n_snp_cup);
ds1_cup_ps__ = rand(n_patient_cup,n_snp_cup)<repmat(frq_s_,[n_patient_cup,1]);
ds2_cup_ps__ = zeros(n_patient_cup,n_snp_cup);
ds2_cup_ps__ = rand(n_patient_cup,n_snp_cup)<repmat(frq_s_,[n_patient_cup,1]);
%%%%%%%%;
if (strcmp(str_lak_vs_dex,'lak'));
if (verbose); disp(sprintf(' %% Warning, str_lak_vs_dex %s not implemented',str_lak_vs_dex)); end;
end;%if (strcmp(str_lak_vs_dex,'lak'));
%%%%%%%%;
if (strcmp(str_lak_vs_dex,'dex'));
%%%%;
% now magnify diplotypes within scA. ;
%%%%;
maf_upb_use = 1.0;
maf_lob_use = 0.0;
f_case = 0.65;
f_ctrl = 0.35;
frq_case_scA_s_ = (frq_s_ >= 0.5) .* (frq_s_ + f_case*(maf_upb_use-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + f_case*(maf_lob_use-frq_s_)) ;
ds1_cup_ps__(rij_case_scA_,cij_scA_) = rand(n_patient_case_scA,n_snp_scA) < repmat(frq_case_scA_s_(cij_scA_),[n_patient_case_scA,1]);
ds1_cup_ps__(rij_case_scA_,cij_scA_) = ds1_cup_ps__(rij_case_scA_(randperm(n_patient_case_scA)),cij_scA_);
ds2_cup_ps__(rij_case_scA_,cij_scA_) = rand(n_patient_case_scA,n_snp_scA) < repmat(frq_case_scA_s_(cij_scA_),[n_patient_case_scA,1]);
ds2_cup_ps__(rij_case_scA_,cij_scA_) = ds2_cup_ps__(rij_case_scA_(randperm(n_patient_case_scA)),cij_scA_);
frq_ctrl_scA_s_ = (frq_s_ >= 0.5) .* (frq_s_ + f_ctrl*(maf_upb_use-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + f_ctrl*(maf_lob_use-frq_s_)) ;
ds1_cup_ps__(rij_ctrl_scA_,cij_scA_) = rand(n_patient_ctrl_scA,n_snp_scA) < repmat(frq_ctrl_scA_s_(cij_scA_),[n_patient_ctrl_scA,1]);
ds1_cup_ps__(rij_ctrl_scA_,cij_scA_) = ds1_cup_ps__(rij_ctrl_scA_(randperm(n_patient_ctrl_scA)),cij_scA_);
ds2_cup_ps__(rij_ctrl_scA_,cij_scA_) = rand(n_patient_ctrl_scA,n_snp_scA) < repmat(frq_ctrl_scA_s_(cij_scA_),[n_patient_ctrl_scA,1]);
ds2_cup_ps__(rij_ctrl_scA_,cij_scA_) = ds2_cup_ps__(rij_ctrl_scA_(randperm(n_patient_ctrl_scA)),cij_scA_);
%%%%;
% now magnify diplotypes within scB. ;
%%%%;
f_case = 0.65;
f_ctrl = 0.35;
frq_case_scB_s_ = (frq_s_ >= 0.5) .* (frq_s_ + f_case*(maf_upb_use-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + f_case*(maf_lob_use-frq_s_)) ;
ds1_cup_ps__(rij_case_scB_,cij_scB_) = rand(n_patient_case_scB,n_snp_scB) < repmat(frq_case_scB_s_(cij_scB_),[n_patient_case_scB,1]);
ds1_cup_ps__(rij_case_scB_,cij_scB_) = ds1_cup_ps__(rij_case_scB_(randperm(n_patient_case_scB)),cij_scB_);
ds2_cup_ps__(rij_case_scB_,cij_scB_) = rand(n_patient_case_scB,n_snp_scB) < repmat(frq_case_scB_s_(cij_scB_),[n_patient_case_scB,1]);
ds2_cup_ps__(rij_case_scB_,cij_scB_) = ds2_cup_ps__(rij_case_scB_(randperm(n_patient_case_scB)),cij_scB_);
frq_ctrl_scB_s_ = (frq_s_ >= 0.5) .* (frq_s_ + f_ctrl*(maf_upb_use-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + f_ctrl*(maf_lob_use-frq_s_)) ;
ds1_cup_ps__(rij_ctrl_scB_,cij_scB_) = rand(n_patient_ctrl_scB,n_snp_scB) < repmat(frq_ctrl_scB_s_(cij_scB_),[n_patient_ctrl_scB,1]);
ds1_cup_ps__(rij_ctrl_scB_,cij_scB_) = ds1_cup_ps__(rij_ctrl_scB_(randperm(n_patient_ctrl_scB)),cij_scB_);
ds2_cup_ps__(rij_ctrl_scB_,cij_scB_) = rand(n_patient_ctrl_scB,n_snp_scB) < repmat(frq_ctrl_scB_s_(cij_scB_),[n_patient_ctrl_scB,1]);
ds2_cup_ps__(rij_ctrl_scB_,cij_scB_) = ds2_cup_ps__(rij_ctrl_scB_(randperm(n_patient_ctrl_scB)),cij_scB_);
%%%%;
end;%if (strcmp(str_lak_vs_dex,'dex'));
%%%%%%%%;
dsg_cup_ps__ = ds1_cup_ps__ + ds2_cup_ps__;
dsg_cup_frq_mss_ = mean(dsg_cup_ps__>=3,1); tmp_z_ = 1 - dsg_cup_frq_mss_;
dsg_cup_frq_nor_ = mean(dsg_cup_ps__==0,1)./tmp_z_;
dsg_cup_frq_xor_ = mean(dsg_cup_ps__==1,1)./tmp_z_;
dsg_cup_frq_and_ = mean(dsg_cup_ps__==2,1)./tmp_z_;
dsg_cup_p_opt_ = dsg_cup_frq_and_ + 0.5*dsg_cup_frq_xor_ ;
dsg_cup_q_opt_ = dsg_cup_frq_nor_ + 0.5*dsg_cup_frq_xor_ ;
tmp_and_ = dsg_cup_frq_and_.*log(dsg_cup_frq_and_./(dsg_cup_p_opt_.^2)); tmp_and_(find(~isfinite(tmp_and_)))=0;
tmp_xor_ = dsg_cup_frq_xor_.*log(dsg_cup_frq_xor_./(2*dsg_cup_p_opt_.*dsg_cup_q_opt_)); tmp_xor_(find(~isfinite(tmp_xor_)))=0;
tmp_nor_ = dsg_cup_frq_nor_.*log(dsg_cup_frq_nor_./(dsg_cup_q_opt_.^2)); tmp_nor_(find(~isfinite(tmp_nor_)))=0;
dsg_cup_I_opt_ = tmp_and_ + tmp_xor_ + tmp_nor_ ;
%%%%%%%%;
rij_case_prm_ = [ ... 
; setdiff(rij_case_scA_,unionall({rij_case_scB_})) ...
; setdiff(rij_case_scB_,unionall({})) ...
; setdiff(rij_case_,unionall({rij_case_scA_,rij_case_scB_})) ...
];
rij_ctrl_prm_ = [ ... 
; setdiff(rij_ctrl_scA_,unionall({rij_ctrl_scB_})) ...
; setdiff(rij_ctrl_scB_,unionall({})) ...
; setdiff(rij_ctrl_,unionall({rij_ctrl_scA_,rij_ctrl_scB_})) ...
];
rij_prm_ = [ rij_case_prm_ ; rij_ctrl_prm_ ];
cij_prm_ = [ ...
; setdiff(cij_scA_,unionall({cij_scB_})) ...
; setdiff(cij_scB_,unionall({})) ...
; setdiff(cij_,unionall({cij_scA_,cij_scB_})) ...
];
%%%%%%%%;
mr_scA_ = zeros(n_patient_cup,1); mr_scA_(union(rij_case_scA_,rij_ctrl_scA_)) = 1;
mr_scB_ = zeros(n_patient_cup,1); mr_scB_(union(rij_case_scB_,rij_ctrl_scB_)) = 1;
mc_cup_scA_ = zeros(n_snp_cup,1); mc_cup_scA_(cij_scA_) = 1;
mc_cup_scB_ = zeros(n_snp_cup,1); mc_cup_scB_(cij_scB_) = 1;
%%%%%%%;

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we visualize the dosage-data (magenta and cyan) for the cases (top) and ctrls (bottom). ;')); end;
if (verbose); disp(sprintf(' %% The first two mds-components are shown on the far left column (green and yellow). ;')); end;
if (verbose); disp(sprintf(' %% Adjacent are the row-masks associated with subcontinents scA, scB (grey). ;')); end;
if (verbose); disp(sprintf(' %% Atop lie the column-masks associated with subcontinents scA, scB (grey). ;')); end;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,512,1024]); 
imagesc_uAZm_1( ...
 mds_pv__(rij_prm_,:) ...
,dsg_cup_ps__(rij_case_prm_,cij_prm_) ...
,dsg_cup_ps__(rij_ctrl_prm_,cij_prm_) ...
,[ mr_scA_(rij_prm_) , mr_scB_(rij_prm_) ] ...
,[ mc_cup_scA_(cij_prm_) , mc_cup_scB_(cij_prm_) ] ...
);

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we convert the dosage-data to b16-data. ;')); end;
if (verbose); disp(sprintf(' %% This encodes each dosage value as a 3-bit binary array. ;')); end;
if (verbose); disp(sprintf(' %% dosage ''0'' is encoded as [ 1 0 0 ]. ;')); end;
if (verbose); disp(sprintf(' %% dosage ''1'' is encoded as [ 0 1 0 ]. ;')); end;
if (verbose); disp(sprintf(' %% dosage ''2'' is encoded as [ 0 0 1 ]. ;')); end;
if (verbose); disp(sprintf(' %% For simplicity, we put everything into a single study directory. ;')); end;
if (verbose); disp(sprintf(' %% (although the dosage-data could instead be split across multiple files). ;')); end;
%%%%%%%%;
dir_trunk = sprintf('%s/dir_test_subcontinents_xxxcluster_fromdisk_ver16',pwd);
if (~exist(dir_trunk,'dir')); disp(sprintf(' %% mkdir %s',dir_trunk)); mkdir(dir_trunk); end;
n_study = 1; nstudy = 0;
study_name = sprintf('study%.2d',nstudy);
study_trunk = sprintf('dir_%s',study_name);
tmp_dir = sprintf('%s/%s',dir_trunk,study_trunk);
if (~exist(tmp_dir,'dir')); disp(sprintf(' %% mkdir %s',tmp_dir)); mkdir(tmp_dir); end;
n_snp = n_snp_cup;
index_nsnp_from_nstudy_ = 0:n_snp-1;
n_patient = n_patient_cup;
index_npatient_from_nstudy_ = 0:n_patient-1;
%%%%%%%%;
fn_fam = sprintf('%s/%s/%s.fam',dir_trunk,study_trunk,study_name);
if (flag_force_create | ~exist(fn_fam,'file'));
fid = fopen(fn_fam,'w');
for npatient=0:n_patient-1;
np = index_npatient_from_nstudy_(1+npatient);
fprintf( ...
 fid ...
,'%s\t%s\t%s\t%s\t%d\t%d\n' ...
,fam_cup_fid_{1+np} ...
,fam_cup_iid_{1+np} ...
,fam_cup_yid_{1+np} ...
,fam_cup_xid_{1+np} ...
,fam_cup_sex_(1+np) ...
,fam_cup_dvx_(1+np) ...
);
end;%for npatient=0:n_patient-1;
fclose(fid);
end;%if (flag_force_create | ~exist(fn_fam,'file'));
%%%%%%%%;
fn_bim = sprintf('%s/%s/%s.bim',dir_trunk,study_trunk,study_name);
if (flag_force_create | ~exist(fn_bim,'file'));
fid = fopen(fn_bim,'w');
for nsnp=0:n_snp-1;
ns = index_nsnp_from_nstudy_(1+nsnp);
fprintf( ...
 fid ...
,'%d\t%s\t%d\t%d\t%s\t%s\n' ...
,bim_cup_khr_(1+ns) ...
,bim_cup_vid_{1+ns} ...
,bim_cup_pos_(1+ns) ...
,bim_cup_bpc_(1+ns) ...
,bim_cup_al1_{1+ns} ...
,bim_cup_al2_{1+ns} ...
);
end;%for nsnp=0:n_snp-1;
fclose(fid);
end;%if (flag_force_create | ~exist(fn_bim,'file'));
%%%%%%%%;
fn_bed = sprintf('%s/%s/%s.bed',dir_trunk,study_trunk,study_name);
if (flag_force_create | ~exist(fn_bed,'file'));
fid = fopen(fn_bed,'w');
tmp_dsg_ps__ = dsg_cup_ps__(1+index_npatient_from_nstudy_,1+index_nsnp_from_nstudy_);
[tmp_bed_ps__,n_patient_bed] = dsg_to_bed_0(n_patient,n_snp,tmp_dsg_ps__);
key1 = uint8(108); key2 = uint8(27); key3 = uint8(1);
fwrite(fid,key1,'uint8'); fwrite(fid,key2,'uint8'); fwrite(fid,key3,'uint8');
fwrite(fid,tmp_bed_ps__,'uint8');
fclose(fid);
end;%if (flag_force_create | ~exist(fn_bed,'file'));
%%%%%%%%;
mds_fidandiid_p_ = cell(n_mds_p,1);
for npatient_cup=0:n_patient_cup-1;
mds_fidandiid_p_{1+npatient_cup} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for npatient_cup=0:n_patient_cup-1;
%%%%%%%%;
n_famex = 0;
index_npatient_cup_from_nfamex_ = randperm(n_patient_cup,n_famex)-1;
famex_fidandiid_p_ = cell(n_famex,1);
for nfamex=0:n_famex-1;
npatient_cup = index_npatient_cup_from_nfamex_(1+nfamex);
famex_fidandiid_p_{1+nfamex} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for nfamex=0:n_famex-1;
%%%%%%%%;
if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% When converting the bed, bim and fam data into a b16 array, ;')); end;
if (verbose); disp(sprintf(' %% we must specify a few parameters. ;')); end;
if (verbose); disp(sprintf(' %% Some of the more important parameter are: ;')); end;
if (verbose); disp(sprintf(' %% maf_cutoff: scalar double: this is the minor-allele-frequency (maf) cutoff for snps. ;')); end;
if (verbose); disp(sprintf(' %%\t Only snps with a strictly higher maf will be considered. ;')); end;
if (verbose); disp(sprintf(' %% ent_cutoff: scalar double: this is the entropy-cutoff for snps. ;')); end;
if (verbose); disp(sprintf(' %%\t We typically expect the dosage-data for any particular snp ;')); end;
if (verbose); disp(sprintf(' %%\t to have a frequency of 0, 1 and 2 values ;')); end;
if (verbose); disp(sprintf(' %%\t (i.e., ''aa'', ''Aa'' and ''AA'' values) ;')); end;
if (verbose); disp(sprintf(' %%\t which is close to pp, 2pq and qq, respectively, ;')); end;
if (verbose); disp(sprintf(' %%\t where p is the associated minor-allele-frequency (maf). ;')); end;
if (verbose); disp(sprintf(' %%\t The entropy (I_opt) measures the deviation from this expectation, ;')); end;
if (verbose); disp(sprintf(' %%\t with higher values indicating a greater deviation. ;')); end;
if (verbose); disp(sprintf(' %%\t Only snps with a strictly lower I_opt will be considered.  ;')); end;
str_output_prefix = 'test';
parameter = struct('type','parameter');
parameter.verbose = 1;
parameter.str_output_prefix = str_output_prefix;
parameter.dir_trunk = dir_trunk;
parameter.ent_cutoff = 0.001; %<-- this is a typical entropy-cutoff. ;
parameter.ent_cutoff = prctile(dsg_cup_I_opt_,99); %<-- for this example we will make sure to include most of the snps. ;
parameter.maf_cutoff = 0.010;
parameter.flag_force_create = flag_force_create;
[ ...
 parameter ...
] = ...
bed_to_b16_flip_ver7( ...
 parameter ...
,n_study ...
,{study_trunk} ...
,{study_name} ...
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
if (verbose>0); disp(sprintf(' %% str_output_prefix_local %s ;\n %% dir_out %s ;\n',str_output_prefix_local,dir_out)); end;

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now for illustration we uncompress the b16 file. ;')); end;
if (verbose); disp(sprintf(' %% This is definitely not recommended for larger data-sets. ;')); end;
fname_b16 = sprintf('%s/%s_A_full_n.b16',dir_out,str_output_prefix_plus_maf);
A_full_n_ = binary_uncompress(fname_b16);
%%%%%%%%;
fname_famext = sprintf('%s/%s_fam.ext',dir_out,str_output_prefix_plus_maf);
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
[~,~,ij_fam_ext_from_cap_] = intersect(fam_cup_fid_,famext_fid_,'stable');
disp(sprintf(' %% famext_fid_ vs fam_cup_fid_: %0.16f',fnorm(ij_fam_ext_from_cap_-transpose(1:n_patient))));
%%%%%%%%;
fname_bimext = sprintf('%s/%s_bim.ext',dir_out,str_output_prefix_plus_maf);
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
[bim_vid_cap_,index_bim_cup_from_cap_,index_bim_u_from_cap_] = intersect(bim_cup_vid_,u_bimext_vid_,'stable');
index_bim_cup_from_cap_ = index_bim_cup_from_cap_ - 1;
index_bim_u_from_cap_ = index_bim_u_from_cap_ - 1;
[~,~,index_bim_cap_from_u_] = intersect(u_bimext_vid_,bim_vid_cap_,'stable');
index_bim_cap_from_u_ = index_bim_cap_from_u_ - 1;
index_bim_ext_from_cap_ = index_bim_ext_from_u_(1+index_bim_u_from_cap_);
index_bim_cup_from_ext_ = index_bim_cup_from_cap_(1+index_bim_cap_from_u_(1+index_bim_u_from_ext_));
ij_bim_cup_from_ext__ = sparse(1+index_bim_cup_from_ext_,1:n_snp_bimext,1,n_snp_cup,n_snp_bimext);
[~,cij_prm_inv_] = sort(cij_prm_,'ascend');
cij_prm_inv_from_ext_ = transpose(cij_prm_inv_) * ij_bim_cup_from_ext__;
[~,cij_ext_] = sort(cij_prm_inv_from_ext_,'ascend');
%%%%;
mc_cup_scA_ = zeros(n_snp_cup,1);
mc_cup_scA_(cij_scA_) = 1;
mc_cup_scA_from_ext_ = transpose(mc_cup_scA_) * ij_bim_cup_from_ext__;
n_snp_ext_scA = sum(mc_cup_scA_from_ext_);
[~,ij_mc_cup_scA_from_ext_] = sort(mc_cup_scA_from_ext_,'ascend');
mc_ext_scA_ = zeros(n_snp_bimext,1);
mc_ext_scA_(ij_mc_cup_scA_from_ext_(end-n_snp_ext_scA+1:end)) = 1;
%%%%;
mc_cup_scB_ = zeros(n_snp_cup,1);
mc_cup_scB_(cij_scB_) = 1;
mc_cup_scB_from_ext_ = transpose(mc_cup_scB_) * ij_bim_cup_from_ext__;
n_snp_ext_scB = sum(mc_cup_scB_from_ext_);
[~,ij_mc_cup_scB_from_ext_] = sort(mc_cup_scB_from_ext_,'ascend');
mc_ext_scB_ = zeros(n_snp_bimext,1);
mc_ext_scB_(ij_mc_cup_scB_from_ext_(end-n_snp_ext_scB+1:end)) = 1;
%%%%;
cij_ex2_ = ...
[ ...
 find(mc_ext_scA_) ...
; find(mc_ext_scB_) ...
; setdiff(transpose(1:n_snp_bimext),find(mc_ext_scA_ + mc_ext_scB_)) ...
];
%%%%;
mr_bcx_ = zeros(n_patient,1);
mr_bcx_(rij_case_scB_) = 1;
mr_bcx_(rij_ctrl_scB_) = 2;
mr_bcx_(rij_case_scA_) = 3;
mr_bcx_(rij_ctrl_scA_) = 4;
%%%%;
mc_bcx_ = zeros(n_snp_bimext,1);
mc_bcx_(find(mc_ext_scA_))=4;
mc_bcx_(find(mc_ext_scB_))=2;
n_snp_bcx = numel(find(mc_bcx_));
n_snp_ext_scA = numel(find(mc_bcx_==4));
n_snp_ext_scB = numel(find(mc_bcx_==2));
%%%%;
if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we visualize the b16-data (magenta and cyan) for the cases (top) and ctrls (bottom). ;')); end;
%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,1024]); fig80s;
imagesc_uAZm_1( ...
 mds_pv__(rij_prm_,:) ...
,A_full_n_(rij_case_prm_,cij_ex2_) ...
,A_full_n_(rij_ctrl_prm_,cij_ex2_) ...
,[ mr_scA_(rij_prm_) , mr_scB_(rij_prm_) ] ...
,[ mc_ext_scA_(cij_ex2_) , mc_ext_scB_(cij_ex2_) ] ...
);
%%%%;

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we run the basic biclustering method on the b16-data. ;')); end;
if (verbose); disp(sprintf(' %% This should take approximately 10 seconds for differentially-expressed biclustering. ;')); end;
if (verbose); disp(sprintf(' %% Some of the important parameters are: ;')); end;
if (verbose); disp(sprintf(' %% n_mds: scalar integer: number of mds-components to correct for. ;')); end;
if (verbose); disp(sprintf(' %%\t Note that the mds-components will be binarized into orthants. ;')); end;
if (verbose); disp(sprintf(' %% n_mds_repl: scalar integer: number of times to replicate the mds-components. ;')); end;
if (verbose); disp(sprintf(' %%\t This specifically refers to rotation+replication across orthants. ;')); end;
if (verbose); disp(sprintf(' %%\t Larger numbers are better, but 1 or 2 is usually sufficient. ;')); end;
if (verbose); disp(sprintf(' %% ij_mds_use_: integer array of size n_mds: 1-based indices of mds-components to correct for. ;')); end;
if (verbose); disp(sprintf(' %% gamma: scalar double: fraction of rows or columns eliminated at each iteration of the bicluster algorithm. ;')); end;
if (verbose); disp(sprintf(' %%\t Smaller numbers are more accurate, but larger numbers require less computation time. ;')); end;
if (verbose); disp(sprintf(' %%\t Empirically, gamma = 0.01 or 0.05 is usually sufficient. ;')); end;
if (verbose); disp(sprintf(' %% flag_force_create: scalar integer: set to 1 to force recreation of all temporary files. ;')); end;
if (verbose); disp(sprintf(' %%\t if set to 0 will load temporary files from disk. ;')); end;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
gamma = 0.05;
dir_code = sprintf('%s/..',pwd);
str_prefix = str_output_prefix_plus_maf;
parameter.str_lak_vs_dex = str_lak_vs_dex;
parameter.str_prefix = str_prefix;
parameter.dir_code = dir_code;
parameter.gamma = gamma;
parameter.n_mds = n_mds_0in;
parameter.n_mds_repl = n_mds_repl;
parameter.ij_mds_use_ = ij_mds_use_;
parameter.flag_force_create = flag_force_create;
[ ...
 parameter ...
] = ...
xxxcluster_fromdisk_uADZSZDA_ver16( ...
 parameter ...
);

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we load the results of this single biclustering run. ;')); end;
parameter.str_out_xdrop_a_s0000 = sprintf('%s/out_xdrop_a.txt',parameter.dir_out_s0000);
[ ...
 xdrop_ ...
] = ...
load_out_xdrop_from_str_ver0( ...
parameter.str_out_xdrop_a_s0000 ...
);
index_rdrop_ = xdrop_.index_rdrop_;
index_cdrop_ = xdrop_.index_cdrop_;

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we display the results of this single biclustering run. ;')); end;
if (verbose); disp(sprintf(' %% The rows and cols are ordered by the biclustering-algorithm. ;')); end;
if (verbose); disp(sprintf(' %% Those rows and cols eliminated first are shown towards the bottom-right. ;')); end;
if (verbose); disp(sprintf(' %% Those rows and cols retained longest are shown towards the top-left. ;')); end;
if (verbose); disp(sprintf(' %% Note that most of scA is retained until the end. ;')); end;
if (verbose); disp(sprintf(' %% That is to say, the biclustering-algorithm successfully focused on scA. ;')); end;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,1024]); fig80s;
tmp_rij_prm_ = [xdrop_.ij_rkeep_;rij_ctrl_prm_];
tmp_cij_prm_ = [xdrop_.ij_ckeep_];
imagesc_uAZm_1( ...
 mds_pv__(tmp_rij_prm_,:) ...
,A_full_n_(xdrop_.ij_rkeep_,xdrop_.ij_ckeep_) ...
,A_full_n_(rij_ctrl_prm_,xdrop_.ij_ckeep_) ...
,[ mr_scA_(tmp_rij_prm_) , mr_scB_(tmp_rij_prm_) ] ...
,[ mc_ext_scA_(tmp_cij_prm_) , mc_ext_scB_(tmp_cij_prm_) ] ...
);

mr_ABC__ = zeros(n_patient,3);
mr_ABC__(rij_case_scA_,1+0) = 1.0;
mr_ABC__(rij_ctrl_scA_,1+0) = 0.5;
mr_ABC__(rij_case_scB_,1+1) = 1.0;
mr_ABC__(rij_ctrl_scB_,1+1) = 0.5;
%mr_ABC__(rij_case_bcC_,1+2) = 1;
A_p_ = mean(A_full_n_>0,1);
[A_p_c_,A_p_0_,AZ_rsum_] = load_A_p_c_from_dir_0(parameter.dir_out_s0000);
alpha_c_ = A_p_c_ - (1-A_p_c_);
D_c_ = sqrt(1./max(0.01,4.0*A_p_c_.*(1-A_p_c_)));
A_full_normalized_n_ = bsxfun(@times,bsxfun(@minus,A_full_n_,transpose(alpha_c_)),transpose(sqrt(D_c_)));
[tmp_U_,tmp_S_,tmp_V_] = svds(A_full_normalized_n_,2);
AZnV_ni0_ = A_full_normalized_n_*tmp_V_;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
scatter(AZnV_ni0_(:,1),AZnV_ni0_(:,2),16,mr_ABC__,'filled','MarkerEdgeColor','k');
axisnotick; title('full data-set (matlab)'); xlabel('PC1'); ylabel('PC2');

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we run a permutation-test. ;')); end;
if (verbose); disp(sprintf(' %% This involves rerunning the biclustering-algorithm on ;')); end;
if (verbose); disp(sprintf(' %% several label-permuted data-sets. ;')); end;
if (verbose); disp(sprintf(' %% The parameter nshuffle indicates the (0-based) shuffle-index. ')); end;
if (verbose); disp(sprintf(' %% When nshuffle==0 we do not shuffle, and process the original data-set. ;')); end;
if (verbose); disp(sprintf(' %% This was already done above (as nshuffle==0 by default). ;')); end;
if (verbose); disp(sprintf(' %% Each run should take as long as the original. ;')); end;
n_shuffle = 32;
for nshuffle=1:n_shuffle;
n_mds_0in = 2; n_mds_repl = 1; ij_mds_use_ = [1:2];
gamma = 0.05;
dir_code = sprintf('%s/..',pwd);
str_prefix = str_output_prefix_plus_maf;
parameter.str_lak_vs_dex = str_lak_vs_dex;
parameter.str_prefix = str_prefix;
parameter.dir_code = dir_code;
parameter.gamma = gamma;
parameter.n_mds = n_mds_0in;
parameter.n_mds_repl = n_mds_repl;
parameter.ij_mds_use_ = ij_mds_use_;
parameter.flag_force_create = flag_force_create;
parameter.nshuffle = nshuffle;
[ ...
 parameter ...
] = ...
xxxcluster_fromdisk_uADZSZDA_ver16( ...
 parameter ...
);
end;%for nshuffle=1:n_shuffle;

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we load the traces from each run. ;')); end;
if (verbose); disp(sprintf(' %% Each trace records the (reweighted) mean-squared-correlation of the remaining array (across iterations).')); end;
trace__ = load_trace__from_dir_ver0(parameter.dir_out_trace);
if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we display the original trace (red) against the label-shuffled traces (black). ;')); end;
if (verbose); disp(sprintf(' %% The horizontal indicates iteration, and the vertical indicates z-score. ;')); end;
if (verbose); disp(sprintf(' %% The maximum negative-log-p-value for the original trace is circled (red).')); end;
figure(1+nf);nf=nf+1;clf;figsml;
linewidth_sml = 0.5;
linewidth_big = 3;
markersize_big = 16;
[tmp_nlpR,ij_nlpR] = max(trace__.nlpR_s0000_);
hold on;
plot(trace__.niter_s0000_,trace__.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(trace__.niter_s0000_,trace__.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(trace__.niter_s0000_(ij_nlpR),trace__.ZR_s0000_(ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([min(trace__.niter_s0000_),max(trace__.niter_s0000_)]); xlabel('iteration');
ylabel('negative-log-p');

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% We use this maximum negative-log-p-value to define a ''dominant bicluster''. ;')); end;
if (verbose); disp(sprintf(' %% We extract this dominant bicluster (i.e., using the rows and cols at that iteration), ;')); end;
if (verbose); disp(sprintf(' %% and then measure the principal-components of that bicluster. ;')); end;
if (verbose); disp(sprintf(' %% By projecting each of the patients onto those principal-components, ;')); end;
if (verbose); disp(sprintf(' %% we can visualize the effect of the biclustering. ;')); end;
if (verbose); disp(sprintf(' %% A scatterplot of these patient-projections is shown below (right). ;')); end;
if (verbose); disp(sprintf(' %% Each patient is colored with a (r,g,b) triplet indicating its membership in (scA,scB). ;')); end;
if (verbose); disp(sprintf(' %% Thus, patients in scA will be: ')); end;
if (verbose); disp(sprintf(' %%\t bright red if they are cases in scA, ')); end;
if (verbose); disp(sprintf(' %%\t dark red if they are ctrls in scA, ')); end;
if (verbose); disp(sprintf(' %%\t bright green if they are cases in scB, ')); end;
if (verbose); disp(sprintf(' %%\t dark green if they are ctrls in scB, ')); end;
if (verbose); disp(sprintf(' %% The top two subplots use matlab to calculate principal-components (using the svd). ;')); end;
if (verbose); disp(sprintf(' %% The bottom two subplots use our ''pca_driver'' code. ;')); end;
if (verbose); disp(sprintf(' %% These results are not exactly the same since our driver does not operate at double precision. ;')); end;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024]);fig80s;
p_row = 2; p_col = 2; np=0;
mr_ABC__ = zeros(n_patient,3);
mr_ABC__(rij_case_scA_,1+0) = 1.0;
mr_ABC__(rij_ctrl_scA_,1+0) = 0.5;
mr_ABC__(rij_case_scB_,1+1) = 1.0;
mr_ABC__(rij_ctrl_scB_,1+1) = 0.5;
%mr_ABC__(rij_case_bcC_,1+2) = 1;
markersize_sml = 4;
markersize_med = 8;
markersize_big = 16;
xdrop_ij_rkeep_ = xdrop_.ij_rkeep_;
xdrop_ij_ckeep_ = xdrop_.ij_ckeep_;
%%%%;
A_p_ = mean(A_full_n_>0,1);
[A_p_c_,A_p_0_,AZ_rsum_] = load_A_p_c_from_dir_0(parameter.dir_out_s0000);
alpha_c_ = A_p_c_ - (1-A_p_c_);
D_c_ = sqrt(1./max(0.01,4.0*A_p_c_.*(1-A_p_c_)));
A_full_normalized_n_ = bsxfun(@times,bsxfun(@minus,A_full_n_,transpose(alpha_c_)),transpose(sqrt(D_c_)));
[tmp_U_,tmp_S_,tmp_V_] = svds(A_full_normalized_n_(xdrop_ij_rkeep_,xdrop_ij_ckeep_),2);
AZnV_ni0_ = A_full_normalized_n_(:,xdrop_ij_ckeep_)*tmp_V_;
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_ni0_(:,1),AZnV_ni0_(:,2),16,mr_ABC__,'filled','MarkerEdgeColor','k');
axisnotick; title('full data-set (matlab)'); xlabel('PC1'); ylabel('PC2');
%%%%;
n_r_tmp = trace__.r_rem_s0000_(ij_nlpR);
n_c_tmp = trace__.c_rem_s0000_(ij_nlpR);
[tmp_U_,tmp_S_,tmp_V_] = svds(A_full_normalized_n_(xdrop_ij_rkeep_(1:n_r_tmp),xdrop_ij_ckeep_(1:n_c_tmp)),2);
AZnV_nix_ = A_full_normalized_n_(:,xdrop_ij_ckeep_(1:n_c_tmp))*tmp_V_;
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_nix_(:,1),AZnV_nix_(:,2),16,mr_ABC__,'filled','MarkerEdgeColor','k');
axisnotick; title('bicluster informed (matlab)'); xlabel('PC1'); ylabel('PC2');
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
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_ni0_driver_(:,1),AZnV_ni0_driver_(:,2),16,mr_ABC__,'filled','MarkerEdgeColor','k');
axisnotick; title('full data-set (driver)'); xlabel('PC1'); ylabel('PC2');
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
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_nix_driver_(:,1),AZnV_nix_driver_(:,2),16,mr_ABC__,'filled','MarkerEdgeColor','k');
axisnotick; title('bicluster-informed (driver)'); xlabel('PC1'); ylabel('PC2');

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we load the principal-component projection from iteration 0, ; ')); end;
if (verbose); disp(sprintf(' %% and use this projection to determine subcontinents. ; ')); end;
if (verbose); disp(sprintf(' %% In this simple example we just threshold the projection. ; ')); end;
if (verbose); disp(sprintf(' %% In more serious applications one could apply k-means clustering ; ')); end;
if (verbose); disp(sprintf(' %% or isosplit to determine the subcontinent memebership. ; ')); end;
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
tmp_threshold = 1;
ij_scA_ = find(AZnV_ni0_driver_(:,1+1)>=tmp_threshold);
ij_scB_ = find(AZnV_ni0_driver_(:,1+1)< tmp_threshold);
mr_scX__ = zeros(n_patient,3);
mr_scX__(ij_scA_,1+0) = 1;
mr_scX__(ij_scB_,1+1) = 1;
subplot(1,2,1);
scatter(AZnV_ni0_driver_(:,1),AZnV_ni0_driver_(:,2),16,mr_ABC__,'filled','MarkerEdgeColor','k');
axisnotick; title('colored by construction'); xlabel('PC1'); ylabel('PC2');
subplot(1,2,2);
scatter(AZnV_ni0_driver_(:,1),AZnV_ni0_driver_(:,2),16,mr_scX__,'filled','MarkerEdgeColor','k');
axisnotick; title('colored by subcontinent'); xlabel('PC1'); ylabel('PC2');

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we set up row masks to pick out the subcontinents. ; ')); end;
%%%%%%%%;
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter.dir_0in,parameter.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter.dir_0in,parameter.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter.dir_0in,parameter.str_prefix);
mr_A_ori_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mr_A_alt_ = mr_A_ori_; %<-- this is the same as the original mask, as it includes all cases. ;
mr_Z_alt_ = mr_Z_ori_; %<-- this is the same as the original mask, as it includes all ctrls. ;
n_bin = 2; %<-- now we set up two ''bins'' (i.e., subcontinents) used to divide the loop-scores. ;
mr_A_alt__ = cell(n_bin,1); %<-- cases for each subcontinent. ;
mr_Z_alt__ = cell(n_bin,1); %<-- ctrls for each subcontinent. ;
mr_scA_ = zeros(n_patient,1); mr_scA_(ij_scA_) = 1;
mr_scB_ = zeros(n_patient,1); mr_scB_(ij_scB_) = 1;
for nbin=0:n_bin-1;
if (nbin==0); %<-- scA. ;
mr_A_alt__{1+nbin} = mr_A_ori_ & mr_scA_;
mr_Z_alt__{1+nbin} = mr_Z_ori_ & mr_scA_;
end;%if (nbin==0); %<-- scA. ;
if (nbin==1); %<-- scB. ;
mr_A_alt__{1+nbin} = mr_A_ori_ & mr_scB_;
mr_Z_alt__{1+nbin} = mr_Z_ori_ & mr_scB_;
end;%if (nbin==1); %<-- scB. ;
end;%for nbin=0:n_bin-1;
%%%%%%%%;
str_mr_0in = 'scX'; %<-- reference for the subcontinent analysis. ;
fname_mr_A_alt_full = sprintf('%s/%s_mr_A_%s_full.b16',parameter.dir_0in,parameter.str_prefix,str_mr_0in);
fname_mr_Z_alt_full = sprintf('%s/%s_mr_Z_%s_full.b16',parameter.dir_0in,parameter.str_prefix,str_mr_0in);
bitj=16;
binary_compress(bitj,mr_A_alt_>0,fname_mr_A_alt_full);
binary_compress(bitj,mr_Z_alt_>0,fname_mr_Z_alt_full);
if (verbose); disp(sprintf(' %% ;')); end;
if (verbose); disp(sprintf(' %% We do the same for each study. ;')); end;
for nbin=0:n_bin-1;
fname_mr_A_alt_ij = sprintf('%s/%s_mr_A_%s_%.2d.b16',parameter.dir_0in,parameter.str_prefix,str_mr_0in,1+nbin);
fname_mr_Z_alt_ij = sprintf('%s/%s_mr_Z_%s_%.2d.b16',parameter.dir_0in,parameter.str_prefix,str_mr_0in,1+nbin);
if (verbose); disp(sprintf(' %% nbin %d/%d <-- ij=1+nbin=%d ...\n %% %s ; \n %% %s ; ',nbin,n_bin,1+nbin,fname_mr_A_alt_ij,fname_mr_Z_alt_ij)); end;
bitj=16;nbin==0;
binary_compress(bitj,mr_A_alt__{1+nbin}>0,fname_mr_A_alt_ij);
binary_compress(bitj,mr_Z_alt__{1+nbin}>0,fname_mr_Z_alt_ij);
end;%for nbin=0:n_bin-1;
%%%%%%%%;

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now for the tedious part: ; ')); end;
if (verbose); disp(sprintf(' %% In order to use the same matlab driver, ; ')); end;
if (verbose); disp(sprintf(' %% we will need to write out parts of the c-input manually. ; ')); end;
if (verbose); disp(sprintf(' %% (I have yet to come up with a cleaner implementation). ; ')); end;
parameter_subcontinents = parameter;
parameter_subcontinents.flag_force_create = 1;
parameter_subcontinents.flag_force_create_fname_0in = 0;
parameter_fname_0in_manual = struct('type','parameter_fname_0in');
if (verbose); disp(sprintf(' %% We will have two bins for the data, ; ')); end;
if (verbose); disp(sprintf(' %% but each bin will draw its values from the same b16 array. ; ')); end;
if (verbose); disp(sprintf(' %% Thus, we need to specify the arrays where the b16 data is located, ; ')); end;
if (verbose); disp(sprintf(' %% once for each bin (i.e., the same name twice). ; ')); end;
if (verbose); disp(sprintf(' %% If we had three bins, we would have to repeat the same name three times. ; ')); end;
dir_0in_plus_prefix = sprintf('%s/%s',parameter.dir_0in,parameter.str_prefix); 
parameter_fname_0in_manual.GLOBAL_A_n_name_ = sprintf('GLOBAL_A_n_name_= %s_A_full_n.b16,%s_A_full_n.b16;\n',dir_0in_plus_prefix,dir_0in_plus_prefix);
parameter_fname_0in_manual.GLOBAL_A_t_name_ = sprintf('GLOBAL_A_t_name_= %s_A_full_t.b16,%s_A_full_t.b16;\n',dir_0in_plus_prefix,dir_0in_plus_prefix);
parameter_fname_0in_manual.GLOBAL_Z_n_name_ = sprintf('GLOBAL_Z_n_name_= %s_A_full_n.b16,%s_A_full_n.b16;\n',dir_0in_plus_prefix,dir_0in_plus_prefix);
parameter_fname_0in_manual.GLOBAL_Z_t_name_ = sprintf('GLOBAL_Z_t_name_= %s_A_full_t.b16,%s_A_full_t.b16;\n',dir_0in_plus_prefix,dir_0in_plus_prefix);
parameter_fname_0in_manual.GLOBAL_T_n_name_ = sprintf('GLOBAL_T_n_name_= %s_T_m2r1_full_n.b16,%s_T_m2r1_full_n.b16;\n',dir_0in_plus_prefix,dir_0in_plus_prefix);
parameter_fname_0in_manual.GLOBAL_T_t_name_ = sprintf('GLOBAL_T_t_name_= %s_T_m2r1_full_t.b16,%s_T_m2r1_full_t.b16;\n',dir_0in_plus_prefix,dir_0in_plus_prefix);
parameter_fname_0in_manual.GLOBAL_S_n_name_ = sprintf('GLOBAL_S_n_name_= %s_T_m2r1_full_n.b16,%s_T_m2r1_full_n.b16;\n',dir_0in_plus_prefix,dir_0in_plus_prefix);
parameter_fname_0in_manual.GLOBAL_S_t_name_ = sprintf('GLOBAL_S_t_name_= %s_T_m2r1_full_t.b16,%s_T_m2r1_full_t.b16;\n',dir_0in_plus_prefix,dir_0in_plus_prefix);
parameter_subcontinents.parameter_fname_0in = parameter_fname_0in_manual;
parameter_subcontinents.n_bin = 2; %<-- two bins. ;
parameter_subcontinents.Ireq = 1; %<-- require only a single bin for significance. ;
parameter_subcontinents.str_mr_0in = str_mr_0in;
if (verbose); disp(sprintf(' %% The other values in the c-input should be filled automatically. ; ')); end;
parameter_subcontinents.nshuffle = 0;
parameter_subcontinents = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_subcontinents);
for nshuffle=1:n_shuffle;
parameter_subcontinents.nshuffle = nshuffle;
parameter_subcontinents.flag_force_create_fname_0in = 0;
parameter_subcontinents.parameter_fname_0in = parameter_fname_0in_manual; %<-- reset fields. ;
[ ...
 parameter_subcontinents ...
] = ...
xxxcluster_fromdisk_uADZSZDA_ver16( ...
 parameter_subcontinents ...
);
end;%for nshuffle=1:n_shuffle;

if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we load the traces from each run of the subcontinent analysis. ;')); end;
if (verbose); disp(sprintf(' %% Each trace records the (reweighted) mean-squared-correlation of the remaining array (across iterations),.')); end;
if (verbose); disp(sprintf(' %% where, roughly speaking, the mean is limited to squared-correlations between members of the same subcontinent. ; ')); end;
trace_subcontinents__ = load_trace__from_dir_ver0(parameter_subcontinents.dir_out_trace);
if (verbose); disp(sprintf(' %% ; ')); end;
if (verbose); disp(sprintf(' %% Now we display the original trace (red) against the label-shuffled traces (black). ;')); end;
if (verbose); disp(sprintf(' %% The horizontal indicates iteration, and the vertical indicates z-score. ;')); end;
if (verbose); disp(sprintf(' %% The maximum negative-log-p-value for the original trace is circled (red).')); end;
if (verbose); disp(sprintf(' %% We compare this result (right) to the previous trace (left). ; ')); end;
if (verbose); disp(sprintf(' %% Note that the results are much more significant. ; ')); end;
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_sml = 0.5;
linewidth_big = 3;
markersize_big = 16;
%%%%%%%%;
subplot(1,2,2);
[tmp_nlpR,ij_nlpR] = max(trace_subcontinents__.nlpR_s0000_);
hold on;
plot(trace_subcontinents__.niter_s0000_,trace_subcontinents__.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(trace_subcontinents__.niter_s0000_,trace_subcontinents__.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(trace_subcontinents__.niter_s0000_(ij_nlpR),trace_subcontinents__.ZR_s0000_(ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([min(trace_subcontinents__.niter_s0000_),max(trace_subcontinents__.niter_s0000_)]); xlabel('iteration');
ylim_ = [-10,+1.125*max(trace_subcontinents__.ZR_s0000_)];
ylim(ylim_);
ylabel('negative-log-p');
title('trace_subcontinents','Interpreter','none');
%%%%%%%%;
subplot(1,2,1);
[tmp_nlpR,ij_nlpR] = max(trace__.nlpR_s0000_);
hold on;
plot(trace__.niter_s0000_,trace__.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(trace__.niter_s0000_,trace__.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(trace__.niter_s0000_(ij_nlpR),trace__.ZR_s0000_(ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([min(trace__.niter_s0000_),max(trace__.niter_s0000_)]); xlabel('iteration');
ylim(ylim_);
ylabel('negative-log-p');
title('trace_original','Interpreter','none');
%%%%%%%%;


