clear;

verbose=1;
if (verbose); disp(sprintf(' %% [entering test_bed_to_b16_flip_ver7]')); end;
rng(0);

%%%%%%%%;
n_snp_cup = 1024*8;
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
%%%%%%%%;
n_patient_cup = 1024*4;
fam_cup_fid_ = cell(n_patient_cup,1);
fam_cup_iid_ = cell(n_patient_cup,1);
fam_cup_yid_ = cell(n_patient_cup,1);
fam_cup_xid_ = cell(n_patient_cup,1);
fam_cup_sex_ = zeros(n_patient_cup,1);
fam_cup_dvx_ = zeros(n_patient_cup,1);
fam_cup_fidandiid_ = cell(n_patient_cup,1);
tmp_sex_ = [1,2,1,2];
tmp_dvx_ = [1,1,2,2];
for npatient_cup=0:n_patient_cup-1;
fam_cup_fid_{1+npatient_cup} = sprintf('fid%.4d',1+npatient_cup);
fam_cup_iid_{1+npatient_cup} = sprintf('iid%.4d',1+npatient_cup);
fam_cup_yid_{1+npatient_cup} = sprintf('yid%.4d',1+npatient_cup);
fam_cup_xid_{1+npatient_cup} = sprintf('xid%.4d',1+npatient_cup);
fam_cup_sex_(1+npatient_cup) = tmp_sex_(1+mod(npatient_cup,4));
fam_cup_dvx_(1+npatient_cup) = tmp_dvx_(1+mod(npatient_cup,4));
end;%for npatient_cup=0:n_patient_cup-1;
for npatient_cup=0:n_patient_cup-1;
fam_cup_fidandiid_{1+npatient_cup} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for npatient_cup=0:n_patient_cup-1;
%%%%%%%%;
frq_s_ = linspace(0,1,n_snp_cup);
ds1_cup_ps__ = zeros(n_patient_cup,n_snp_cup);
ds1_cup_ps__ = rand(n_patient_cup,n_snp_cup)<repmat(frq_s_,[n_patient_cup,1]);
ds2_cup_ps__ = zeros(n_patient_cup,n_snp_cup);
ds2_cup_ps__ = rand(n_patient_cup,n_snp_cup)<repmat(frq_s_,[n_patient_cup,1]);
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

%%%%%%%%;
dir_trunk = sprintf('%s/dir_test_bed_to_b16_flip_ver7',pwd);
if (~exist(dir_trunk,'dir')); disp(sprintf(' %% mkdir %s',dir_trunk)); mkdir(dir_trunk); end;
n_study = 3;
study_trunk_s_ = cell(n_study,1);
study_name_s_ = cell(n_study,1);
for nstudy=0:n_study-1;
study_name = sprintf('study%.2d',nstudy);
study_trunk = sprintf('dir_%s',study_name);
study_name_s_{1+nstudy} = study_name;
study_trunk_s_{1+nstudy} = study_trunk;
tmp_dir = sprintf('%s/%s',dir_trunk,study_trunk);
if (~exist(tmp_dir,'dir')); disp(sprintf(' %% mkdir %s',tmp_dir)); mkdir(tmp_dir); end;
end;%for nstudy=0:n_study-1;
%%%%%%%%;
n_snp_s_ = zeros(n_study,1);
index_nsnp_from_nstudy__ = cell(n_study,1);
for nstudy=0:n_study-1;
n_snp = ceil( 0.5^(1/n_study) * n_snp_cup);
index_nsnp_from_nstudy_ = randperm(n_snp_cup,n_snp) - 1;
index_nsnp_from_nstudy__{1+nstudy} = index_nsnp_from_nstudy_;
n_snp_s_(1+nstudy) = n_snp;
end;%for nstudy=0:n_study-1;
%%%%%%%%;
n_patient_s_ = zeros(n_study,1);
index_npatient_from_nstudy__ = cell(n_study,1);
n_patient_sum=0;
for nstudy=0:n_study-1;
n_patient = floor(n_patient_cup/n_study);
index_npatient_from_nstudy_ = n_patient_sum + [0:n_patient-1];
index_npatient_from_nstudy__{1+nstudy} = index_npatient_from_nstudy_;
n_patient_sum = n_patient_sum + n_patient;
n_patient_s_(1+nstudy) = n_patient;
end;%for nstudy=0:n_study-1;
%%%%%%%%;
fn_fam_s_ = cell(n_study,1);
for nstudy=0:n_study-1;
fn_fam = sprintf('%s/%s/%s.fam',dir_trunk,study_trunk_s_{1+nstudy},study_name_s_{1+nstudy});
fid = fopen(fn_fam,'w');
n_patient = n_patient_s_(1+nstudy);
index_npatient_from_nstudy_ = index_npatient_from_nstudy__{1+nstudy};
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
fn_fam_s_{1+nstudy} = fn_fam;
end;%for nstudy=0:n_study-1;
%%%%%%%%;
fn_bim_s_ = cell(n_study,1);
for nstudy=0:n_study-1;
fn_bim = sprintf('%s/%s/%s.bim',dir_trunk,study_trunk_s_{1+nstudy},study_name_s_{1+nstudy});
fid = fopen(fn_bim,'w');
n_snp = n_snp_s_(1+nstudy);
index_nsnp_from_nstudy_ = index_nsnp_from_nstudy__{1+nstudy};
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
fn_bim_s_{1+nstudy} = fn_bim;
end;%for nstudy=0:n_study-1;
%%%%%%%%;
fn_bed_s_ = cell(n_study,1);
for nstudy=0:n_study-1;
fn_bed = sprintf('%s/%s/%s.bed',dir_trunk,study_trunk_s_{1+nstudy},study_name_s_{1+nstudy});
fid = fopen(fn_bed,'w');
n_snp = n_snp_s_(1+nstudy);
index_nsnp_from_nstudy_ = index_nsnp_from_nstudy__{1+nstudy};
n_patient = n_patient_s_(1+nstudy);
index_npatient_from_nstudy_ = index_npatient_from_nstudy__{1+nstudy};
tmp_dsg_ps__ = dsg_cup_ps__(1+index_npatient_from_nstudy_,1+index_nsnp_from_nstudy_);
[tmp_bed_ps__,n_patient_bed] = dsg_to_bed_0(n_patient,n_snp,tmp_dsg_ps__);
key1 = uint8(108); key2 = uint8(27); key3 = uint8(1);
fwrite(fid,key1,'uint8'); fwrite(fid,key2,'uint8'); fwrite(fid,key3,'uint8');
fwrite(fid,tmp_bed_ps__,'uint8');
fclose(fid);
fn_bed_s_{1+nstudy} = fn_bed;
end;%for nstudy=0:n_study-1;
%%%%%%%%;
n_mds_p = n_patient_cup;
n_mds_v = 2;
mds_pv__ = randn(n_mds_p,n_mds_v);
mds_fidandiid_p_ = cell(n_mds_p,1);
for npatient_cup=0:n_patient_cup-1;
mds_fidandiid_p_{1+npatient_cup} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for npatient_cup=0:n_patient_cup-1;
%%%%%%%%;
n_snp_cap = numel(intersectall(index_nsnp_from_nstudy__));
n_famex = min(128,n_patient_cup/4);
index_npatient_cup_from_nfamex_ = randperm(n_patient_cup,n_famex)-1;
famex_fidandiid_p_ = cell(n_famex,1);
for nfamex=0:n_famex-1;
npatient_cup = index_npatient_cup_from_nfamex_(1+nfamex);
famex_fidandiid_p_{1+nfamex} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for nfamex=0:n_famex-1;
%%%%%%%%;
str_output_prefix = 'test';
parameter = struct('type','parameter');
parameter.verbose = 1;
parameter.str_output_prefix = str_output_prefix;
parameter.dir_trunk = dir_trunk;
parameter.ent_cutoff = 0.001;
parameter.maf_cutoff = 0.010;
[ ...
 parameter ...
] = ...
bed_to_b16_flip_ver7( ...
 parameter ...
,n_study ...
,study_trunk_s_ ...
,study_name_s_ ...
,n_mds_p ...
,n_mds_v ...
,mds_pv__ ...
,mds_fidandiid_p_ ...
,n_famex ...
,famex_fidandiid_p_ ...
);
%%%%%%%%;
str_tmp = sprintf('%s_maf%.2d',str_output_prefix,floor(100*parameter.maf_cutoff));
str_output_prefix_local = sprintf('%s_',str_tmp);
dir_out = sprintf('%s/dir_%s',parameter.dir_trunk,str_tmp);
if (verbose>0); disp(sprintf('%% str_output_prefix_local %s ;\n%% dir_out %s ;\n',str_output_prefix_local,dir_out)); end;


fn_An_full = sprintf('%s/%sA_full_n.b16',dir_out,str_output_prefix_local); An_full__ = binary_uncompress(fn_An_full);
fn_At_full = sprintf('%s/%sA_full_t.b16',dir_out,str_output_prefix_local); At_full__ = binary_uncompress(fn_At_full);
if (verbose>0); disp(sprintf(' %% An_full__ vs At_full__: %0.16f',fnorm(An_full__-transpose(At_full__)))); end;
fn_famext = sprintf('%s/%sfam.ext',dir_out,str_output_prefix_local);
[ ...
 n_patient_ext ...
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
fn_famext ...
);
%%%%%%%%;
fn_bimext = sprintf('%s/%sbim.ext',dir_out,str_output_prefix_local);
[ ...
 n_snp_ext ...
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
fn_bimext ...
);
%%%%%%%%;
[fam_fidandiid_cap_,index_fam_cup_from_cap_,index_fam_ext_from_cap_] = intersect(fam_cup_fidandiid_,famext_fidandiid_,'stable');
index_fam_cup_from_cap_ = index_fam_cup_from_cap_ - 1;
index_fam_ext_from_cap_ = index_fam_ext_from_cap_ - 1;
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
error_sum=0; n_keep=0; n_flip=0;
for nsnp_ext=0:n_snp_ext-1;
nsnp_cup = index_bim_cup_from_ext_(1+nsnp_ext);
assert(strcmp(bim_cup_vid_(1+nsnp_cup),bimext_vid_(1+nsnp_ext)));
tmp_bim_ext_al1 = bimext_al1_(1+nsnp_ext);
tmp_bim_ext_al2 = bimext_al2_(1+nsnp_ext);
tmp_bim_ext_alt = bimext_alt_{1+nsnp_ext};
tmp_bim_cup_al1 = bim_cup_al1_{1+nsnp_cup};
tmp_bim_cup_al2 = bim_cup_al2_{1+nsnp_cup};
tmp_b16_full_ = An_full__(1+index_fam_ext_from_cap_,1+nsnp_ext)>0;
tmp_dsg_cup_ = dsg_cup_ps__(1+index_fam_cup_from_cap_,1+nsnp_cup);
tmp_b16_cup_ = zeros(n_patient_ext,1);
tmp_flip = ( (tmp_bim_ext_al1=='A') | (tmp_bim_ext_al1=='T') );
if tmp_flip==0;
if (~isempty(strfind(tmp_bim_ext_alt,'nor'))); tmp_b16_cup_ = tmp_dsg_cup_==0; end;
if (~isempty(strfind(tmp_bim_ext_alt,'xor'))); tmp_b16_cup_ = tmp_dsg_cup_==1; end;
if (~isempty(strfind(tmp_bim_ext_alt,'and'))); tmp_b16_cup_ = tmp_dsg_cup_==2; end;
n_keep = n_keep + 1;
end;%if tmp_flip==0;
if tmp_flip==1;
if (~isempty(strfind(tmp_bim_ext_alt,'nor'))); tmp_b16_cup_ = tmp_dsg_cup_==2; end;
if (~isempty(strfind(tmp_bim_ext_alt,'xor'))); tmp_b16_cup_ = tmp_dsg_cup_==1; end;
if (~isempty(strfind(tmp_bim_ext_alt,'and'))); tmp_b16_cup_ = tmp_dsg_cup_==0; end;
n_flip = n_flip + 1;
end;%if tmp_flip==1;
tmp_error = fnorm(tmp_b16_cup_ - tmp_b16_full_);
error_sum = error_sum + tmp_error;
end;%for nsnp_ext=0:n_snp_ext-1;
disp(sprintf(' %% n_keep %0.4d n_flip %0.4d error_sum: %0.16f',n_keep,n_flip,error_sum));
disp(sprintf(' %% max ent: %0.6f vs ent_cutoff %0.6f',max(bimext_ent_),parameter.ent_cutoff));
disp(sprintf(' %% min maf: %0.6f vs maf_cutoff %0.6f',min(bimext_maf_),parameter.maf_cutoff));

if (verbose); disp(sprintf(' %% [finished test_bed_to_b16_flip_ver7]')); end;







