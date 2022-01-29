function ...
test_stripped_xxxcluster_fromdisk_ver16( ...
 global_parameter ...
);

na=0;
if (nargin<1+na); global_parameter=[]; end; na=na+1;
if isempty(global_parameter); global_parameter=struct('type','parameter'); end;
if ~isfield(global_parameter,'flag_verbose'); global_parameter.flag_verbose = 1; end;
flag_verbose = global_parameter.flag_verbose;
if ~isfield(global_parameter,'flag_force_create'); global_parameter.flag_force_create = 0; end;
flag_force_create = global_parameter.flag_force_create;

flag_disp = flag_verbose;
if (flag_verbose); disp(sprintf(' %% [entering test_stripped_xxxcluster_fromdisk_ver16]')); end;
setup_0;
rng(0); nf=0;

if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% This runs a stripped-down version of test_xxxcluster_fromdisk_ver16. ;')); end;
if (flag_verbose); disp(sprintf(' %% The main difference is that the data is never loaded as doubles. ;')); end;
if (flag_verbose); disp(sprintf(' %% This saves quite a bit of memory. ;')); end;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Consequently, you should be able to insert your own bed, bim and fam-files into this pipeline. ;')); end;

if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% First we generate a very simple example of dosage-data ;')); end;
if (flag_verbose); disp(sprintf(' %% which exhibits the following biclusters: ;')); end;
if (flag_verbose); disp(sprintf(' %% bcA: a case-specific bicluster which is balanced in mds-space (which we want to find) ;')); end;
if (flag_verbose); disp(sprintf(' %% bcB: a non-specific bicluster spanning the cases and controls (which we want to ignore) ;')); end;
if (flag_verbose); disp(sprintf(' %% bcC: a case-specific bicluster which is localized in mds-space (which we want to ignore) ;')); end;
str_lak_vs_dex = 'dex';
r_fraction_bcA = 0.125;
r_fraction_bcB = 0.250;
r_fraction_bcC = 0.125;
c_fraction_bcA = 0.0625;
c_fraction_bcB = 0.1250;
c_fraction_bcC = 0.0625;
maf_lob = 0.45;
maf_upb = 0.55;

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
cij_ = transpose(1:n_snp_cup);
cij_tmp_ = cij_;
cij_bcA_ = cij_tmp_(randperm(n_snp_cup,floor(n_snp_cup*c_fraction_bcA)));
n_snp_bcA = numel(cij_bcA_);
%cij_tmp_ = setdiff(cij_,cij_bcA_);
cij_tmp_ = cij_;
cij_bcB_ = cij_tmp_(randperm(numel(cij_tmp_),floor(n_snp_cup*c_fraction_bcB)));
n_snp_bcB = numel(cij_bcB_);
%cij_tmp_ = setdiff(cij_,union(cij_bcA_,cij_bcB_));
cij_tmp_ = cij_;
cij_bcC_ = cij_tmp_(randperm(numel(cij_tmp_),floor(n_snp_cup*c_fraction_bcC)));
n_snp_bcC = numel(cij_bcC_);
%%%%%%%%;
n_patient_cup = 1024*5;
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
rij_case_bcA_ = rij_case_tmp_(randperm(numel(rij_case_tmp_),floor(r_fraction_bcA*n_patient_case)));
rij_ctrl_bcA_ = [];
n_patient_case_bcA = numel(rij_case_bcA_);
n_patient_ctrl_bcA = numel(rij_ctrl_bcA_);
%rij_case_tmp_ = setdiff(rij_case_,rij_case_bcA_);
rij_case_tmp_ = rij_case_;
%rij_ctrl_tmp_ = setdiff(rij_ctrl_,rij_ctrl_bcA_);
rij_ctrl_tmp_ = rij_ctrl_;
rij_case_bcB_ = rij_case_tmp_(randperm(numel(rij_case_tmp_),floor(r_fraction_bcB*n_patient_case)));
rij_ctrl_bcB_ = rij_ctrl_tmp_(randperm(numel(rij_ctrl_tmp_),floor(r_fraction_bcB*n_patient_ctrl)));
n_patient_case_bcB = numel(rij_case_bcB_);
n_patient_ctrl_bcB = numel(rij_ctrl_bcB_);
%rij_case_tmp_ = setdiff(rij_case_,union(rij_case_bcA_,rij_case_bcB_));
rij_case_tmp_ = rij_case_;
ij_tmp_ = find( (mds_w_p_(rij_case_tmp_) > 0*pi/6) & (mds_w_p_(rij_case_tmp_) < 3*pi/6) & (mds_r_p_(rij_case_tmp_) > 0.05 ) );
n_patient_case_bcC = min(floor(r_fraction_bcC*n_patient_case),numel(ij_tmp_));
rij_case_bcC_ = rij_case_tmp_(ij_tmp_(1:n_patient_case_bcC));
rij_ctrl_bcC_ = [];
n_patient_case_bcC = numel(rij_case_bcC_);
n_patient_ctrl_bcC = numel(rij_ctrl_bcC_);
rij_case_bcA_cmp_ = setdiff(rij_case_,rij_case_bcA_); n_patient_case_bcA_cmp = numel(rij_case_bcA_cmp_);
rij_case_bcB_cmp_ = setdiff(rij_case_,rij_case_bcB_); n_patient_case_bcB_cmp = numel(rij_case_bcB_cmp_);
rij_ctrl_bcB_cmp_ = setdiff(rij_ctrl_,rij_ctrl_bcB_); n_patient_ctrl_bcB_cmp = numel(rij_ctrl_bcB_cmp_);
rij_case_bcC_cmp_ = setdiff(rij_case_,rij_case_bcC_); n_patient_case_bcC_cmp = numel(rij_case_bcC_cmp_);
rij_ctrl_bcC_cmp_ = setdiff(rij_ctrl_,rij_ctrl_bcC_); n_patient_ctrl_bcC_cmp = numel(rij_ctrl_bcC_cmp_);
%%%%%%%%;
frq_s_ = linspace(maf_lob,maf_upb,n_snp_cup);
ds1_cup_ps__ = zeros(n_patient_cup,n_snp_cup);
ds1_cup_ps__ = rand(n_patient_cup,n_snp_cup)<repmat(frq_s_,[n_patient_cup,1]);
ds2_cup_ps__ = zeros(n_patient_cup,n_snp_cup);
ds2_cup_ps__ = rand(n_patient_cup,n_snp_cup)<repmat(frq_s_,[n_patient_cup,1]);
if (strcmp(str_lak_vs_dex,'lak'));
%%%%;
% now correlate diplotypes within bcA. ;
%%%%;
ds1_cup_ps__(rij_case_bcA_,cij_bcA_) = repmat(transpose(linspace(0,1,n_patient_case_bcA)),[1,n_snp_bcA]) < repmat(frq_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds1_cup_ps__(rij_case_bcA_,cij_bcA_) = ds1_cup_ps__(rij_case_bcA_(randperm(n_patient_case_bcA)),cij_bcA_);
ds2_cup_ps__(rij_case_bcA_,cij_bcA_) = repmat(transpose(linspace(0,1,n_patient_case_bcA)),[1,n_snp_bcA]) < repmat(frq_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds2_cup_ps__(rij_case_bcA_,cij_bcA_) = ds2_cup_ps__(rij_case_bcA_(randperm(n_patient_case_bcA)),cij_bcA_);
%%%%;
% now correlate diplotypes within bcB. ;
%%%%;
ds1_cup_ps__(rij_case_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_case_bcB)),[1,n_snp_bcB]) < repmat(frq_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds1_cup_ps__(rij_case_bcB_,cij_bcB_) = ds1_cup_ps__(rij_case_bcB_(randperm(n_patient_case_bcB)),cij_bcB_);
ds2_cup_ps__(rij_case_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_case_bcB)),[1,n_snp_bcB]) < repmat(frq_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds2_cup_ps__(rij_case_bcB_,cij_bcB_) = ds2_cup_ps__(rij_case_bcB_(randperm(n_patient_case_bcB)),cij_bcB_);
ds1_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_ctrl_bcB)),[1,n_snp_bcB]) < repmat(frq_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds1_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = ds1_cup_ps__(rij_ctrl_bcB_(randperm(n_patient_ctrl_bcB)),cij_bcB_);
ds2_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_ctrl_bcB)),[1,n_snp_bcB]) < repmat(frq_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds2_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = ds2_cup_ps__(rij_ctrl_bcB_(randperm(n_patient_ctrl_bcB)),cij_bcB_);
%%%%;
% now correlate diplotypes within bcC. ;
%%%%;
ds1_cup_ps__(rij_case_bcC_,cij_bcC_) = repmat(transpose(linspace(0,1,n_patient_case_bcC)),[1,n_snp_bcC]) < repmat(frq_s_(cij_bcC_),[n_patient_case_bcC,1]);
ds1_cup_ps__(rij_case_bcC_,cij_bcC_) = ds1_cup_ps__(rij_case_bcC_(randperm(n_patient_case_bcC)),cij_bcC_);
ds2_cup_ps__(rij_case_bcC_,cij_bcC_) = repmat(transpose(linspace(0,1,n_patient_case_bcC)),[1,n_snp_bcC]) < repmat(frq_s_(cij_bcC_),[n_patient_case_bcC,1]);
ds2_cup_ps__(rij_case_bcC_,cij_bcC_) = ds2_cup_ps__(rij_case_bcC_(randperm(n_patient_case_bcC)),cij_bcC_);
end;%if (strcmp(str_lak_vs_dex,'lak'));
%%%%%%%%;
if (strcmp(str_lak_vs_dex,'dex'));
%%%%;
% now lift diplotypes within bcA. ;
%%%%;
maf_upb_use = 1.0;
maf_lob_use = 0.0;
frq_case_bcA_s_ = (frq_s_ >= 0.5) .* (frq_s_ + 0.5*(maf_upb_use-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + 0.5*(maf_lob_use-frq_s_)) ;
frq_case_bcA_cmp_s_ = (frq_s_*n_patient_case - frq_case_bcA_s_*n_patient_case_bcA)/(n_patient_case - n_patient_case_bcA);
%ds1_cup_ps__(rij_case_bcA_,cij_bcA_) = repmat(transpose(linspace(0,1,n_patient_case_bcA)),[1,n_snp_bcA]) < repmat(frq_case_bcA_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds1_cup_ps__(rij_case_bcA_,cij_bcA_) = rand(n_patient_case_bcA,n_snp_bcA) < repmat(frq_case_bcA_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds1_cup_ps__(rij_case_bcA_,cij_bcA_) = ds1_cup_ps__(rij_case_bcA_(randperm(n_patient_case_bcA)),cij_bcA_);
%ds2_cup_ps__(rij_case_bcA_,cij_bcA_) = repmat(transpose(linspace(0,1,n_patient_case_bcA)),[1,n_snp_bcA]) < repmat(frq_case_bcA_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds2_cup_ps__(rij_case_bcA_,cij_bcA_) = rand(n_patient_case_bcA,n_snp_bcA) < repmat(frq_case_bcA_s_(cij_bcA_),[n_patient_case_bcA,1]);
ds2_cup_ps__(rij_case_bcA_,cij_bcA_) = ds2_cup_ps__(rij_case_bcA_(randperm(n_patient_case_bcA)),cij_bcA_);
ds1_cup_ps__(rij_case_bcA_cmp_,cij_bcA_) = rand(n_patient_case_bcA_cmp,n_snp_bcA) < repmat(frq_case_bcA_cmp_s_(cij_bcA_),[n_patient_case_bcA_cmp,1]);
ds1_cup_ps__(rij_case_bcA_cmp_,cij_bcA_) = ds1_cup_ps__(rij_case_bcA_cmp_(randperm(n_patient_case_bcA_cmp)),cij_bcA_);
ds2_cup_ps__(rij_case_bcA_cmp_,cij_bcA_) = rand(n_patient_case_bcA_cmp,n_snp_bcA) < repmat(frq_case_bcA_cmp_s_(cij_bcA_),[n_patient_case_bcA_cmp,1]);
ds2_cup_ps__(rij_case_bcA_cmp_,cij_bcA_) = ds2_cup_ps__(rij_case_bcA_cmp_(randperm(n_patient_case_bcA_cmp)),cij_bcA_);
%%%%;
% now lift diplotypes within bcB. ;
%%%%;
frq_case_bcB_s_ = (frq_s_ >= 0.5) .* (frq_s_ + 0.5*(maf_upb_use-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + 0.5*(maf_lob_use-frq_s_)) ;
frq_case_bcB_cmp_s_ = (frq_s_*n_patient_case - frq_case_bcB_s_*n_patient_case_bcB)/(n_patient_case - n_patient_case_bcB);
%ds1_cup_ps__(rij_case_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_case_bcB)),[1,n_snp_bcB]) < repmat(frq_case_bcB_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds1_cup_ps__(rij_case_bcB_,cij_bcB_) = rand(n_patient_case_bcB,n_snp_bcB) < repmat(frq_case_bcB_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds1_cup_ps__(rij_case_bcB_,cij_bcB_) = ds1_cup_ps__(rij_case_bcB_(randperm(n_patient_case_bcB)),cij_bcB_);
%ds2_cup_ps__(rij_case_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_case_bcB)),[1,n_snp_bcB]) < repmat(frq_case_bcB_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds2_cup_ps__(rij_case_bcB_,cij_bcB_) = rand(n_patient_case_bcB,n_snp_bcB) < repmat(frq_case_bcB_s_(cij_bcB_),[n_patient_case_bcB,1]);
ds2_cup_ps__(rij_case_bcB_,cij_bcB_) = ds2_cup_ps__(rij_case_bcB_(randperm(n_patient_case_bcB)),cij_bcB_);
ds1_cup_ps__(rij_case_bcB_cmp_,cij_bcB_) = rand(n_patient_case_bcB_cmp,n_snp_bcB) < repmat(frq_case_bcB_cmp_s_(cij_bcB_),[n_patient_case_bcB_cmp,1]);
ds1_cup_ps__(rij_case_bcB_cmp_,cij_bcB_) = ds1_cup_ps__(rij_case_bcB_cmp_(randperm(n_patient_case_bcB_cmp)),cij_bcB_);
ds2_cup_ps__(rij_case_bcB_cmp_,cij_bcB_) = rand(n_patient_case_bcB_cmp,n_snp_bcB) < repmat(frq_case_bcB_cmp_s_(cij_bcB_),[n_patient_case_bcB_cmp,1]);
ds2_cup_ps__(rij_case_bcB_cmp_,cij_bcB_) = ds2_cup_ps__(rij_case_bcB_cmp_(randperm(n_patient_case_bcB_cmp)),cij_bcB_);
frq_ctrl_bcB_s_ = (frq_s_ >= 0.5) .* (frq_s_ + 0.5*(maf_upb_use-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + 0.5*(maf_lob_use-frq_s_)) ;
frq_ctrl_bcB_cmp_s_ = (frq_s_*n_patient_ctrl - frq_ctrl_bcB_s_*n_patient_ctrl_bcB)/(n_patient_ctrl - n_patient_ctrl_bcB);
%ds1_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_ctrl_bcB)),[1,n_snp_bcB]) < repmat(frq_ctrl_bcB_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds1_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = rand(n_patient_ctrl_bcB,n_snp_bcB) < repmat(frq_ctrl_bcB_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds1_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = ds1_cup_ps__(rij_ctrl_bcB_(randperm(n_patient_ctrl_bcB)),cij_bcB_);
%ds2_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = repmat(transpose(linspace(0,1,n_patient_ctrl_bcB)),[1,n_snp_bcB]) < repmat(frq_ctrl_bcB_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds2_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = rand(n_patient_ctrl_bcB,n_snp_bcB) < repmat(frq_ctrl_bcB_s_(cij_bcB_),[n_patient_ctrl_bcB,1]);
ds2_cup_ps__(rij_ctrl_bcB_,cij_bcB_) = ds2_cup_ps__(rij_ctrl_bcB_(randperm(n_patient_ctrl_bcB)),cij_bcB_);
ds1_cup_ps__(rij_ctrl_bcB_cmp_,cij_bcB_) = rand(n_patient_ctrl_bcB_cmp,n_snp_bcB) < repmat(frq_ctrl_bcB_cmp_s_(cij_bcB_),[n_patient_ctrl_bcB_cmp,1]);
ds1_cup_ps__(rij_ctrl_bcB_cmp_,cij_bcB_) = ds1_cup_ps__(rij_ctrl_bcB_cmp_(randperm(n_patient_ctrl_bcB_cmp)),cij_bcB_);
ds2_cup_ps__(rij_ctrl_bcB_cmp_,cij_bcB_) = rand(n_patient_ctrl_bcB_cmp,n_snp_bcB) < repmat(frq_ctrl_bcB_cmp_s_(cij_bcB_),[n_patient_ctrl_bcB_cmp,1]);
ds2_cup_ps__(rij_ctrl_bcB_cmp_,cij_bcB_) = ds2_cup_ps__(rij_ctrl_bcB_cmp_(randperm(n_patient_ctrl_bcB_cmp)),cij_bcB_);
%%%%;
% now lift diplotypes within bcC. ;
%%%%;
maf_upb_use = 1.0;
maf_lob_use = 0.0;
frq_case_bcC_s_ = (frq_s_ >= 0.5) .* (frq_s_ + 0.5*(maf_upb_use-frq_s_)) + (frq_s_ < 0.5) .* (frq_s_ + 0.5*(maf_lob_use-frq_s_)) ;
frq_case_bcC_cmp_s_ = (frq_s_*n_patient_case - frq_case_bcC_s_*n_patient_case_bcC)/(n_patient_case - n_patient_case_bcC);
%ds1_cup_ps__(rij_case_bcC_,cij_bcC_) = repmat(transpose(linspace(0,1,n_patient_case_bcC)),[1,n_snp_bcC]) < repmat(frq_case_bcC_s_(cij_bcC_),[n_patient_case_bcC,1]);
ds1_cup_ps__(rij_case_bcC_,cij_bcC_) = rand(n_patient_case_bcC,n_snp_bcC) < repmat(frq_case_bcC_s_(cij_bcC_),[n_patient_case_bcC,1]);
ds1_cup_ps__(rij_case_bcC_,cij_bcC_) = ds1_cup_ps__(rij_case_bcC_(randperm(n_patient_case_bcC)),cij_bcC_);
%ds2_cup_ps__(rij_case_bcC_,cij_bcC_) = repmat(transpose(linspace(0,1,n_patient_case_bcC)),[1,n_snp_bcC]) < repmat(frq_case_bcC_s_(cij_bcC_),[n_patient_case_bcC,1]);
ds2_cup_ps__(rij_case_bcC_,cij_bcC_) = rand(n_patient_case_bcC,n_snp_bcC) < repmat(frq_case_bcC_s_(cij_bcC_),[n_patient_case_bcC,1]);
ds2_cup_ps__(rij_case_bcC_,cij_bcC_) = ds2_cup_ps__(rij_case_bcC_(randperm(n_patient_case_bcC)),cij_bcC_);
ds1_cup_ps__(rij_case_bcC_cmp_,cij_bcC_) = rand(n_patient_case_bcC_cmp,n_snp_bcC) < repmat(frq_case_bcC_cmp_s_(cij_bcC_),[n_patient_case_bcC_cmp,1]);
ds1_cup_ps__(rij_case_bcC_cmp_,cij_bcC_) = ds1_cup_ps__(rij_case_bcC_cmp_(randperm(n_patient_case_bcC_cmp)),cij_bcC_);
ds2_cup_ps__(rij_case_bcC_cmp_,cij_bcC_) = rand(n_patient_case_bcC_cmp,n_snp_bcC) < repmat(frq_case_bcC_cmp_s_(cij_bcC_),[n_patient_case_bcC_cmp,1]);
ds2_cup_ps__(rij_case_bcC_cmp_,cij_bcC_) = ds2_cup_ps__(rij_case_bcC_cmp_(randperm(n_patient_case_bcC_cmp)),cij_bcC_);
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
; setdiff(rij_case_bcA_,unionall({rij_case_bcB_,rij_case_bcC_})) ...
; setdiff(rij_case_bcB_,unionall({rij_case_bcC_})) ...
; setdiff(rij_case_bcC_,unionall({})) ...
; setdiff(rij_case_,unionall({rij_case_bcA_,rij_case_bcB_,rij_case_bcC_})) ...
];
rij_ctrl_prm_ = [ ... 
; setdiff(rij_ctrl_bcA_,unionall({rij_ctrl_bcB_,rij_ctrl_bcC_})) ...
; setdiff(rij_ctrl_bcB_,unionall({rij_ctrl_bcC_})) ...
; setdiff(rij_ctrl_bcC_,unionall({})) ...
; setdiff(rij_ctrl_,unionall({rij_ctrl_bcA_,rij_ctrl_bcB_,rij_ctrl_bcC_})) ...
];
rij_prm_ = [ rij_case_prm_ ; rij_ctrl_prm_ ];
cij_prm_ = [ ...
; setdiff(cij_bcA_,unionall({cij_bcB_,cij_bcC_})) ...
; setdiff(cij_bcB_,unionall({cij_bcC_})) ...
; setdiff(cij_bcC_,unionall({})) ...
; setdiff(cij_,unionall({cij_bcA_,cij_bcB_,cij_bcC_})) ...
];
%%%%%%%%;
mr_bcA_ = zeros(n_patient_cup,1); mr_bcA_(union(rij_case_bcA_,rij_ctrl_bcA_)) = 1;
mr_bcB_ = zeros(n_patient_cup,1); mr_bcB_(union(rij_case_bcB_,rij_ctrl_bcB_)) = 1;
mr_bcC_ = zeros(n_patient_cup,1); mr_bcC_(union(rij_case_bcC_,rij_ctrl_bcC_)) = 1;
mc_cup_bcA_ = zeros(n_snp_cup,1); mc_cup_bcA_(cij_bcA_) = 1;
mc_cup_bcB_ = zeros(n_snp_cup,1); mc_cup_bcB_(cij_bcB_) = 1;
mc_cup_bcC_ = zeros(n_snp_cup,1); mc_cup_bcC_(cij_bcC_) = 1;
%%%%%%%;

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we visualize the dosage-data (magenta and cyan) for the cases (top) and ctrls (bottom). ;')); end;
if (flag_verbose); disp(sprintf(' %% The first two mds-components are shown on the far left column (green and yellow). ;')); end;
if (flag_verbose); disp(sprintf(' %% Adjacent are the row-masks associated with biclusters bcA, bcB and bcC (grey). ;')); end;
if (flag_verbose); disp(sprintf(' %% Atop lie the column-masks associated with biclusters bcA, bcB and bcC (grey). ;')); end;
if (flag_verbose); disp(sprintf(' %% Note that bcA and bcB are both ''global'' (i.e., mds-balanced), ;')); end;
if (flag_verbose); disp(sprintf(' %% while bcC is ''local'' (i.e., biased with respect to the mds-components). ;')); end;
if (flag_verbose); disp(sprintf(' %% Note also that bcA and bcC are both ''case-specific'', ;')); end;
if (flag_verbose); disp(sprintf(' %% while bcB is ''non-specific'' (i.e., extends across cases and ctrls). ;')); end;
if flag_disp;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,512,1024]); 
imagesc_uAZm_1( ...
 mds_pv__(rij_prm_,:) ...
,dsg_cup_ps__(rij_case_prm_,cij_prm_) ...
,dsg_cup_ps__(rij_ctrl_prm_,cij_prm_) ...
,[ mr_bcA_(rij_prm_) , mr_bcB_(rij_prm_) , mr_bcC_(rij_prm_) ] ...
,[ mc_cup_bcA_(cij_prm_) , mc_cup_bcB_(cij_prm_) , mc_cup_bcC_(cij_prm_) ] ...
);
end;%if flag_disp;

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we convert the dosage-data to b16-data. ;')); end;
if (flag_verbose); disp(sprintf(' %% This encodes each dosage value as a 3-bit binary array. ;')); end;
if (flag_verbose); disp(sprintf(' %% dosage ''0'' is encoded as [ 1 0 0 ]. ;')); end;
if (flag_verbose); disp(sprintf(' %% dosage ''1'' is encoded as [ 0 1 0 ]. ;')); end;
if (flag_verbose); disp(sprintf(' %% dosage ''2'' is encoded as [ 0 0 1 ]. ;')); end;
if (flag_verbose); disp(sprintf(' %% For simplicity, we put everything into a single study directory. ;')); end;
if (flag_verbose); disp(sprintf(' %% (although the dosage-data could instead be split across multiple files). ;')); end;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% This is where you can put in your own bed, bim and fam-files. ;')); end;
if (flag_verbose); disp(sprintf(' %% To do so you need to change the paths to these files below. ;')); end;
if (flag_verbose); disp(sprintf(' %% This involves changing: ;')); end;
if (flag_verbose); disp(sprintf(' %% study_trunk: the path to the directory where the study files are contained. ;')); end;
if (flag_verbose); disp(sprintf(' %% study_name: the prefix for the study files. ;')); end;
if (flag_verbose); disp(sprintf(' %%\t Assumes that the bed, bim and fam-files ;')); end;
if (flag_verbose); disp(sprintf(' %%\t  each start with the study_name as a prefix. ;')); end;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% If you do indeed insert your own bed, bim and fam-files, ;')); end;
if (flag_verbose); disp(sprintf(' %% you likely have your own mds-components you wish to correct for. ;')); end;
if (flag_verbose); disp(sprintf(' %% To input your own mds-components you need to change both: ;')); end;
if (flag_verbose); disp(sprintf(' %%\t n_mds_p: integer number of patients in mds_pv__ ;')); end;
if (flag_verbose); disp(sprintf(' %%\t n_mds_v: integer number of variables in mds_pv__ ;')); end;
if (flag_verbose); disp(sprintf(' %%\t mds_pv__: a double-array of size (n_mds_p,n_mds_v) storing the mds-components ;')); end;
if (flag_verbose); disp(sprintf(' %%\t mds_fidandiid_p_: a cell array of size n_mds_p storing strings which are ;')); end;
if (flag_verbose); disp(sprintf(' %%\t\t the concatenation of fid, ''&'', and iid for each of the n_mds_p patients. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note that n_mds_v can be zero, and mds_pv__ can be an array of size (n_mds_p,0). ;')); end;
%%%%%%%%;
dir_trunk = sprintf('%s/dir_test_xxxcluster_fromdisk_ver16',pwd);
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
fn_mds = sprintf('%s/%s/%s_mds_csv.txt',dir_trunk,study_trunk,study_name);
if (flag_force_create | ~exist(fn_mds,'file'));
writematrix(mds_pv__,fn_mds,'filetype','text','delimiter',',');
end;%if (flag_force_create | ~exist(fn_mds,'file'));
fn_mds = sprintf('%s/%s/%s_mds_tsv.txt',dir_trunk,study_trunk,study_name);
if (flag_force_create | ~exist(fn_mds,'file'));
writematrix(mds_pv__,fn_mds,'filetype','text','delimiter','\t');
end;%if (flag_force_create | ~exist(fn_mds,'file'));
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
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% This is where you would specify the identity of the patients listed in mds_pv__. ;')); end;
mds_fidandiid_p_ = cell(n_mds_p,1);
for npatient_cup=0:n_patient_cup-1;
mds_fidandiid_p_{1+npatient_cup} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for npatient_cup=0:n_patient_cup-1;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% You can also set certain people to be excluded (similar to a famex file). ;')); end;
if (flag_verbose); disp(sprintf(' %% For now we skip this step. ;')); end;
n_famex = 0;
index_npatient_cup_from_nfamex_ = randperm(n_patient_cup,n_famex)-1;
famex_fidandiid_p_ = cell(n_famex,1);
for nfamex=0:n_famex-1;
npatient_cup = index_npatient_cup_from_nfamex_(1+nfamex);
famex_fidandiid_p_{1+nfamex} = sprintf('%s&%s',fam_cup_fid_{1+npatient_cup},fam_cup_iid_{1+npatient_cup});
end;%for nfamex=0:n_famex-1;
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
str_output_prefix = 'test';
parameter = struct('type','parameter');
parameter.flag_verbose = max(0,flag_verbose-1);
parameter.str_output_prefix = str_output_prefix;
parameter.dir_trunk = dir_trunk;
parameter.ent_cutoff = 0.001; %<-- this is a typical entropy-cutoff. ;
parameter.ent_cutoff = prctile(dsg_cup_I_opt_,99); %<-- for this example we will make sure to include most of the snps. ;
if (flag_verbose); disp(sprintf(' %% Note: using parameter.ent_cutoff: %0.16f',parameter.ent_cutoff)); end;
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
if (flag_verbose>0); disp(sprintf(' %% str_output_prefix_local %s ;\n %% dir_out %s ;\n',str_output_prefix_local,dir_out)); end;
%%%%%%%%;

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we do not uncompress the b16 file. ;')); end;
if (flag_verbose); disp(sprintf(' %% This step is too slow for larger data-sets. ;')); end;
if (flag_verbose); disp(sprintf(' %% We can still load the fam.ext and bim.ext files. ;')); end;
fname_b16 = sprintf('%s/%s_A_full_n.b16',dir_out,str_output_prefix_plus_maf);
%A_full_n_ = binary_uncompress(fname_b16);
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
if (flag_verbose); disp(sprintf(' %% famext_fid_ vs fam_cup_fid_: %0.16f',fnorm(ij_fam_ext_from_cap_-transpose(1:n_patient)))); end;
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
mc_cup_bcA_ = zeros(n_snp_cup,1);
mc_cup_bcA_(cij_bcA_) = 1;
mc_cup_bcA_from_ext_ = transpose(mc_cup_bcA_) * ij_bim_cup_from_ext__;
n_snp_ext_bcA = sum(mc_cup_bcA_from_ext_);
[~,ij_mc_cup_bcA_from_ext_] = sort(mc_cup_bcA_from_ext_,'ascend');
mc_ext_bcA_ = zeros(n_snp_bimext,1);
mc_ext_bcA_(ij_mc_cup_bcA_from_ext_(end-n_snp_ext_bcA+1:end)) = 1;
%%%%;
mc_cup_bcB_ = zeros(n_snp_cup,1);
mc_cup_bcB_(cij_bcB_) = 1;
mc_cup_bcB_from_ext_ = transpose(mc_cup_bcB_) * ij_bim_cup_from_ext__;
n_snp_ext_bcB = sum(mc_cup_bcB_from_ext_);
[~,ij_mc_cup_bcB_from_ext_] = sort(mc_cup_bcB_from_ext_,'ascend');
mc_ext_bcB_ = zeros(n_snp_bimext,1);
mc_ext_bcB_(ij_mc_cup_bcB_from_ext_(end-n_snp_ext_bcB+1:end)) = 1;
%%%%;
mc_cup_bcC_ = zeros(n_snp_cup,1);
mc_cup_bcC_(cij_bcC_) = 1;
mc_cup_bcC_from_ext_ = transpose(mc_cup_bcC_) * ij_bim_cup_from_ext__;
n_snp_ext_bcC = sum(mc_cup_bcC_from_ext_);
[~,ij_mc_cup_bcC_from_ext_] = sort(mc_cup_bcC_from_ext_,'ascend');
mc_ext_bcC_ = zeros(n_snp_bimext,1);
mc_ext_bcC_(ij_mc_cup_bcC_from_ext_(end-n_snp_ext_bcC+1:end)) = 1;
%%%%;
cij_ex2_ = ...
[ ...
 find(mc_ext_bcA_) ...
; find(mc_ext_bcB_) ...
; find(mc_ext_bcC_) ...
; setdiff(transpose(1:n_snp_bimext),find(mc_ext_bcA_ + mc_ext_bcB_ + mc_ext_bcC_)) ...
];
%%%%;
mr_bcx_ = zeros(n_patient,1);
mr_bcx_(rij_case_bcC_) = 1;
mr_bcx_(rij_case_bcB_) = 2;
mr_bcx_(rij_ctrl_bcB_) = 2;
mr_bcx_(rij_case_bcA_) = 3;
%%%%;
mc_bcx_ = zeros(n_snp_bimext,1);
mc_bcx_(find(mc_ext_bcA_))=3;
mc_bcx_(find(mc_ext_bcB_))=2;
mc_bcx_(find(mc_ext_bcC_))=1;
n_snp_bcx = numel(find(mc_bcx_));
n_snp_ext_bcA = numel(find(mc_bcx_==3));
n_snp_ext_bcB = numel(find(mc_bcx_==2));
n_snp_ext_bcC = numel(find(mc_bcx_==1));
%%%%;
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% We do not visualize the b16-data. ;')); end;
%%%%;
%{
if flag_disp;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,1024]); fig80s;
imagesc_uAZm_1( ...
 mds_pv__(rij_prm_,:) ...
,A_full_n_(rij_case_prm_,cij_ex2_) ...
,A_full_n_(rij_ctrl_prm_,cij_ex2_) ...
,[ mr_bcA_(rij_prm_) , mr_bcB_(rij_prm_) , mr_bcC_(rij_prm_) ] ...
,[ mc_ext_bcA_(cij_ex2_) , mc_ext_bcB_(cij_ex2_) , mc_ext_bcC_(cij_ex2_) ] ...
);
 end;%if flag_disp;
 %}
%%%%;

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we run the basic biclustering method on the b16-data. ;')); end;
if (flag_verbose); disp(sprintf(' %% This should take approximately 10 seconds for differentially-expressed biclustering. ;')); end;
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
parameter.flag_verbose = max(0,flag_verbose-1);
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

%{
if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we display the results of this single biclustering run. ;')); end;
if (flag_verbose); disp(sprintf(' %% The rows and cols are ordered by the biclustering-algorithm. ;')); end;
if (flag_verbose); disp(sprintf(' %% Those rows and cols eliminated first are shown towards the bottom-right. ;')); end;
if (flag_verbose); disp(sprintf(' %% Those rows and cols retained longest are shown towards the top-left. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note that most of bcA is retained until the end. ;')); end;
if (flag_verbose); disp(sprintf(' %% That is to say, the biclustering-algorithm successfully focused on bcA. ;')); end;
if flag_disp;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,1024]); fig80s;
tmp_rij_prm_ = [xdrop_.ij_rkeep_;rij_ctrl_prm_];
tmp_cij_prm_ = [xdrop_.ij_ckeep_];
imagesc_uAZm_1( ...
 mds_pv__(tmp_rij_prm_,:) ...
,A_full_n_(xdrop_.ij_rkeep_,xdrop_.ij_ckeep_) ...
,A_full_n_(rij_ctrl_prm_,xdrop_.ij_ckeep_) ...
,[ mr_bcA_(tmp_rij_prm_) , mr_bcB_(tmp_rij_prm_) , mr_bcC_(tmp_rij_prm_) ] ...
,[ mc_ext_bcA_(tmp_cij_prm_) , mc_ext_bcB_(tmp_cij_prm_) , mc_ext_bcC_(tmp_cij_prm_) ] ...
);
 end;%if flag_disp;
 %}

if (flag_verbose); disp(sprintf(' %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now we run a permutation-test. ;')); end;
if (flag_verbose); disp(sprintf(' %% This involves rerunning the biclustering-algorithm on ;')); end;
if (flag_verbose); disp(sprintf(' %% several label-permuted data-sets. ;')); end;
if (flag_verbose); disp(sprintf(' %% The parameter nshuffle indicates the (0-based) shuffle-index. ')); end;
if (flag_verbose); disp(sprintf(' %% When nshuffle==0 we do not shuffle, and process the original data-set. ;')); end;
if (flag_verbose); disp(sprintf(' %% This was already done above (as nshuffle==0 by default). ;')); end;
if (flag_verbose); disp(sprintf(' %% Each run should take as long as the original. ;')); end;
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
if (flag_verbose); disp(sprintf(' %% Each patient is colored with a (r,g,b) triplet indicating its membership in (bcA,bcB,bcC). ;')); end;
if (flag_verbose); disp(sprintf(' %% Thus, patients in bcA will be: ')); end;
if (flag_verbose); disp(sprintf(' %%\t red if they are only in bcA, ')); end;
if (flag_verbose); disp(sprintf(' %%\t yellow if they are in bcA+bcB, ')); end;
if (flag_verbose); disp(sprintf(' %%\t magenta if they are in bcA+bcC, ')); end;
if (flag_verbose); disp(sprintf(' %%\t white if they are in bcA+bcB+bcC. ')); end;
if (flag_verbose); disp(sprintf(' %% Similarly, patients in bcB will be mostly green, ;')); end;
if (flag_verbose); disp(sprintf(' %% while patients in bcC will be mostly blue. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note that the bicluster-informed projection singles out bcA, ;')); end;
if (flag_verbose); disp(sprintf(' %% while an uninformed projection (using the principal-components of the full data-set) ;')); end;
if (flag_verbose); disp(sprintf(' %% does not single out bcA, but instead is confounded by the nonspecific bicluster bcB, ;')); end;
if (flag_verbose); disp(sprintf(' %% as well as the mds-biased bicluster bcC (c.f. principal-components 1 and 2, respectively). ;')); end;
if (flag_verbose); disp(sprintf(' %% We use our ''pca_driver'' code to calculate the principal-components (using power iteration). ;')); end;
if flag_disp;
figure(1+nf);nf=nf+1;clf;
%set(gcf,'Position',1+[0,0,1024,1024]);fig80s;
%p_row = 2; p_col = 2; np=0;
figmed;fig80s;
end;%if flag_disp;
p_row = 1; p_col = 2; np=0;
mr_ABC__ = zeros(n_patient,3);
mr_ABC__(rij_case_bcA_,1+0) = 1;
mr_ABC__(rij_case_bcB_,1+1) = 1;
mr_ABC__(rij_ctrl_bcB_,1+1) = 1;
mr_ABC__(rij_case_bcC_,1+2) = 1;
markersize_sml = 4;
markersize_med = 8;
markersize_big = 16;
xdrop_ij_rkeep_ = xdrop_.ij_rkeep_;
xdrop_ij_ckeep_ = xdrop_.ij_ckeep_;
%%%%;
%A_p_ = mean(A_full_n_>0,1);
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
scatter(AZnV_ni0_driver_(:,1),AZnV_ni0_driver_(:,2),16,mr_ABC__,'filled','MarkerEdgeColor','k');
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
scatter(AZnV_nix_driver_(:,1),AZnV_nix_driver_(:,2),16,mr_ABC__,'filled','MarkerEdgeColor','k');
axisnotick; title('bicluster-informed (driver)'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;

if (flag_verbose); disp(sprintf(' %% [finished test_stripped_xxxcluster_fromdisk_ver16]')); end;





