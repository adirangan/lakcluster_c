function lisa = lisa_struct_prs_scatterplot_ver1(specification);
% try: ;
%{
  
  clear;
  dir_bip32 = '/data/rangan/dir_bcc/dir_BIP32_loo_data';
  dir_icuk = sprintf('%s/dir_icuk_scores',dir_bip32);
  dir_icuk_ex = sprintf('%s/dir_icuk_ex_bc1',dir_bip32);
  dir_icuk_on_icuk = sprintf('%s/dir_icuk_on_icuk_loo_proriles',dir_bip32);
  dir_SCZ_MDD = sprintf('%s/dir_SCZ_MDD_PRS_profile_scores',dir_bip32);
  dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PGC_20190328');
  dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c_dev');
  flag_dex_vs_lak = 'dex'; %<-- differentially expressed clustering. ;
  cl_num_arm1 = 4; cl_num_arm2 = 1; %<-- train on platform 4, replicate on platform 1. ;
  cl_num_arm1_mc_A = []; cl_num_arm1_mr_A_full_ = []; cl_num_arm1_mr_Z_full_ = []; cl_num_arm1_string = []; 
  cl_num_arm2_mc_A = []; cl_num_arm2_mr_A_full_ = []; cl_num_arm2_mr_Z_full_ = []; cl_num_arm2_string = []; 
  flag_reverse = 0; %<-- forward bicluster (i.e., case-specific). ;
  gamma = [0.004]; %<-- gamma is the fraction eliminated per iteration. 000 implies a single patient eliminated per iteration. ;
  B_MLT = 34; n_mds = 20; mr_string = '';mc_string = ''; %<-- accurate to 2^(-34), 20 total mds components (but only 2 used). ; No special mc_string. ;
  n_maf = 5; n_cov = 2; %<-- minor-allele-frequency cutoff 25-50, 2 covariates (mds-components) used, repeated twice. ;
  n_scramble = 0; n_shuffle = 0; %<-- no previous bicluster extracted/scrambled first, no random shuffling. ;
  flag_rerun=1; %<-- regenerate mat-file.; 
  pca_b_mlt = 44; pca_tolerance = 1e-2; pca_rank = 2;
  prs_prefix = 'ori';
  n_iteration_stride = 25/2;
  ni_r0_arm2 = 175;
  ni_r1_arm2 = 163;
  specification = struct();
  specification.dir_bip32 = dir_bip32;
  specification.dir_icuk = dir_icuk;
  specification.dir_icuk_ex = dir_icuk_ex;
  specification.dir_icuk_on_icuk = dir_icuk_on_icuk;
  specification.dir_SCZ_MDD = dir_SCZ_MDD;
  specification.dir_trunk = dir_trunk;
  specification.dir_code = dir_code;
  specification.mr_string = mr_string;
specification.mc_string = mc_string;
  specification.cl_num_arm1 = cl_num_arm1;
  specification.cl_num_arm1_mc_A = cl_num_arm1_mc_A;
  specification.cl_num_arm1_mr_A_full_ = cl_num_arm1_mr_A_full_;
  specification.cl_num_arm1_mr_Z_full_ = cl_num_arm1_mr_Z_full_;
  specification.cl_num_arm1_string = cl_num_arm1_string;
  specification.cl_num_arm2 = cl_num_arm2;
  specification.cl_num_arm2_mc_A = cl_num_arm2_mc_A;
  specification.cl_num_arm2_mr_A_full_ = cl_num_arm2_mr_A_full_;
  specification.cl_num_arm2_mr_Z_full_ = cl_num_arm2_mr_Z_full_;
  specification.cl_num_arm2_string = cl_num_arm2_string;
  specification.flag_dex_vs_lak = flag_dex_vs_lak;
  specification.gamma = gamma;
  specification.B_MLT = B_MLT;
  specification.n_mds = n_mds;
  specification.flag_reverse = flag_reverse;
  specification.n_maf = n_maf;
  specification.n_cov = n_cov;
  specification.n_scramble = n_scramble;
  specification.n_shuffle = n_shuffle;
  specification.flag_rerun = flag_rerun;
  specification.pca_b_mlt = pca_b_mlt;
  specification.pca_tolerance = pca_tolerance;
  specification.pca_rank = pca_rank;
  specification.n_iteration_stride = n_iteration_stride;
  specification.ni_r0_arm2 = ni_r0_arm2;
  specification.ni_r1_arm2 = ni_r1_arm2;
  specification.prs_prefix = prs_prefix;
  lisa_struct_prs_scatterplot_ver1(specification) ;

  %}

setup;
dir_bip32 = specification.dir_bip32;
dir_icuk = specification.dir_icuk;
dir_icuk_ex = specification.dir_icuk_ex;
dir_icuk_on_icuk = specification.dir_icuk_on_icuk;
dir_SCZ_MDD = dir_SCZ_MDD;
dir_trunk = specification.dir_trunk;
dir_code = specification.dir_code;
mr_string = specification.mr_string;
mc_string = specification.mc_string;
cl_num_arm1 = specification.cl_num_arm1;
cl_num_arm1_mc_A = specification.cl_num_arm1_mc_A;
cl_num_arm1_mr_A_full_ = specification.cl_num_arm1_mr_A_full_;
cl_num_arm1_mr_Z_full_ = specification.cl_num_arm1_mr_Z_full_;
cl_num_arm1_string = specification.cl_num_arm1_string;
cl_num_arm2 = specification.cl_num_arm2;
cl_num_arm2_mc_A = specification.cl_num_arm2_mc_A;
cl_num_arm2_mr_A_full_ = specification.cl_num_arm2_mr_A_full_;
cl_num_arm2_mr_Z_full_ = specification.cl_num_arm2_mr_Z_full_;
cl_num_arm2_string = specification.cl_num_arm2_string;
flag_dex_vs_lak = specification.flag_dex_vs_lak;
gamma = specification.gamma;
B_MLT = specification.B_MLT;
n_mds = specification.n_mds;
flag_reverse = specification.flag_reverse;
n_maf = specification.n_maf;
n_cov = specification.n_cov;
n_scramble = specification.n_scramble;
n_shuffle = specification.n_shuffle;
flag_rerun = specification.flag_rerun;
pca_b_mlt = specification.pca_b_mlt;
pca_tolerance = specification.pca_tolerance;
pca_rank = specification.pca_rank;
n_iteration_stride = specification.n_iteration_stride;
ni_r0_arm2 = specification.ni_r0_arm2;
ni_r1_arm2 = specification.ni_r1_arm2;
prs_prefix = specification.prs_prefix;

%%%%%%%%;
lisa_arm1_r0 = lisa_struct_make_ver0(mr_string,mc_string,cl_num_arm1,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,0,n_shuffle) ;
lisa_arm1_r0 = lisa_struct_prefix_ver0(lisa_arm1_r0,dir_code,dir_trunk); 
lisa_arm1_r0.nshuffle = 0;  lisa_arm1_r0 = lisa_struct_names_ver0(lisa_arm1_r0); 
lisa_arm1_r0 = lisa_struct_xdrop_ver0(lisa_arm1_r0); lisa_arm1_r0 = lisa_struct_mdsfam_ver0(lisa_arm1_r0); 
%lisa_arm1_r0 = lisa_struct_bim_ver0(lisa_arm1_r0); %<-- this is large and takes a while to load. ;
lisa_arm1_r0 = lisa_struct_mx_ver0(lisa_arm1_r0); lisa_arm1_r0 = lisa_struct_studyindex_ver0(lisa_arm1_r0); 
lisa_arm1_r0 = lisa_struct_trace_ver0(lisa_arm1_r0);
%%%%%%%%;
lisa_arm1_r1 = lisa_struct_make_ver0(mr_string,mc_string,cl_num_arm1,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,1,n_shuffle) ;
lisa_arm1_r1 = lisa_struct_prefix_ver0(lisa_arm1_r1,dir_code,dir_trunk); 
lisa_arm1_r1.nshuffle = 0;  lisa_arm1_r1 = lisa_struct_names_ver0(lisa_arm1_r1); 
lisa_arm1_r1 = lisa_struct_xdrop_ver0(lisa_arm1_r1); lisa_arm1_r1 = lisa_struct_mdsfam_ver0(lisa_arm1_r1); 
%lisa_arm1_r1 = lisa_struct_bim_ver0(lisa_arm1_r1); %<-- this is large and takes a while to load. ;
lisa_arm1_r1 = lisa_struct_mx_ver0(lisa_arm1_r1); lisa_arm1_r1 = lisa_struct_studyindex_ver0(lisa_arm1_r1); 
lisa_arm1_r1 = lisa_struct_trace_ver0(lisa_arm1_r1);
%%%%%%%%;
lisa_arm2 = lisa_struct_make_ver0(mr_string,mc_string,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,0,n_shuffle) ;
lisa_arm2 = lisa_struct_prefix_ver0(lisa_arm2,dir_code,dir_trunk); 
lisa_arm2.nshuffle = 0;  lisa_arm2 = lisa_struct_names_ver0(lisa_arm2); 
lisa_arm2 = lisa_struct_xdrop_ver0(lisa_arm2); lisa_arm2 = lisa_struct_mdsfam_ver0(lisa_arm2); 
%lisa_arm2 = lisa_struct_bim_ver0(lisa_arm2); %<-- this is large and takes a while to load. ;
lisa_arm2 = lisa_struct_mx_ver0(lisa_arm2); lisa_arm2 = lisa_struct_studyindex_ver0(lisa_arm2); 
lisa_arm2 = lisa_struct_trace_ver0(lisa_arm2);
%%%%%%%%;

fname_tmp = sprintf('%s/MDD_SCZ_possible_control_overlap.txt',dir_bip32); fcheck(fname_tmp);
fid=fopen(fname_tmp);
MDD_SCZ_overlap_ = textscan(fid,'%s%d%d%d%d','headerlines',1);
fclose(fid);
MDD_SCZ_overlap_name = 1;
MDD_SCZ_overlap_MDD_exclude = 2;
MDD_SCZ_overlap_MDD_notsure = 3;
MDD_SCZ_overlap_SCZ_exclude = 4;
MDD_SCZ_overlap_SCZ_notsure = 5;
%%%%%%%%;
verbose = 1; n_figure = 1;

if strcmp(prs_prefix,'ori'); str_prs_ = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
end;
if strcmp(prs_prefix,'756'); str_prs_ = {'n3','n2','n1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-9.5','1e-9','1e-8.5','1e-8','1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
end;
if strcmp(prs_prefix,'469'); str_prs_ = {'n3','n2','n1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-9.5','1e-9','1e-8.5','1e-8','1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
end;
if strcmp(prs_prefix,'319'); str_prs_ = {'n3','n2','n1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-9.5','1e-9','1e-8.5','1e-8','1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
end;
if strcmp(prs_prefix,'MDD_exclude'); str_prs_ = {'n3','n2','n1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-9.5','1e-9','1e-8.5','1e-8','1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
end;
if strcmp(prs_prefix,'SCZ_exclude'); str_prs_ = {'n3','n2','n1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-9.5','1e-9','1e-8.5','1e-8','1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
end;
if strcmp(prs_prefix,'MDD_notsure'); str_prs_ = {'n3','n2','n1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-9.5','1e-9','1e-8.5','1e-8','1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
end;
if strcmp(prs_prefix,'SCZ_notsure'); str_prs_ = {'n3','n2','n1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-9.5','1e-9','1e-8.5','1e-8','1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
end;
if strcmp(prs_prefix,'icuk_on_icuk'); str_prs_ = {'n3','n2','n1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-9.5','1e-9','1e-8.5','1e-8','1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
end;

lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,lisa_arm2.cl_num,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver1',lisa_arm1_r0.dir_out_s0_pca)); end;
lisa_arm1_r0.dir_out_s0_pca_jpg = sprintf('%s/dir_jpg',lisa_arm1_r0.dir_out_s0_pca);
if (~exist(lisa_arm1_r0.dir_out_s0_pca_jpg,'dir')); disp(sprintf(' %% creating %s',lisa_arm1_r0.dir_out_s0_pca_jpg)); mkdir(lisa_arm1_r0.dir_out_s0_pca_jpg); end;
lisa_arm1_r1.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r1.dir_out_s0,lisa_arm1_r1.cl_num,cl_num_arm1_string,lisa_arm2.cl_num,cl_num_arm2_string);
if (~exist(lisa_arm1_r1.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver1',lisa_arm1_r1.dir_out_s0_pca)); end;
lisa_arm1_r1.dir_out_s0_pca_jpg = sprintf('%s/dir_jpg',lisa_arm1_r1.dir_out_s0_pca);
if (~exist(lisa_arm1_r1.dir_out_s0_pca_jpg,'dir')); disp(sprintf(' %% creating %s',lisa_arm1_r1.dir_out_s0_pca_jpg)); mkdir(lisa_arm1_r1.dir_out_s0_pca_jpg); end;

%%%%%%%%;
dir_out_s0_pca_r0_arm1 = lisa_arm1_r0.dir_out_s0_pca;
dir_out_s0_pca_r1_arm1 = lisa_arm1_r1.dir_out_s0_pca;
%%%%%%%%%%%%%%%%;
% load prs scores. ;
%%%%%%%%%%%%%%%%;
n_prs = length(str_prs_);
n_patient_prs_ = zeros(n_prs,1);
prs_fam_ = cell(n_prs,1);
prs_iid_ = cell(n_prs,1);
prs_phe_ = cell(n_prs,1);
prs_sco_ = cell(n_prs,1);
prs_inc_ = cell(n_prs,1);
iprs_ = cell(n_prs,1);
%%%%%%%%;
for nprs=1:n_prs;
%%%%%%%%;
n_patient_prs_(nprs)=0;
for nstudy_arm2=1:lisa_arm2.n_study;
tmp_study_name = lisa_arm2.study_name_{nstudy_arm2};
if (tmp_study_name(end)~='c'); tmp_study_name = tmp_study_name(1:end-1); end;
if (strcmp(tmp_study_name(1:4),'bip_')); tmp_study_name = tmp_study_name(5:end); end;
tmp_study_name = tmp_study_name(1:4);
if strcmp(prs_prefix,'ori'); tmp_fname = sprintf('%s/icuk_dos_bip_%s_eur_sr-qc.hg19.ch.fl.out.dosage_imp_scores.txt.S%s.profile',dir_icuk,tmp_study_name,str_prs_{nprs}); end;
if strcmp(prs_prefix,'756'); tmp_fname = sprintf('%s/icuk_geno_ex_bicl1_%s_on_dos_%s.S%s.profile',dir_icuk_ex,prs_prefix,tmp_study_name,str_prs_{nprs}); end;
if strcmp(prs_prefix,'469'); tmp_fname = sprintf('%s/icuk_geno_ex_bicl1_%s_on_dos_%s.S%s.profile',dir_icuk_ex,prs_prefix,tmp_study_name,str_prs_{nprs}); end;
if strcmp(prs_prefix,'319'); tmp_fname = sprintf('%s/icuk_geno_ex_bicl1_%s_on_dos_%s.S%s.profile',dir_icuk_ex,prs_prefix,tmp_study_name,str_prs_{nprs}); end;
if strcmp(prs_prefix,'MDD_exclude'); tmp_fname = sprintf('%s/%s_imp_%s_eur_sr-qc.hg19.ch.fl.out.dosage_imp_scores.txt.S%s.profile',dir_SCZ_MDD,prs_prefix(1:3),tmp_study_name,str_prs_{nprs}); end;
if strcmp(prs_prefix,'SCZ_exclude'); tmp_fname = sprintf('%s/%s_imp_%s_eur_sr-qc.hg19.ch.fl.out.dosage_imp_scores.txt.S%s.profile',dir_SCZ_MDD,prs_prefix(1:3),tmp_study_name,str_prs_{nprs}); end;
if strcmp(prs_prefix,'MDD_notsure'); tmp_fname = sprintf('%s/%s_imp_%s_eur_sr-qc.hg19.ch.fl.out.dosage_imp_scores.txt.S%s.profile',dir_SCZ_MDD,prs_prefix(1:3),tmp_study_name,str_prs_{nprs}); end;
if strcmp(prs_prefix,'SCZ_notsure'); tmp_fname = sprintf('%s/%s_imp_%s_eur_sr-qc.hg19.ch.fl.out.dosage_imp_scores.txt.S%s.profile',dir_SCZ_MDD,prs_prefix(1:3),tmp_study_name,str_prs_{nprs}); end;
if strcmp(prs_prefix,'icuk_on_icuk'); tmp_fname = sprintf('%s/daner_PGC_BIP32b_mds7a.loo.no.%s.gz.clumped.xmhcicuk_eur_sr-qc.hg19.ch.fl.out.dosage_leave_out_icuk_scores.txt.S%s.profile',dir_icuk_on_icuk,tmp_study_name,str_prs_{nprs}); end;
fcheck(tmp_fname); 
fid = fopen(tmp_fname);
tmp_prs_ = textscan(fid,'%s%s%d%f','headerlines',1);
fclose(fid);
tmp_n = length(tmp_prs_{1});
flag_exclude=0;
if strcmp(prs_prefix,'MDD_exclude'); tmp_ij = find(strcmp(MDD_SCZ_overlap_{MDD_SCZ_overlap_name},tmp_study_name)); flag_exclude = MDD_SCZ_overlap_{MDD_SCZ_overlap_MDD_exclude}(tmp_ij); end;
if strcmp(prs_prefix,'SCZ_exclude'); tmp_ij = find(strcmp(MDD_SCZ_overlap_{MDD_SCZ_overlap_name},tmp_study_name)); flag_exclude = MDD_SCZ_overlap_{MDD_SCZ_overlap_SCZ_exclude}(tmp_ij); end;
if strcmp(prs_prefix,'MDD_notsure'); tmp_ij = find(strcmp(MDD_SCZ_overlap_{MDD_SCZ_overlap_name},tmp_study_name)); flag_exclude = MDD_SCZ_overlap_{MDD_SCZ_overlap_MDD_exclude}(tmp_ij) | MDD_SCZ_overlap_{MDD_SCZ_overlap_MDD_notsure}(tmp_ij); end;
if strcmp(prs_prefix,'SCZ_notsure'); tmp_ij = find(strcmp(MDD_SCZ_overlap_{MDD_SCZ_overlap_name},tmp_study_name)); flag_exclude = MDD_SCZ_overlap_{MDD_SCZ_overlap_SCZ_exclude}(tmp_ij) | MDD_SCZ_overlap_{MDD_SCZ_overlap_SCZ_notsure}(tmp_ij); end;
if flag_exclude; disp(sprintf(' %% excluding %s',tmp_study_name)); end;
ni = 1;
prs_fam_{nprs}(n_patient_prs_(nprs) + (1:tmp_n)) = tmp_prs_{ni}(1:tmp_n); ni=ni+1;
prs_iid_{nprs}(n_patient_prs_(nprs) + (1:tmp_n)) = tmp_prs_{ni}(1:tmp_n); ni=ni+1;
prs_phe_{nprs}(n_patient_prs_(nprs) + (1:tmp_n)) = tmp_prs_{ni}(1:tmp_n); ni=ni+1;
prs_sco_{nprs}(n_patient_prs_(nprs) + (1:tmp_n)) = tmp_prs_{ni}(1:tmp_n); ni=ni+1;
prs_inc_{nprs}(n_patient_prs_(nprs) + (1:tmp_n)) = ones(1,tmp_n)*(~flag_exclude);
n_patient_prs_(nprs) = n_patient_prs_(nprs) + tmp_n;
end;%for nstudy_arm2=1:lisa_arm2.n_study;
prs_name_{nprs} = cell(n_patient_prs_(nprs),1); 
for np=1:n_patient_prs_(nprs); prs_name_{nprs}{np} = sprintf('%s%s%s',prs_fam_{nprs}{np},'&',prs_iid_{nprs}{np}); end;%for np=1:n_patient_prs_(nprs);
if (flag_reverse==1); prs_sco_{nprs} = -prs_sco_{nprs}; end;
%%%%%%%%%%%%%%%%;
disp(sprintf(' %% nprs %d: %s: intersection %d',nprs,tmp_fname,length(intersect(prs_name_{nprs},lisa_arm2.fam_name_))));
%%%%%%%%%%%%%%%%;
[~,ifam_,iprs_{nprs}] = intersect(lisa_arm2.fam_name_,prs_name_{nprs},'stable'); if (length(ifam_)<length(lisa_arm2.fam_name_)); disp('Warning! fam_name_ not a subset of prs_name_{nprs}'); end; 
tmp_f_ = cast(lisa_arm2.fam_dvx_,'double'); tmp_p_ = transpose(cast(prs_phe_{nprs}(iprs_{nprs}),'double'));
tmp_ij_ = find( ( (tmp_f_==1) | (tmp_f_==2) ) & ( (tmp_p_==1) | (tmp_p_==2) ) );
disp(sprintf('nprs %d: overlap %d/%d: phe error: %0.16f',nprs,length(tmp_ij_),length(tmp_f_),norm(tmp_f_(tmp_ij_)-tmp_p_(tmp_ij_),'fro')));
%%%%%%%%;
end;%for nprs=1:n_prs;

%%%%%%%%%%%%%%%%;
% calculate auc and p-value for arm2;
%%%%%%%%%%%%%%%%;
prs_min_ = zeros(n_prs,1); prs_max_ = zeros(n_prs,1);
for nprs=1:n_prs;
disp(sprintf(' %% nprs %d: str %s',nprs,str_prs_{nprs}));
find_x_ = find(lisa_arm2.mr_X_ & transpose(prs_inc_{nprs})); find_d_ = find(lisa_arm2.mr_D_ & transpose(prs_inc_{nprs}));
prs_min = min(prs_sco_{nprs}); prs_max = max(prs_sco_{nprs});
prs_min_(nprs) = prs_min; prs_max_(nprs) = prs_max;
end;%for nprs=1:n_prs;
tmp_fname = sprintf('%s/pca_ni%d_tst%d_k2_B44_V_arm2_.mda',dir_out_s0_pca_r0_arm1,ni_r0_arm2,cl_num_arm2); fcheck(tmp_fname); V_r0_arm2_ = mda_read_r8(tmp_fname);
tmp_fname = sprintf('%s/pca_proj_arm2_ni%d_tst%d_k2_B44_AnV_.mda',dir_out_s0_pca_r0_arm1,ni_r0_arm2,cl_num_arm2); fcheck(tmp_fname); DnV_r0_arm2_ = mda_read_r8(tmp_fname);
tmp_fname = sprintf('%s/pca_proj_arm2_ni%d_tst%d_k2_B44_ZnV_.mda',dir_out_s0_pca_r0_arm1,ni_r0_arm2,cl_num_arm2); fcheck(tmp_fname); XnV_r0_arm2_ = mda_read_r8(tmp_fname);
tmp_mds_x_ = lisa_arm2.mds_sort_(find_x_,[1:6,19]); tmp_mds_d_ = lisa_arm2.mds_sort_(find_d_,[1:6,19]);
tmp_r0_auc19 = cauc_0(XnV_r0_arm2_(find_x_,1),DnV_r0_arm2_(find_d_,1),tmp_mds_x_,tmp_mds_d_);
if (tmp_r0_auc19<0.5); disp(sprintf(' %% flipping nir0_arm2_ %d',nir0_arm2_)); DnV_r0_arm2_ = -DnV_r0_arm2_; XnV_r0_arm2_ = -XnV_r0_arm2_; end;
tmp_fname = sprintf('%s/pca_ni%d_tst%d_k2_B44_V_arm2_.mda',dir_out_s0_pca_r1_arm1,ni_r1_arm2,cl_num_arm2); fcheck(tmp_fname); V_r1_arm2_ = mda_read_r8(tmp_fname);
tmp_fname = sprintf('%s/pca_proj_arm2_ni%d_tst%d_k2_B44_AnV_.mda',dir_out_s0_pca_r1_arm1,ni_r1_arm2,cl_num_arm2); fcheck(tmp_fname); DnV_r1_arm2_ = mda_read_r8(tmp_fname);
tmp_fname = sprintf('%s/pca_proj_arm2_ni%d_tst%d_k2_B44_ZnV_.mda',dir_out_s0_pca_r1_arm1,ni_r1_arm2,cl_num_arm2); fcheck(tmp_fname); XnV_r1_arm2_ = mda_read_r8(tmp_fname);
tmp_mds_x_ = lisa_arm2.mds_sort_(find_x_,[1:6,19]); tmp_mds_d_ = lisa_arm2.mds_sort_(find_d_,[1:6,19]);
tmp_r1_auc19 = cauc_0(XnV_r1_arm2_(find_x_,1),DnV_r1_arm2_(find_d_,1),tmp_mds_x_,tmp_mds_d_);
if(tmp_r1_auc19<0.5); disp(sprintf(' %% flipping nir1_arm2_ %d',nir1_arm2_)); DnV_r1_arm2_ = -DnV_r1_arm2_; XnV_r1_arm2_ = -XnV_r1_arm2_; end;

%%%%%%%%%%%%%%%%;
% calculate cauc and p-value for arm2 under top+bottom thresholding;
%%%%%%%%%%%%%%%%;
fname_mat = sprintf('%s/prs_scatterplot_ver1_%s_trn%d%s_tst%d%s.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if ( flag_rerun==0 &  exist(fname_mat,'file')); 
disp(sprintf(' %% found %s, loading',fname_mat));
load(fname_mat,...
     'n_rx_t','prctile_','prctile_rem_',...
     'auc19_prs_t_',...
     'p_x_rem_','p_d_rem_',...
     'rx_min_','rx_max_','rx_t__',...
     'n_omega','omega_');
end;%if ( exist(fname_mat,'file')); 
if ( flag_rerun==1 | ~exist(fname_mat,'file')); 
disp(sprintf(' %% could not find %s, creating',fname_mat));
n_omega = 24; omega_ = linspace(0,pi/2,n_omega);
n_rx_t = 20;
nprs=n_prs; find_x_ = find(lisa_arm2.mr_X_ & transpose(prs_inc_{nprs})); find_d_ = find(lisa_arm2.mr_D_ & transpose(prs_inc_{nprs}));
auc19_prs_t_ = zeros(n_omega,1+n_rx_t,n_prs);
p_x_rem_ = zeros(n_omega,1+n_rx_t);
p_d_rem_ = zeros(n_omega,1+n_rx_t);
rx_min_ = zeros(n_omega,1);
rx_max_ = zeros(n_omega,1);
rx_t__ = zeros(n_omega,1+2*n_rx_t);
%%%%%%%%;
tmp_Xx_ = XnV_r0_arm2_(find_x_,1); tmp_Dx_ = DnV_r0_arm2_(find_d_,1); tmp_Ax_ = [tmp_Xx_;tmp_Dx_]; tmp_avg = mean(tmp_Ax_); tmp_std = std(tmp_Ax_);
tmp_Xx_ = (tmp_Xx_ - tmp_avg)/tmp_std; tmp_Dx_ = (tmp_Dx_ - tmp_avg)/tmp_std;
tmp_Xy_ = XnV_r1_arm2_(find_x_,1); tmp_Dy_ = DnV_r1_arm2_(find_d_,1); tmp_Ay_ = [tmp_Xy_;tmp_Dy_]; tmp_avg = mean(tmp_Ay_); tmp_std = std(tmp_Ay_);
tmp_Xy_ = (tmp_Xy_ - tmp_avg)/tmp_std; tmp_Dy_ = (tmp_Dy_ - tmp_avg)/tmp_std;
for nomega=1:n_omega;
omega = omega_(nomega);
disp(sprintf(' %% nomega %d/%d omega %d',nomega,n_omega,round(omega*180/pi)));
tmp_Xw_ = cos(omega)*tmp_Xx_ + sin(omega)*tmp_Xy_;
tmp_Dw_ = cos(omega)*tmp_Dx_ + sin(omega)*tmp_Dy_;
rx_min_(nomega) = prctile([tmp_Xw_;tmp_Dw_], 5);
rx_max_(nomega) = prctile([tmp_Xw_;tmp_Dw_], 95);
prctile_ = linspace(0,100,1+2*n_rx_t);
prctile_rem_ = linspace(100,0,1+n_rx_t);
rx_t__(nomega,:) = prctile([tmp_Xw_;tmp_Dw_] , prctile_ );
for nrx_t=0:n_rx_t;
rx_t_min = rx_t__(nomega,1 + n_rx_t - nrx_t); 
rx_t_max = rx_t__(nomega,1 + n_rx_t + nrx_t); 
find_x_t_ = find( ( tmp_Xw_<=rx_t_min) | ( tmp_Xw_>=rx_t_max ) );
find_d_t_ = find( ( tmp_Dw_<=rx_t_min) | ( tmp_Dw_>=rx_t_max ) );
p_x_rem_(nomega,1+nrx_t) = length(find_x_t_)/length(find_x_);
p_d_rem_(nomega,1+nrx_t) = length(find_d_t_)/length(find_d_);
for nprs=1:n_prs;
tmp_prs_x_t_ = prs_sco_{nprs}(iprs_{nprs}(find_x_(find_x_t_))); tmp_prs_d_t_ = prs_sco_{nprs}(iprs_{nprs}(find_d_(find_d_t_)));
tmp_mds_x_t_ = lisa_arm2.mds_sort_(find_x_(find_x_t_),[1:6,19]); tmp_mds_d_t_ = lisa_arm2.mds_sort_(find_d_(find_d_t_),[1:6,19]);
if (length(find_x_t_)==0 & length(find_d_t_)==0); auc19_prs_t_(nomega,1+nrx_t,nprs) = 0.5; end;
if (length(find_x_t_)==0 & length(find_d_t_)> 0); auc19_prs_t_(nomega,1+nrx_t,nprs) = 1.0; end;
if (length(find_x_t_)> 0 & length(find_d_t_)==0); auc19_prs_t_(nomega,1+nrx_t,nprs) = 0.0; end;
if (length(find_x_t_)> 0 & length(find_d_t_)> 0); tmp_cauc = cauc_0(tmp_prs_x_t_,tmp_prs_d_t_,tmp_mds_x_t_,tmp_mds_d_t_); auc19_prs_t_(nomega,1+nrx_t,nprs) = tmp_cauc; end;
end;%for nprs=1:n_prs;
end;%for nrx_t=0:n_rx_t;
end;%for nomega=1:n_omega;
%%%%%%%%;
disp(sprintf(' %% saving %s',fname_mat)); 
save(fname_mat,...
     'n_rx_t','prctile_','prctile_rem_',...
     'auc19_prs_t_',...
     'p_x_rem_','p_d_rem_',...
     'rx_min_','rx_max_','rx_t__',...
     'n_omega','omega_');
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%%%%%%%%%;
% plot cauc for arm2 under top+bottom thresholding. ;
%%%%%%%%%%%%%%%%;
flag_replot = 1;
if strcmp(prs_prefix,'ori'); nprs_sub_ = [7:16]; end;
if strcmp(prs_prefix,'756'); nprs_sub_ = [11:20]; end;
if strcmp(prs_prefix,'469'); nprs_sub_ = [11:20]; end;
if strcmp(prs_prefix,'319'); nprs_sub_ = [11:20]; end;
if strcmp(prs_prefix,'MDD_exclude'); nprs_sub_ = [11:20]; end;
if strcmp(prs_prefix,'SCZ_exclude'); nprs_sub_ = [11:20]; end;
if strcmp(prs_prefix,'MDD_notsure'); nprs_sub_ = [11:20]; end;
if strcmp(prs_prefix,'SCZ_notsure'); nprs_sub_ = [11:20]; end;
if strcmp(prs_prefix,'icuk_on_icuk'); nprs_sub_ = [11:20]; end;
n_tick = 11; 
if strcmp(prs_prefix,'ori'); cauclim_ = [0.500,0.625]; end;
if strcmp(prs_prefix,'756'); cauclim_ = [0.500,0.625]; end;
if strcmp(prs_prefix,'469'); cauclim_ = [0.500,0.625]; end;
if strcmp(prs_prefix,'319'); cauclim_ = [0.500,0.625]; end;
if strcmp(prs_prefix,'MDD_exclude'); cauclim_ = [0.500,0.700]; end;
if strcmp(prs_prefix,'SCZ_exclude'); cauclim_ = [0.500,0.750]; end;
if strcmp(prs_prefix,'MDD_notsure'); cauclim_ = [0.500,0.700]; end;
if strcmp(prs_prefix,'SCZ_notsure'); cauclim_ = [0.500,0.750]; end;
if strcmp(prs_prefix,'icuk_on_icuk'); cauclim_ = [0.500,0.750]; end;
if flag_reverse==0; str_reverse = 'D'; end;
if flag_reverse==1; str_reverse = 'X'; end;
lisa_arm1_r0.dir_out_s0_prs_jpg = sprintf('%s/dir_prs_jpg',lisa_arm1_r0.dir_out_s0_pca);
fname_fig = sprintf('%s/prs_tbthreshold_ver1_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.jpg',lisa_arm1_r0.dir_out_s0_prs_jpg,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
if (~flag_replot & exist(fname_fig,'file')); disp(sprintf(' %% %s found, not creating',fname_fig)); end;
if ( flag_replot | ~exist(fname_fig,'file')); disp(sprintf(' %% %s not found, creating',fname_fig)); 
flag_disp=1;
if flag_disp;
figure(n_figure); clf; n_figure = n_figure+1;
colormap(colormap_beach());
for nprs_sub=1:length(nprs_sub_);
nprs = nprs_sub_(nprs_sub);
%%%%%%%%;
subplot(4,length(nprs_sub_)/2,nprs_sub + 1*length(nprs_sub_));
imagesc(squeeze(auc19_prs_t_(:,1:end-1,nprs)),cauclim_);
tmp_tick_rx_ = round(linspace(1,1+n_rx_t,n_tick)); 
xlabel('rx threshold (p_rem)','interpreter','none'); 
set(gca,'XTick',tmp_tick_rx_,'XTickLabel',prctile_rem_(1:end-1));
ylabel('omega'); 
set(gca,'YTick',1:2:n_omega,'YTickLabel',round(omega_(1:2:n_omega)*180/pi));
xtickangle(90);
title(sprintf('%s snp-p<%s cauc [%0.2f %0.2f]',str_reverse,str_prs_cutoff_{nprs},cauclim_));
%%%%%%%%;
subplot(4,length(nprs_sub_)/2,nprs_sub + 0*length(nprs_sub_));
hold on;
plot(0:n_rx_t-1,squeeze(auc19_prs_t_(1  ,1:end-1,nprs)),'r.-','MarkerSize',25,'LineWidth',2); 
plot(0:n_rx_t-1,squeeze(auc19_prs_t_(end,1:end-1,nprs)),'b.-','MarkerSize',25,'LineWidth',2); 
hold off;
xlim([0,n_rx_t]); ylim(cauclim_); legend('bc1','bc2','Location','NorthWest');
tmp_tick_ = round(linspace(1,1+n_rx_t,n_tick)); %tmp_tick_ = tmp_tick_(end:-1:1);
tmp_tick_val_ = tmp_tick_; 
tmp_tick_label_ = num2str(transpose(prctile_rem_(tmp_tick_)),'%d');
set(gca,'XTick',tmp_tick_val_,'XTickLabel',tmp_tick_label_);
xlabel('rx threshold (p_rem)','interpreter','none'); 
xtickangle(90);
ylabel('cauc');
%title(sprintf('%s S%s cauc',str_reverse,str_prs_{nprs}));
title(sprintf('%s snp-p<%s cauc',str_reverse,str_prs_cutoff_{nprs}));
%%%%%%%%;
end;%for nprs=1:n_prs;
set(gcf,'Position',1+[0,0,1024*2,1024]);
print('-djpeg',fname_fig);
disp(sprintf(' %% writing to %s',fname_fig));
end;%if flag_disp;
end;%if (~exist(fname_fig,'file'));

%%%%%%%%%%%%%%%%;
% plot scatterplots. ;
%%%%%%%%%%%%%%%%;
c_X_ = [0.0,0.5,1.0];
c_D_ = [1.0,0.5,0.0];
size_marker = 10;
n_h = 32; 
DnV_r0_arm2_min = min(min(DnV_r0_arm2_(find_d_,1)),min(DnV_r0_arm2_(find_d_,2))); DnV_r0_arm2_max = max(max(DnV_r0_arm2_(find_d_,1)),max(DnV_r0_arm2_(find_d_,2)));
XnV_r0_arm2_min = min(min(XnV_r0_arm2_(find_x_,1)),min(XnV_r0_arm2_(find_x_,2))); XnV_r0_arm2_max = max(max(XnV_r0_arm2_(find_x_,1)),max(XnV_r0_arm2_(find_x_,2)));
V_r0_arm2_min = min(DnV_r0_arm2_min,XnV_r0_arm2_min); V_r0_arm2_max = max(DnV_r0_arm2_max,XnV_r0_arm2_max);
DnV_r1_arm2_min = min(min(DnV_r1_arm2_(find_d_,1)),min(DnV_r1_arm2_(find_d_,2))); DnV_r1_arm2_max = max(max(DnV_r1_arm2_(find_d_,1)),max(DnV_r1_arm2_(find_d_,2)));
XnV_r1_arm2_min = min(min(XnV_r1_arm2_(find_x_,1)),min(XnV_r1_arm2_(find_x_,2))); XnV_r1_arm2_max = max(max(XnV_r1_arm2_(find_x_,1)),max(XnV_r1_arm2_(find_x_,2)));
V_r1_arm2_min = min(DnV_r1_arm2_min,XnV_r1_arm2_min); V_r1_arm2_max = max(DnV_r1_arm2_max,XnV_r1_arm2_max);
%
%%%%%%%%;
figure();clf; %<-- compare bc1 pca1 with  bc1 pca2 and prs ;
%%%%%%%%;
subplot(2,2,1);hold on;
plot3(XnV_r0_arm2_(find_x_,1),XnV_r0_arm2_(find_x_,2),prs_sco_{nprs}(iprs_{nprs}(find_x_)),'.','Color',c_X_,'MarkerSize',size_marker); 
plot3(DnV_r0_arm2_(find_d_,1),DnV_r0_arm2_(find_d_,2),prs_sco_{nprs}(iprs_{nprs}(find_d_)),'.','Color',c_D_,'MarkerSize',size_marker); 
hold off;
xlabel('r0 pca1'); ylabel('r0 pca2'); zlabel('prs');
xlim([V_r0_arm2_min,V_r0_arm2_max]);
ylim([V_r0_arm2_min,V_r0_arm2_max]);
zlim([prs_min,prs_max]);
axis vis3d;
%%%%%%%%;
subplot(2,2,2); hold on;
tmp_xl_ = [V_r0_arm2_min,V_r0_arm2_max];
tmp_yl_ = [V_r0_arm2_min,V_r0_arm2_max];
tmp_x_ = XnV_r0_arm2_(find_x_,1);
tmp_y_ = XnV_r0_arm2_(find_x_,2);
tmp_X_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_x_ = DnV_r0_arm2_(find_d_,1);
tmp_y_ = DnV_r0_arm2_(find_d_,2);
tmp_D_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
colormap(colormap_pm());
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel('r0 pca1');ylabel('r0 pca2');
set(gca,'XTick',[],'YTick',[]);
%%%%%%%%;
subplot(2,2,3); hold on;
tmp_xl_ = [V_r0_arm2_min,V_r0_arm2_max];
tmp_yl_ = [prs_min,prs_max];
tmp_x_ = XnV_r0_arm2_(find_x_,1);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_x_));
tmp_X_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_x_ = DnV_r0_arm2_(find_d_,1);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_d_));
tmp_D_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
colormap(colormap_pm());
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel('r0 pca1');ylabel('prs');
set(gca,'XTick',[],'YTick',[]);
%%%%%%%%;
subplot(2,2,4); hold on;
tmp_xl_ = [V_r0_arm2_min,V_r0_arm2_max];
tmp_yl_ = [prs_min,prs_max];
tmp_x_ = XnV_r0_arm2_(find_x_,2);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_x_));
tmp_X_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_x_ = DnV_r0_arm2_(find_d_,2);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_d_));
tmp_D_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
colormap(colormap_pm());
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel('r0 pca2');ylabel('prs');
set(gca,'XTick',[],'YTick',[]);
%%%%%%%%;
%
%%%%%%%%;
figure();clf; %<-- compare bc2 pca1 with  bc2 pca2 and prs ;
%%%%%%%%;
subplot(2,2,1);hold on;
plot3(XnV_r1_arm2_(find_x_,1),XnV_r1_arm2_(find_x_,2),prs_sco_{nprs}(iprs_{nprs}(find_x_)),'.','Color',c_X_,'MarkerSize',size_marker); 
plot3(DnV_r1_arm2_(find_d_,1),DnV_r1_arm2_(find_d_,2),prs_sco_{nprs}(iprs_{nprs}(find_d_)),'.','Color',c_D_,'MarkerSize',size_marker); 
hold off;
xlabel('r1 pca1'); ylabel('r1 pca2'); zlabel('prs');
xlim([V_r1_arm2_min,V_r1_arm2_max]);
ylim([V_r1_arm2_min,V_r1_arm2_max]);
zlim([prs_min,prs_max]);
axis vis3d;
%%%%%%%%;
subplot(2,2,2); hold on;
tmp_xl_ = [V_r1_arm2_min,V_r1_arm2_max];
tmp_yl_ = [V_r1_arm2_min,V_r1_arm2_max];
tmp_x_ = XnV_r1_arm2_(find_x_,1);
tmp_y_ = XnV_r1_arm2_(find_x_,2);
tmp_X_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_x_ = DnV_r1_arm2_(find_d_,1);
tmp_y_ = DnV_r1_arm2_(find_d_,2);
tmp_D_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
colormap(colormap_pm());
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel('r1 pca1');ylabel('r1 pca2');
set(gca,'XTick',[],'YTick',[]);
%%%%%%%%;
subplot(2,2,3); hold on;
tmp_xl_ = [V_r1_arm2_min,V_r1_arm2_max];
tmp_yl_ = [prs_min,prs_max];
tmp_x_ = XnV_r1_arm2_(find_x_,1);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_x_));
tmp_X_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_x_ = DnV_r1_arm2_(find_d_,1);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_d_));
tmp_D_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
colormap(colormap_pm());
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel('r1 pca1');ylabel('prs');
set(gca,'XTick',[],'YTick',[]);
%%%%%%%%;
subplot(2,2,4); hold on;
tmp_xl_ = [V_r1_arm2_min,V_r1_arm2_max];
tmp_yl_ = [prs_min,prs_max];
tmp_x_ = XnV_r1_arm2_(find_x_,2);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_x_));
tmp_X_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_x_ = DnV_r1_arm2_(find_d_,2);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_d_));
tmp_D_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
colormap(colormap_pm());
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel('r1 pca2');ylabel('prs');
set(gca,'XTick',[],'YTick',[]);
%%%%%%%%;
%
%%%%%%%%;
figure();clf; %<-- compare bc1 pca1 with  bc2 pca1 and prs ;
%%%%%%%%;
subplot(2,2,1);hold on;
plot3(XnV_r0_arm2_(find_x_,1),XnV_r1_arm2_(find_x_,1),prs_sco_{nprs}(iprs_{nprs}(find_x_)),'.','Color',c_X_,'MarkerSize',size_marker); 
plot3(DnV_r0_arm2_(find_d_,1),DnV_r1_arm2_(find_d_,1),prs_sco_{nprs}(iprs_{nprs}(find_d_)),'.','Color',c_D_,'MarkerSize',size_marker); 
hold off;
xlabel('r0 pca1'); ylabel('r1 pca1'); zlabel('prs');
xlim([V_r0_arm2_min,V_r0_arm2_max]);
ylim([V_r1_arm2_min,V_r1_arm2_max]);
zlim([prs_min,prs_max]);
axis vis3d;
%%%%%%%%;
subplot(2,2,2); hold on;
tmp_xl_ = [V_r0_arm2_min,V_r0_arm2_max];
tmp_yl_ = [V_r1_arm2_min,V_r1_arm2_max];
tmp_x_ = XnV_r0_arm2_(find_x_,1);
tmp_y_ = XnV_r1_arm2_(find_x_,1);
tmp_X_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_x_ = DnV_r0_arm2_(find_d_,1);
tmp_y_ = DnV_r1_arm2_(find_d_,1);
tmp_D_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
colormap(colormap_pm());
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel('r0 pca1');ylabel('r1 pca1');
set(gca,'XTick',[],'YTick',[]);
%%%%%%%%;
subplot(2,2,3); hold on;
tmp_xl_ = [V_r0_arm2_min,V_r0_arm2_max];
tmp_yl_ = [prs_min,prs_max];
tmp_x_ = XnV_r0_arm2_(find_x_,1);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_x_));
tmp_X_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_x_ = DnV_r0_arm2_(find_d_,1);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_d_));
tmp_D_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
colormap(colormap_pm());
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel('r0 pca1');ylabel('prs');
set(gca,'XTick',[],'YTick',[]);
%%%%%%%%;
subplot(2,2,4); hold on;
tmp_xl_ = [V_r1_arm2_min,V_r1_arm2_max];
tmp_yl_ = [prs_min,prs_max];
tmp_x_ = XnV_r1_arm2_(find_x_,1);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_x_));
tmp_X_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_x_ = DnV_r1_arm2_(find_d_,1);
tmp_y_ = prs_sco_{nprs}(iprs_{nprs}(find_d_));
tmp_D_ = hist2d_0(tmp_x_,tmp_y_,n_h,n_h,tmp_xl_,tmp_yl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
colormap(colormap_pm());
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel('r1 pca1');ylabel('prs');
set(gca,'XTick',[],'YTick',[]);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Note that bc1 pca1 and bc2 pca1 seem to strongly separate the cases from controls. ;
% Here we try and plot a linear combination of bc1 pca1 and bc2 pca1 versus the prs-score. ;
% omega=0 corresponds to bc1 pca1, and omega=pi/2 corresponds to bc2 pca1. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_omega = 24; omega_ = linspace(0,pi/2,n_omega);
r2_00m_ = zeros(n_omega,n_prs);
r2_0zm_ = zeros(n_omega,n_prs);
r2_w0m_ = zeros(n_omega,n_prs);
r2_wzm_ = zeros(n_omega,n_prs);
for nprs=1:n_prs;
prs_min = prs_min_(nprs); prs_max = prs_max_(nprs);
figure;clf;
for nomega=1:n_omega;
omega = omega_(nomega);
tmp_xl_ = [V_r0_arm2_min,V_r0_arm2_max];
tmp_yl_ = [V_r1_arm2_min,V_r1_arm2_max];
tmp_wl_ = [0,1];
tmp_zl_ = [prs_min,prs_max];
tmp_Xx_ = (XnV_r0_arm2_(find_x_,1) - tmp_xl_(1)) / diff(tmp_xl_);
tmp_Xy_ = (XnV_r1_arm2_(find_x_,1) - tmp_yl_(1)) / diff(tmp_yl_);
tmp_Xw_ = cos(omega)*tmp_Xx_ + sin(omega)*tmp_Xy_;
tmp_Xz_ = transpose(prs_sco_{nprs}(iprs_{nprs}(find_x_)));
tmp_X_ = hist2d_0(tmp_Xw_,tmp_Xz_,n_h,n_h,tmp_wl_,tmp_zl_);
tmp_Dx_ = (DnV_r0_arm2_(find_d_,1) - tmp_xl_(1)) / diff(tmp_xl_);
tmp_Dy_ = (DnV_r1_arm2_(find_d_,1) - tmp_yl_(1)) / diff(tmp_yl_);
tmp_Dw_ = cos(omega)*tmp_Dx_ + sin(omega)*tmp_Dy_;
tmp_Dz_ = transpose(prs_sco_{nprs}(iprs_{nprs}(find_d_)));
tmp_D_ = hist2d_0(tmp_Dw_,tmp_Dz_,n_h,n_h,tmp_wl_,tmp_zl_);
tmp_E_ = tmp_D_/sum(tmp_D_(:)) - tmp_X_/sum(tmp_X_(:));
tmp_outcome_ = [ones(size(tmp_Dw_,1),1) ; 1 + ones(size(tmp_Xw_,1),1)];
%[beta19_arm2,dev19_arm2,stats19_arm2] = mnrfit([tmp_Dw_ , tmp_mds_d_ ; tmp_Xw_ , tmp_mds_x_ ] , tmp_outcome_);
%[beta19_arm2,dev19_arm2,stats19_arm2] = mnrfit([tmp_Dw_ , tmp_Dz_ , tmp_mds_d_ ; tmp_Xw_ , tmp_Xz_ , tmp_mds_x_ ] , tmp_outcome_);
%[beta19_arm2,dev19_arm2,stats19_arm2] = mnrfit([tmp_Dw_ , tmp_Dz_ ; tmp_Xw_ , tmp_Xz_ ] , tmp_outcome_);
%disp(sprintf(' %% nomega %d/%d omega %0.2f : beta [%+0.3f %+0.3f] p [%0.8f %0.8f]',nomega,n_omega,omega*180/pi,beta19_arm2([2,3]),stats19_arm2.p([2,3])));
r2_00m_(nomega,nprs) = nagelkerke_r2_ver0([tmp_mds_d_ ; tmp_mds_x_ ] , tmp_outcome_);
r2_0zm_(nomega,nprs) = nagelkerke_r2_ver0([tmp_Dz_ , tmp_mds_d_ ; tmp_Xz_ , tmp_mds_x_ ] , tmp_outcome_);
r2_w0m_(nomega,nprs) = nagelkerke_r2_ver0([tmp_Dw_ , tmp_mds_d_ ; tmp_Xw_ , tmp_mds_x_ ] , tmp_outcome_);
r2_wzm_(nomega,nprs) = nagelkerke_r2_ver0([tmp_Dw_ , tmp_Dz_ , tmp_mds_d_ ; tmp_Xw_ , tmp_Xz_ , tmp_mds_x_ ] , tmp_outcome_);
colormap(colormap_pm());
subplot(4,6,nomega);
imagesc(tmp_E_,[-1,1]*3/n_h^2); axis image; 
xlabel(sprintf('\\omega %0.2f',omega*180/pi));ylabel('prs');
set(gca,'XTick',[],'YTick',[]);
end;%for nomega=1:n_omega;
end;%for nprs=1:n_prs;
%%%%%%%%;

subplot(1,2,1);
r2_ratio_ = (r2_wzm_ - r2_00m_)./(r2_0zm_ - r2_00m_);
colormap(colormap_beach()); imagesc(r2_ratio_,[0.0,2.0]); colorbar;
xlabel('prs cutoff'); ylabel('omega');
set(gca,'XTick',1:n_prs,'XTickLabel',str_prs_,'YTick',1:n_omega,'YTickLabel',round(omega_*180/pi));
xtickangle(90);
title('nagelkerke pseudo-r^2 ratio (wzm-00m)/(0zm-00m)');
subplot(1,2,2);
r2_ratio_ = (r2_wzm_)./(r2_00m_);
colormap(colormap_beach()); imagesc(r2_ratio_,[0.0,2.00]); colorbar;
xlabel('prs cutoff'); ylabel('omega');
set(gca,'XTick',1:n_prs,'XTickLabel',str_prs_,'YTick',1:n_omega,'YTickLabel',round(omega_*180/pi));
xtickangle(90);
title('nagelkerke pseudo-r^2 ratio (wzm)/(00m)');
set(gcf,'Position',1+[0,0,1024*2,1024*1]);
lisa_arm1_r0.dir_out_s0_prs_jpg = sprintf('%s/dir_prs_jpg',lisa_arm1_r0.dir_out_s0_pca);
if (~exist(lisa_arm1_r0.dir_out_s0_prs_jpg,'dir')); disp(sprintf(' %% mkdir %s',lisa_arm1_r0.dir_out_s0_prs_jpg)); mkdir(lisa_arm1_r0.dir_out_s0_prs_jpg); end;
tmp_infix = sprintf('trn%d%s_tst%d%s_ni%d_ni%d',lisa_arm1_r0.cl_num,cl_num_arm1_string,lisa_arm2.cl_num,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_fname = sprintf('%s/%s_%s_prs_nagelkerk_ver0',lisa_arm1_r0.dir_out_s0_prs_jpg,lisa_arm1_r0.string_name_s0,tmp_infix);
print('-djpeg',sprintf('%s.jpg',tmp_fname));
print('-depsc',sprintf('%s.eps',tmp_fname));










