function lisa_prs_collect_ver10(prs_prefix,cl_num_arm1,cl_num_arm2,flag_reverse_use,flag_rerun,flag_replot);
% try: ;
%{

  flag_rerun = 0; flag_replot = 1;
  %prs_prefix_ = {'ori','756','469','319','MDD_exclude','SCZ_exclude','MDD_notsure','SCZ_notsure'};
  %prs_prefix_ = {'ori','756'};
  %prs_prefix_ = {'469','319'};
  prs_prefix_ = {'MDD_exclude','SCZ_exclude'};
  %prs_prefix_ = {'MDD_notsure','SCZ_notsure'};
  for nd=1:length(prs_prefix_); for cl_num_arm1=[4]; for cl_num_arm2=[1,2,3]; for flag_reverse_use=[0:1];
  lisa_prs_collect_ver10(prs_prefix_{nd},cl_num_arm1,cl_num_arm2,flag_reverse_use,flag_rerun,flag_replot);
  end;end;end;end;%for nd=1:length(prs_prefix_); for cl_num_arm1=[4]; for cl_num_arm2=[1,2,3]; for flag_reverse_use=[0:1];
  
  lisa_prs_collect_ver10('icuk_on_icuk',1,4,0,0,1);

  %}

if nargin<1; cl_num_arm1 = 4; end;
if nargin<2; cl_num_arm2 = 1; end;
if nargin<3; flag_reverse_use = 0; end;
% clear; prs_prefix = 'MDD_notsure'; flag_rerun=0; flag_replot=1; cl_num_arm1 = 4; cl_num_arm2 = 2; flag_reverse_use = 0;
% clear; prs_prefix = 'icuk_on_icuk'; flag_rerun=0; flag_replot=1; cl_num_arm1 = 1; cl_num_arm2 = 4; flag_reverse_use = 0;
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));
addpath('/data/rangan/dir_bcc/dir_code_022316');
addpath('/data/rangan/dir_bcc/dir_PGC_20180304/dir_m');
dir_trunk = '/data/rangan/dir_bcc/dir_PGC_20180304';
dir_bip32 = '/data/rangan/dir_bcc/dir_BIP32_loo_data';
dir_icuk = sprintf('%s/dir_icuk_scores',dir_bip32);
dir_icuk_ex = sprintf('%s/dir_icuk_ex_bc1',dir_bip32);
dir_icuk_on_icuk = sprintf('%s/dir_icuk_on_icuk_loo_proriles',dir_bip32);
dir_SCZ_MDD = sprintf('%s/dir_SCZ_MDD_PRS_profile_scores',dir_bip32);
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
B_MLT = 34; n_mds = 20; n_scramble = 1; n_shuffle = 64;
flag_dex_vs_lak = {'dex'};
if (strcmp(flag_dex_vs_lak{1},'dex')); gamma = 0.004; mc_string = ''; end; 
if (strcmp(flag_dex_vs_lak{1},'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
flag_reverse = flag_reverse_use;
cl_num = cl_num_arm1;
n_maf = 5; n_cov = 2;
%%%%%%%%;
lisa_setprefix_ver2 ;
lisa_setnames_ver2 ;
lisa_xdropextract_ver2 ;
lisa_loadmdsfambim_ver2 ;
lisa_mxcheck_ver2 ;
lisa_studyindex_ver2 ;
%%%%%%%%;
bim_arm1_ = bim_;
fam_arm1_ = fam_;
fam_name_arm1_ = fam_name_;
mds_sort_arm1_ = mds_sort_;
string_prefix_arm1 = string_prefix;
%string_name_arm1 = string_name;
string_name_r1_arm1 = string_name;
string_name_r0_arm1 = string_name_r0;
%dir_out_s0000_arm1 = dir_out_s0000;
dir_out_s0000_r1_arm1 = dir_out_s0000;
dir_out_s0000_r0_arm1 = dir_out_s0000_r0;
%dir_out_trace_arm1 = dir_out_trace;
dir_out_trace_r1_arm1 = dir_out_trace;
dir_out_trace_r0_arm1 = dir_out_trace_r0;
%%%%%%%%;
clear bim_ fam_ mds_sort_ mds_;
%%%%%%%%;
verbose = 1;
B_MLT = 34; n_mds = 20; n_scramble = 1; n_shuffle = 64;
flag_dex_vs_lak = {'dex'};
if (strcmp(flag_dex_vs_lak{1},'dex')); gamma = 0.004; mc_string = ''; end; 
if (strcmp(flag_dex_vs_lak{1},'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
flag_reverse = flag_reverse_use;
cl_num = cl_num_arm2;
n_maf = 5; n_cov = 2;
%%%%%%%%;
lisa_setprefix_ver2 ;
lisa_setnames_ver2 ;
lisa_xdropextract_ver2 ;
lisa_loadmdsfambim_ver2 ;
lisa_mxcheck_ver2 ;
lisa_studyindex_ver2 ;
%%%%%%%%;
study_trunk_arm2_ = study_trunk_;
study_name_arm2_ = study_name_;
n_study_arm2 = n_study;
bim_arm2_ = bim_;
fam_arm2_ = fam_;
fam_name_arm2_ = fam_name_;
mds_sort_arm2_ = mds_sort_;
string_prefix_arm2 = string_prefix;
%string_name_arm2 = string_name;
string_name_r1_arm2 = string_name;
string_name_r0_arm2 = string_name_r0;
%dir_out_s0000_arm2 = dir_out_s0000;
dir_out_s0000_r1_arm2 = dir_out_s0000;
dir_out_s0000_r0_arm2 = dir_out_s0000_r0;
%dir_out_trace_arm2 = dir_out_trace;
dir_out_trace_r1_arm2 = dir_out_trace;
dir_out_trace_r0_arm2 = dir_out_trace_r0;
%%%%%%%%;
clear bim_ fam_ mds_sort_ mds_;
%%%%%%%%;

%%%%%%%%%%%%%%%%;
% copy fam file. ;
%%%%%%%%%%%%%%%%;
ni=1;
fam_fam_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_iid_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_pat_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_mat_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_sex_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_phe_arm2_ = fam_arm2_{ni}; ni=ni+1;
fam_dir_arm2_ = fam_arm2_{ni}; ni=ni+1;

%%%%%%%%%%%%%%%%;
% extract mc;
%%%%%%%%%%%%%%%%;
if (flag_reverse_use==0);
fname = sprintf('%s/%s_mr_A_full.b16',dir__in,string_prefix); fcheck(fname);
mr_D_ = tutorial_binary_uncompress(fname)>0;
fname = sprintf('%s/%s_mr_Z_full.b16',dir__in,string_prefix); fcheck(fname);
mr_X_ = tutorial_binary_uncompress(fname)>0;
end;%if (flag_reverse_use==0);
if (flag_reverse_use==1);
fname = sprintf('%s/%s_mr_Z_full.b16',dir__in,string_prefix); fcheck(fname);
mr_D_ = tutorial_binary_uncompress(fname)>0;
fname = sprintf('%s/%s_mr_A_full.b16',dir__in,string_prefix); fcheck(fname);
mr_X_ = tutorial_binary_uncompress(fname)>0;
end;%if (flag_reverse_use==1);
%%%%%%%%%%%%%%%%;
if (flag_reverse_use==0);
disp(sprintf('mr_D error: %0.16f',norm(cast(mr_D_,'double')-cast(fam_phe_arm2_==2,'double'))));
disp(sprintf('mr_X error: %0.16f',norm(cast(mr_X_,'double')-cast(fam_phe_arm2_==1,'double'))));
end;%if (flag_reverse_use==0);
if (flag_reverse_use==1);
disp(sprintf('mr_D error: %0.16f',norm(cast(mr_D_,'double')-cast(fam_phe_arm2_==1,'double'))));
disp(sprintf('mr_X error: %0.16f',norm(cast(mr_X_,'double')-cast(fam_phe_arm2_==2,'double'))));
end;%if (flag_reverse_use==1);

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

%%%%%%%%;
%flag_rerun = 0;
dir_out_s0000_pca_r0_arm1 = sprintf('%s/dir_pca',dir_out_s0000_r0_arm1);
dir_out_s0000_pca_r1_arm1 = sprintf('%s/dir_pca',dir_out_s0000_r1_arm1);
fname_mat = sprintf('%s/prs_collect_ver10_%s_trn%d_tst%d.mat',dir_out_s0000_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm2);
if (~flag_rerun & exist(fname_mat,'file')); disp(sprintf(' %% %s found, loading',fname_mat)); load(fname_mat); end;
if ( flag_rerun | ~exist(fname_mat,'file')); disp(sprintf(' %% %s not found, creating',fname_mat)); 
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
auc00_prs_ = zeros(n_prs,1);
auc02_prs_ = zeros(n_prs,1);
auc06_prs_ = zeros(n_prs,1);
auc19_prs_ = zeros(n_prs,1);
auc20_prs_ = zeros(n_prs,1);
logp_auc00_prs_ = zeros(n_prs,1);
logp_auc02_prs_ = zeros(n_prs,1);
logp_auc06_prs_ = zeros(n_prs,1);
logp_auc19_prs_ = zeros(n_prs,1);
logp_auc20_prs_ = zeros(n_prs,1);
%%%%%%%%;
for nprs=1:n_prs;
%%%%%%%%;
n_patient_prs_(nprs)=0;
for nstudy_arm2=1:n_study_arm2;
tmp_study_name = study_name_arm2_{nstudy_arm2};
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
end;%for nstudy_arm2=1:n_study_arm2;
prs_name_{nprs} = cell(n_patient_prs_(nprs),1); 
for np=1:n_patient_prs_(nprs); prs_name_{nprs}{np} = sprintf('%s%s%s',prs_fam_{nprs}{np},'&',prs_iid_{nprs}{np}); end;%for np=1:n_patient_prs_(nprs);
if (flag_reverse_use==1); prs_sco_{nprs} = -prs_sco_{nprs}; end;
%%%%%%%%%%%%%%%%;
disp(sprintf(' %% nprs %d: %s: intersection %d',nprs,tmp_fname,length(intersect(prs_name_{nprs},fam_name_arm2_))));
%%%%%%%%%%%%%%%%;
[~,ifam_,iprs_{nprs}] = intersect(fam_name_arm2_,prs_name_{nprs},'stable'); if (length(ifam_)<length(fam_name_arm2_)); disp('Warning! fam_name_ not a subset of prs_name_{nprs}'); end; 
tmp_f_ = cast(fam_phe_arm2_,'double'); tmp_p_ = transpose(cast(prs_phe_{nprs}(iprs_{nprs}),'double'));
disp(sprintf('phe error: %0.16f',norm(tmp_f_-tmp_p_,'fro')));
%%%%%%%%%%%%%%%%;
% build libraries for cauc p-value ;
%%%%%%%%%%%%%%%%;
find_x_ = find(mr_X_ & transpose(prs_inc_{nprs})); find_d_ = find(mr_D_ & transpose(prs_inc_{nprs}));
[cauc02_avg_,cauc02_min_] = p_cauc_0(mds_sort_arm2_(find_x_,1: 2),mds_sort_arm2_(find_d_,1: 2),zeros(length(find_x_),1),zeros(length(find_d_),1),1024*8);
[cauc06_avg_,cauc06_min_] = p_cauc_0(mds_sort_arm2_(find_x_,1: 6),mds_sort_arm2_(find_d_,1: 6),zeros(length(find_x_),1),zeros(length(find_d_),1),1024*8);
[cauc19_avg_,cauc19_min_] = p_cauc_0(mds_sort_arm2_(find_x_,[1:6,19]),mds_sort_arm2_(find_d_,[1:6,19]),zeros(length(find_x_),1),zeros(length(find_d_),1),1024*8);
[cauc20_avg_,cauc20_min_] = p_cauc_0(mds_sort_arm2_(find_x_,1:20),mds_sort_arm2_(find_d_,1:20),zeros(length(find_x_),1),zeros(length(find_d_),1),1024*8);
%%%%%%%%%%%%%%%%;
% calculate auc and p-value for prs;
%%%%%%%%%%%%%%%%;
auc00_prs_(nprs) = auc_0(prs_sco_{nprs}(iprs_{nprs}(find_x_)),prs_sco_{nprs}(iprs_{nprs}(find_d_)));
logp_auc00_prs_(nprs) = logp_auc_0(auc00_prs_(nprs),min(length(find_x_),length(find_d_)));
auc02_prs_(nprs) = cauc_0(prs_sco_{nprs}(iprs_{nprs}(find_x_)),prs_sco_{nprs}(iprs_{nprs}(find_d_)),mds_sort_arm2_(find_x_,1: 2),mds_sort_arm2_(find_d_,1: 2));
logp_auc02_prs_(nprs) = log(max(1,length(find(cauc02_avg_>auc02_prs_(nprs)))/length(cauc02_avg_)));
auc06_prs_(nprs) = cauc_0(prs_sco_{nprs}(iprs_{nprs}(find_x_)),prs_sco_{nprs}(iprs_{nprs}(find_d_)),mds_sort_arm2_(find_x_,1: 6),mds_sort_arm2_(find_d_,1: 6));
logp_auc06_prs_(nprs) = log(max(1,length(find(cauc06_avg_>auc06_prs_(nprs)))/length(cauc06_avg_)));
auc19_prs_(nprs) = cauc_0(prs_sco_{nprs}(iprs_{nprs}(find_x_)),prs_sco_{nprs}(iprs_{nprs}(find_d_)),mds_sort_arm2_(find_x_,[1:6,19]),mds_sort_arm2_(find_d_,[1:6,19]));
logp_auc19_prs_(nprs) = log(max(1,length(find(cauc19_avg_>auc19_prs_(nprs)))/length(cauc19_avg_)));
auc20_prs_(nprs) = cauc_0(prs_sco_{nprs}(iprs_{nprs}(find_x_)),prs_sco_{nprs}(iprs_{nprs}(find_d_)),mds_sort_arm2_(find_x_,1:20),mds_sort_arm2_(find_d_,1:20));
logp_auc20_prs_(nprs) = log(max(1,length(find(cauc20_avg_>auc20_prs_(nprs)))/length(cauc20_avg_)));
%%%%%%%%;
end;%for nprs=1:n_prs;

%%%%%%%%%%%%%%%%;
% calculate auc and p-value for arm2;
%%%%%%%%%%%%%%%%;
prs_min_ = zeros(n_prs,1); prs_max_ = zeros(n_prs,1);
for nprs=1:n_prs; prs_min_(nprs) = min(prs_sco_{nprs}); prs_max_(nprs) = max(prs_sco_{nprs}); end;%for nprs=1:n_prs;
n_r0_t = 21; n_r1_t = 21;
if (flag_reverse_use==0);
if (cl_num_arm1==1); ni_r0_arm2_ = [350]; ni_r1_arm2_ = [100]; end;
if (cl_num_arm1==2); ni_r0_arm2_ = [350]; ni_r1_arm2_ = [350]; end;
if (cl_num_arm1==3); ni_r0_arm2_ = [400]; ni_r1_arm2_ = [425]; end;
if (cl_num_arm1==4); ni_r0_arm2_ = [275,350,450]; ni_r1_arm2_ = [275,350]; end;
end;%if (flag_reverse_use==0);
if (flag_reverse_use==1);
if (cl_num_arm1==1); ni_r0_arm2_ = [400]; ni_r1_arm2_ = [150]; end;
if (cl_num_arm1==2); ni_r0_arm2_ = [450]; ni_r1_arm2_ = [425]; end;
if (cl_num_arm1==3); ni_r0_arm2_ = [450]; ni_r1_arm2_ = [500]; end;
if (cl_num_arm1==4); ni_r0_arm2_ = [300]; ni_r1_arm2_ = [350]; end;
end;%if (flag_reverse_use==1);
auc19_prs_t_ = zeros(n_r0_t,n_r1_t,n_prs,length(ni_r0_arm2_),length(ni_r1_arm2_));
p_x_rem_ = zeros(n_r0_t,n_r1_t,length(ni_r0_arm2_),length(ni_r1_arm2_));
p_d_rem_ = zeros(n_r0_t,n_r1_t,length(ni_r0_arm2_),length(ni_r1_arm2_));
for nir0_arm2_ = 1:length(ni_r0_arm2_); 
ni_r0_arm2 = ni_r0_arm2_(nir0_arm2_);
tmp_fname = sprintf('%s/pca_ni%d_tst%d_k2_B44_V_arm2_.mda',dir_out_s0000_pca_r0_arm1,ni_r0_arm2,cl_num_arm2); fcheck(tmp_fname); V_r0_arm2_ = mda_read_r8(tmp_fname);
tmp_fname = sprintf('%s/pca_proj_arm2_ni%d_tst%d_k2_B44_AnV_.mda',dir_out_s0000_pca_r0_arm1,ni_r0_arm2,cl_num_arm2); fcheck(tmp_fname); DnV_r0_arm2_ = mda_read_r8(tmp_fname);
tmp_fname = sprintf('%s/pca_proj_arm2_ni%d_tst%d_k2_B44_ZnV_.mda',dir_out_s0000_pca_r0_arm1,ni_r0_arm2,cl_num_arm2); fcheck(tmp_fname); XnV_r0_arm2_ = mda_read_r8(tmp_fname);
tmp_mds_x_ = mds_sort_arm2_(find_x_,[1:6,19]); tmp_mds_d_ = mds_sort_arm2_(find_d_,[1:6,19]);
tmp_r0_auc19 = cauc_0(XnV_r0_arm2_(find_x_,1),DnV_r0_arm2_(find_d_,1),tmp_mds_x_,tmp_mds_d_);
if (tmp_r0_auc19<0.5); disp(sprintf(' %% flipping nir0_arm2_ %d',nir0_arm2_)); DnV_r0_arm2_ = -DnV_r0_arm2_; XnV_r0_arm2_ = -XnV_r0_arm2_; end;
%tmp_r0_min = prctile([DnV_r0_arm2_(find_d_,1);XnV_r0_arm2_(find_x_,1)], 0);
%tmp_r0_max = prctile([DnV_r0_arm2_(find_d_,1);XnV_r0_arm2_(find_x_,1)], 95);
tmp_r0_min = prctile(DnV_r0_arm2_(find_d_,1), 0);
tmp_r0_max = prctile(DnV_r0_arm2_(find_d_,1), 95);
tmp_r0_t_ = prctile(DnV_r0_arm2_(find_d_,1), linspace(0,100,n_r0_t) );
for nir1_arm2_ = 1:length(ni_r1_arm2_);
ni_r1_arm2 = ni_r1_arm2_(nir1_arm2_);
tmp_fname = sprintf('%s/pca_ni%d_tst%d_k2_B44_V_arm2_.mda',dir_out_s0000_pca_r1_arm1,ni_r1_arm2,cl_num_arm2); fcheck(tmp_fname); V_r1_arm2_ = mda_read_r8(tmp_fname);
tmp_fname = sprintf('%s/pca_proj_arm2_ni%d_tst%d_k2_B44_AnV_.mda',dir_out_s0000_pca_r1_arm1,ni_r1_arm2,cl_num_arm2); fcheck(tmp_fname); DnV_r1_arm2_ = mda_read_r8(tmp_fname);
tmp_fname = sprintf('%s/pca_proj_arm2_ni%d_tst%d_k2_B44_ZnV_.mda',dir_out_s0000_pca_r1_arm1,ni_r1_arm2,cl_num_arm2); fcheck(tmp_fname); XnV_r1_arm2_ = mda_read_r8(tmp_fname);
tmp_mds_x_ = mds_sort_arm2_(find_x_,[1:6,19]); tmp_mds_d_ = mds_sort_arm2_(find_d_,[1:6,19]);
tmp_r1_auc19 = cauc_0(XnV_r1_arm2_(find_x_,1),DnV_r1_arm2_(find_d_,1),tmp_mds_x_,tmp_mds_d_);
if(tmp_r1_auc19<0.5); disp(sprintf(' %% flipping nir1_arm2_ %d',nir1_arm2_)); DnV_r1_arm2_ = -DnV_r1_arm2_; XnV_r1_arm2_ = -XnV_r1_arm2_; end;
%tmp_r1_min = prctile([DnV_r1_arm2_(find_d_,1);XnV_r1_arm2_(find_x_,1)], 0);
%tmp_r1_max = prctile([DnV_r1_arm2_(find_d_,1);XnV_r1_arm2_(find_x_,1)], 95);
tmp_r1_min = prctile(DnV_r1_arm2_(find_d_,1), 0);
tmp_r1_max = prctile(DnV_r1_arm2_(find_d_,1), 95);
tmp_r1_t_ = prctile(DnV_r1_arm2_(find_d_,1), linspace(0,100,n_r1_t) );
%%%%%%%%;
for nr0_t=1:n_r0_t; for nr1_t=1:n_r1_t;
tmp_r0_t = tmp_r0_t_(nr0_t); tmp_r1_t = tmp_r1_t_(nr1_t); 
find_xt_ = find((XnV_r0_arm2_(find_x_,1)>=tmp_r0_t) & (XnV_r1_arm2_(find_x_,1)>=tmp_r1_t)); 
find_dt_ = find((DnV_r0_arm2_(find_d_,1)>=tmp_r0_t) & (DnV_r1_arm2_(find_d_,1)>=tmp_r1_t));
p_x_rem_(nr0_t,nr1_t,nir0_arm2_,nir1_arm2_) = length(find_xt_)/length(find_x_);
p_d_rem_(nr0_t,nr1_t,nir0_arm2_,nir1_arm2_) = length(find_dt_)/length(find_d_);
for nprs=1:n_prs;
tmp_prs_x_ = prs_sco_{nprs}(iprs_{nprs}(find_x_(find_xt_))); tmp_prs_d_ = prs_sco_{nprs}(iprs_{nprs}(find_d_(find_dt_)));
tmp_mds_x_ = mds_sort_arm2_(find_x_(find_xt_),[1:6,19]); tmp_mds_d_ = mds_sort_arm2_(find_d_(find_dt_),[1:6,19]);
if (length(find_xt_)==0 & length(find_dt_)==0); auc19_prs_t_(nr0_t,nr1_t,nprs,nir0_arm2_,nir1_arm2_) = 0.5; end;
if (length(find_xt_)==0 & length(find_dt_)> 0); auc19_prs_t_(nr0_t,nr1_t,nprs,nir0_arm2_,nir1_arm2_) = 1.0; end;
if (length(find_xt_)> 0 & length(find_dt_)==0); auc19_prs_t_(nr0_t,nr1_t,nprs,nir0_arm2_,nir1_arm2_) = 0.0; end;
if (length(find_xt_)> 0 & length(find_dt_)> 0); auc19_prs_t_(nr0_t,nr1_t,nprs,nir0_arm2_,nir1_arm2_) = cauc_0(tmp_prs_x_,tmp_prs_d_,tmp_mds_x_,tmp_mds_d_); end;
end;%for nprs=1:n_prs;
end;end;%for nr0_t=1:n_r0_t;for nr1_t=1:n_r1_t;
%%%%%%%%;
end;end;%for nir0_arm2_ = 1:length(ni_r0_arm2_); for nir1_arm2_ = 1:length(ni_r1_arm2_);
%%%%%%%%;
lp_auc00_ = zeros(n_r0_t,n_r1_t,n_prs,length(ni_r0_arm2_),length(ni_r1_arm2_));
for nir0_arm2_ = 1:length(ni_r0_arm2_); 
ni_r0_arm2 = ni_r0_arm2_(nir0_arm2_);
for nir1_arm2_ = 1:length(ni_r1_arm2_);
ni_r1_arm2 = ni_r1_arm2_(nir1_arm2_);
for nr0_t=1:n_r0_t; for nr1_t=1:n_r1_t;
tmp_r0_t = tmp_r0_t_(nr0_t); tmp_r1_t = tmp_r1_t_(nr1_t); 
for nprs=1:n_prs;
tmp_n = min(p_x_rem_(nr0_t,nr1_t,nir0_arm2_,nir1_arm2_)*length(find_x_),p_d_rem_(nr0_t,nr1_t,nir0_arm2_,nir1_arm2_)*length(find_d_));
tmp_a = auc19_prs_t_(nr0_t,nr1_t,nprs,nir0_arm2_,nir1_arm2_);
lp_auc00_(nr0_t,nr1_t,nprs,nir0_arm2_,nir1_arm2_) = logp_auc_0(tmp_a,tmp_n);
end;%for nprs=1:n_prs;
end;end;%for nr0_t=1:n_r0_t;for nr1_t=1:n_r1_t;
end;end;%for nir0_arm2_ = 1:length(ni_r0_arm2_); for nir1_arm2_ = 1:length(ni_r1_arm2_);
disp(sprintf(' %% saving %s',fname_mat)); save(fname_mat);
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
% plot output. ;
%%%%%%%%;
%flag_replot = 0;
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
if flag_reverse_use==0; str_reverse = 'D'; end;
if flag_reverse_use==1; str_reverse = 'X'; end;
for nir0_arm2_ = 1:length(ni_r0_arm2_); 
ni_r0_arm2 = ni_r0_arm2_(nir0_arm2_);
for nir1_arm2_ = 1:length(ni_r1_arm2_);
ni_r1_arm2 = ni_r1_arm2_(nir1_arm2_);
fname_fig = sprintf('%s/prs_collect_ver10_%s_trn%d_tst%d_ni%.3d_ni%.3d.jpg',dir_out_s0000_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm2,ni_r0_arm2,ni_r1_arm2);
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
imagesc(squeeze(auc19_prs_t_(:,:,nprs,nir0_arm2_,nir1_arm2_)),cauclim_);
tmp_tick_r1_ = round(linspace(1,n_r1_t,n_tick)); xlabel('r1 threshold (p_rem)'); set(gca,'XTick',tmp_tick_r1_,'XTickLabel',0.01*round(100*p_d_rem_(1,tmp_tick_r1_,nir0_arm2_,nir1_arm2_)));
tmp_tick_r0_ = round(linspace(1,n_r0_t,n_tick)); ylabel('r0 threshold (p_rem)'); set(gca,'YTick',tmp_tick_r0_,'YTickLabel',0.01*round(100*p_d_rem_(tmp_tick_r0_,1,nir0_arm2_,nir1_arm2_)));
xtickangle(90);
%title(sprintf('%s S%s cauc [%0.2f %0.2f]',str_reverse,str_prs_{nprs},cauclim_));
title(sprintf('%s snp-p<%s cauc [%0.2f %0.2f]',str_reverse,str_prs_cutoff_{nprs},cauclim_));
%%%%%%%%;
subplot(4,length(nprs_sub_)/2,nprs_sub + 0*length(nprs_sub_));
hold on;
plot(1:n_r0_t,squeeze(auc19_prs_t_(1,:,nprs,nir0_arm2_,nir1_arm2_)),'b.-','MarkerSize',25,'LineWidth',2); 
plot(1:n_r1_t,squeeze(auc19_prs_t_(:,1,nprs,nir0_arm2_,nir1_arm2_)),'r.-','MarkerSize',25,'LineWidth',2); 
xlim([1,n_r0_t]); ylim(cauclim_); legend('bc2','bc1');
tmp_tick_ = round(linspace(1,n_r0_t,n_tick)); %tmp_tick_ = tmp_tick_(end:-1:1);
tmp_tick_val_ = tmp_tick_; %tmp_tick_val_ = p_d_rem_(tmp_tick_,1,nir0_arm2_,nir1_arm2_);
tmp_tick_label_ = num2str(0.01*round(100*p_d_rem_(tmp_tick_,1,nir0_arm2_,nir1_arm2_)),'%0.2f');
set(gca,'XTick',tmp_tick_val_,'XTickLabel',tmp_tick_label_);
xlabel('r threshold (p_rem)'); 
xtickangle(90);
ylabel('cauc');
%title(sprintf('%s S%s cauc',str_reverse,str_prs_{nprs}));
title(sprintf('%s snp-p<%s cauc',str_reverse,str_prs_cutoff_{nprs}));
%%%%%%%%;
%subplot(3,length(nprs_sub_),nprs_sub + 1*length(nprs_sub_));
%imagesc(squeeze(p_d_rem_(:,:,nir0_arm2_,nir1_arm2_)),[0.0,1.0]);
%tmp_tick_r1_ = round(linspace(1,n_r1_t,n_tick)); xlabel('r1 threshold (p_rem)'); set(gca,'XTick',tmp_tick_r1_,'XTickLabel',0.01*round(100*p_d_rem_(1,tmp_tick_r1_,nir0_arm2_,nir1_arm2_)));
%tmp_tick_r0_ = round(linspace(1,n_r0_t,n_tick)); ylabel('r0 threshold (p_rem)'); set(gca,'YTick',tmp_tick_r0_,'YTickLabel',0.01*round(100*p_d_rem_(tmp_tick_r0_,1,nir0_arm2_,nir1_arm2_)));
%xtickangle(90);
%title(sprintf('%s S%s fraction remaining',str_reverse,str_prs_cutoff_{nprs}));
%%%%%%%%;
%subplot(3,length(nprs_sub_),nprs_sub + 2*length(nprs_sub_));
%imagesc(-(1/log(10))*squeeze(lp_auc00_(:,:,nprs,nir0_arm2_,nir1_arm2_)),[0,10]);
%tmp_tick_r1_ = round(linspace(1,n_r1_t,n_tick)); xlabel('r1 threshold (p_rem)'); set(gca,'XTick',tmp_tick_r1_,'XTickLabel',0.01*round(100*p_d_rem_(1,tmp_tick_r1_,nir0_arm2_,nir1_arm2_)));
%tmp_tick_r0_ = round(linspace(1,n_r0_t,n_tick)); ylabel('r0 threshold (p_rem)'); set(gca,'YTick',tmp_tick_r0_,'YTickLabel',0.01*round(100*p_d_rem_(tmp_tick_r0_,1,nir0_arm2_,nir1_arm2_)));
%xtickangle(90);
%title(sprintf('%s S%s -log10(pvalue)',str_reverse,str_prs_cutoff_{nprs}));
%%%%%%%%;
end;%for nprs=1:n_prs;
set(gcf,'Position',1+[0,0,1024*2,1024]);
print('-djpeg',fname_fig);
disp(sprintf(' %% writing to %s',fname_fig));
end;%if flag_disp;
end;%if (~exist(fname_fig,'file'));
end;end;%for nir0_arm2_ = 1:length(ni_r0_arm2_); for nir1_arm2_ = 1:length(ni_r1_arm2_);












