%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run lisa_setprefix_ver2 first. ;
% run lisa_setnames_ver2 first. ;
% run lisa_xdropextract_ver2 first. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% load mds and fam and bim file. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fname_bim = sprintf('%s/%s_bim.ex2',dir__in,string_prefix); fcheck(fname_bim);
[n_snp,bim__id_,bim_al1_,bim_al2_,bim_alt_,bim_ent_,bim_frq_,bim_maf_,bim_name_,bim_] = lisa_bimex2_ver0(fname_bim);
fcheck(sprintf('%s/mds.mat',dir_trunk));
load(sprintf('%s/mds.mat',dir_trunk));
fname_fam = sprintf('%s/%s_fam.ex2',dir__in,string_prefix); fcheck(fname_fam);
[n_patient,fam_fam_,fam_sam_,fam_pat_,fam_mat_,fam_sex_,fam_dvx_,fam_dir_,fam_name_,fam_] = lisa_famext_ver0(fname_fam); n_patient_F = n_patient;
[~,ifam_,imds_] = intersect(fam_name_,mds_name_,'stable'); if (length(ifam_)<n_patient_F); disp('Warning! fam_name_ not a subset of mds_name_'); end; mds_sort_ = mds_(imds_,:);
