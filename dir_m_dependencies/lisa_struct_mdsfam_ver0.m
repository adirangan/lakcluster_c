function l = lisa_struct_mdsfam_ver0(l);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% load mds and fam file. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%fname_bim = sprintf('%s/%s_bim.ex2',dir__in,string_prefix); fcheck(fname_bim);
%[n_snp,bim__id_,bim_al1_,bim_al2_,bim_alt_,bim_ent_,bim_frq_,bim_maf_,bim_name_,bim_] = lisa_bimex2_ver0(fname_bim);
l.mds_str = sprintf('m%dr%d',length(l.mds_used_),l.mds_repl);
l.fname_mds_kappa_squared = sprintf('%s/%s_T_%s_kappa.txt',l.dir__in,l.string_prefix,l.mds_str);
l.T_n_cols = 1+length(l.mds_used_)*l.mds_repl;
fcheck(sprintf('%s/mds.mat',l.dir_trunk));
load(sprintf('%s/mds.mat',l.dir_trunk),'mds_name_','mds_');
l.fname_fam = sprintf('%s/%s_fam.ex2',l.dir__in,l.string_prefix); fcheck(l.fname_fam);
[l.n_patient,l.fam_fam_,l.fam_sam_,l.fam_pat_,l.fam_mat_,l.fam_sex_,l.fam_dvx_,l.fam_dir_,l.fam_name_,l.fam_] = lisa_famext_ver0(l.fname_fam);
[~,ifam_,imds_] = intersect(l.fam_name_,mds_name_,'stable'); if (length(ifam_)<l.n_patient); disp('Warning! fam_name_ not a subset of mds_name_'); end; l.mds_sort_ = mds_(imds_,:);
