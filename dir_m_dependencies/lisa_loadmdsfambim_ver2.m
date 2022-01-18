%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run lisa_setprefix_ver2 first. ;
% run lisa_setnames_ver2 first. ;
% run lisa_xdropextract_ver2 first. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% load mds and fam and bim file. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fcheck(sprintf('%s/%s_bim.ext',dir__in,string_prefix));
fid = fopen(sprintf('%s/%s_bim.ext',dir__in,string_prefix)); bim_ = textscan(fid,'%d\t%s\t%d\t%d\t%c\t%c\t%s\t%f\t%f\t%f\t%f\n','headerlines',0); fclose(fid);
n_snp = length(bim_{1});
fcheck(sprintf('%s/mds.mat',dir_trunk));
load(sprintf('%s/mds.mat',dir_trunk));
fcheck(sprintf('%s/%s_fam.ext',dir__in,string_prefix));
fid = fopen(sprintf('%s/%s_fam.ext',dir__in,string_prefix)); fam_ = textscan(fid,'%s %s %s %s %d %d %s','headerlines',0); fclose(fid);
n_patient_F = length(fam_{1}); fam_name_ = cell(n_patient_F,1); for np=1:n_patient_F; fam_name_{np} = sprintf('%s%s%s',fam_{1}{np},'&',fam_{2}{np}); end;%for np=1:n_patient_F;
[~,ifam_,imds_] = intersect(fam_name_,mds_name_,'stable'); if (length(ifam_)<n_patient_F); disp('Warning! fam_name_ not a subset of mds_name_'); end; mds_sort_ = mds_(imds_,:);
