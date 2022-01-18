%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run lisa_setprefix_ver2 first. ;
% run lisa_setnames_ver2 first. ;
% run lisa_xdropextract_ver2 first. ;
% run lisa_loadmdsfambim_ver2 first. ;
% run lisa_mxcheck_ver2 first. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% determine which studies correspond to which patients. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
study_path_ = cell(n_study,1);
study_index_ = zeros(n_patient_F,1);
for ns=1:n_study;
study_path_{ns} = sprintf('/home/grouther/%s/%s.fam',study_trunk_{ns},study_name_{ns});
study_index_(find(strcmp(fam_{7},study_path_{ns}))) = ns;
end;%for ns=1:n_study;
