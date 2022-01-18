function [n_patient,fam_fam_,fam_sam_,fam_pat_,fam_mat_,fam_sex_,fam_dvx_,fam_dir_,fam_name_,fam_] = lisa_famext_ver0(fname);
% extract fields from fam.ext file. ;

fcheck(fname);
fid = fopen(fname); 
fam_ = textscan(fid,'%s %s %s %s %d %d %s','headerlines',0); 
fclose(fid);

n_patient = length(fam_{1});
assert(n_patient==wc_0(fname));

fam_fam_ = fam_{1}; %<-- family id. ;
fam_sam_ = fam_{2}; %<-- sample id. ;
fam_pat_ = fam_{3}; %<-- paternal id. ;
fam_mat_ = fam_{4}; %<-- maternal id. ;
fam_sex_ = fam_{5}; %<-- sex (1=male, 2=female, otherwise unknown). ;
fam_dvx_ = fam_{6}; %<-- case-control-status (0=unknown, 1=unaffected, 2=affected). ;
fam_dir_ = fam_{7}; %<-- directory and filename of original fam file. ;

if nargout>8;
fam_name_ = cell(n_patient,1); 
for np=1:n_patient; 
fam_name_{np} = sprintf('%s%s%s',fam_{1}{np},'&',fam_{2}{np}); 
end;%for np=1:n_patient;
end;%if nargout>8;
