function ...
[ ...
 n_patient ...
,fam_fid_ ...
,fam_iid_ ...
,fam_yid_ ...
,fam_xid_ ...
,fam_sex_ ...
,fam_dvx_ ...
,fam_dir_ ...
,fam_fidandiid_ ...
,fam_ ...
  ] = ...
load_famext_ver1( ...
 fname ...
);
% extract fields from fam.ext file. ;

fcheck(fname);
fid = fopen(fname); 
fam_ = textscan(fid,'%s\t%s\t%s\t%s\t%d\t%d\t%s','headerlines',0); 
fclose(fid);

n_patient = length(fam_{1});
assert(n_patient==wc_0(fname));

fam_fid_ = fam_{1+0}; %<-- family id. ;
fam_iid_ = fam_{1+1}; %<-- sample id. ;
fam_xid_ = fam_{1+2}; %<-- paternal id. ;
fam_yid_ = fam_{1+3}; %<-- maternal id. ;
fam_sex_ = fam_{1+4}; %<-- sex (1=male, 2=female, otherwise unknown). ;
fam_dvx_ = fam_{1+5}; %<-- case-control-status (0=unknown, 1=unaffected, 2=affected). ;
fam_dir_ = fam_{1+6}; %<-- directory and filename of original fam file. ;

if nargout>8;
fam_fidandiid_ = cell(n_patient,1); 
for np=0:n_patient-1; 
fam_fidandiid_{1+np} = sprintf('%s%s%s',fam_fid_{1+np},'&',fam_iid_{1+np}); 
end;%for np=0:n_patient-1;
end;%if nargout>8;
