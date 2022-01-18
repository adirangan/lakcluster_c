function test_loader_helper_wrap_2(dir_trunk,n_iteration,n_Sample_Label,str_Sample_Label_,C_VariableName_LC_,E_rank_,I_rank_,xeta_,n_rank);
u_Sample_Label_ = unique(str_Sample_Label_);
n_CCOV = size(xeta_.C_rank_,2);
n_I_GENE = xeta_.n_I_GENE;
n_E_GENE = xeta_.n_E_GENE;

for nrank=1:n_rank;

disp(sprintf(' %% nrank %d/%d',nrank,n_rank));

n_zeta_E = min(nrank,xeta_.n_zeta_E);
n_beta_E = min(nrank,xeta_.n_beta_E);
n_zeta_I = min(nrank,xeta_.n_zeta_I);
n_beta_I = min(nrank,xeta_.n_beta_I);
n_zeta_EI = min(nrank,xeta_.n_zeta_EI);
n_beta_EI = min(nrank,xeta_.n_beta_EI);

%%%%%%%%;
str_prefix = sprintf('E_%s_r%d',xeta_.infix,nrank);
fname_pre = sprintf('%s/dir_mat/dir_AB/AB_%s',dir_trunk,str_prefix);
fname_mat = sprintf('%s.mat',fname_pre);
fname_tmp = sprintf('%s.tmp',fname_pre);
if ( exist(fname_mat,'file') |  exist(fname_tmp,'file'));
disp(sprintf(' %% %s found, not creating',fname_pre));
end;%if ( exist(fname_mat,'file') |  exist(fname_tmp,'file'));
if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
disp(sprintf(' %% %s not found, creating',fname_pre));
save(fname_tmp,'fname_mat');
try; 
test_loader_helper_2(dir_trunk,str_prefix,n_iteration,n_Sample_Label,n_CCOV,n_E_GENE,str_Sample_Label_,xeta_.C_rank_,E_rank_,xeta_.zeta_E_un_,xeta_.zeta_E_vn_,n_zeta_E,xeta_.beta_E_un_,xeta_.beta_E_vn_,n_beta_E);
catch; disp(sprintf(' %% WARNING: error creating %s',fname_mat)); end;%try;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
%%%%%%%%;

%%%%%%%%;
str_prefix = sprintf('I_%s_r%d',xeta_.infix,nrank);
fname_pre = sprintf('%s/dir_mat/dir_AB/AB_%s',dir_trunk,str_prefix);
fname_mat = sprintf('%s.mat',fname_pre);
fname_tmp = sprintf('%s.tmp',fname_pre);
if ( exist(fname_mat,'file') |  exist(fname_tmp,'file'));
disp(sprintf(' %% %s found, not creating',fname_pre));
end;%if ( exist(fname_mat,'file') |  exist(fname_tmp,'file'));
if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
disp(sprintf(' %% %s not found, creating',fname_pre));
save(fname_tmp,'fname_mat');
try; 
test_loader_helper_2(dir_trunk,str_prefix,n_iteration,n_Sample_Label,n_CCOV,n_I_GENE,str_Sample_Label_,xeta_.C_rank_,I_rank_,xeta_.zeta_I_un_,xeta_.zeta_I_vn_,n_zeta_I,xeta_.beta_I_un_,xeta_.beta_I_vn_,n_beta_I);
catch; disp(sprintf(' %% WARNING: error creating %s',fname_mat)); end;%try;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
%%%%%%%%;

%%%%%%%%;
str_prefix = sprintf('EI_%s_r%d',xeta_.infix,nrank);
fname_pre = sprintf('%s/dir_mat/dir_AB/AB_%s',dir_trunk,str_prefix);
fname_mat = sprintf('%s.mat',fname_pre);
fname_tmp = sprintf('%s.tmp',fname_pre);
if ( exist(fname_mat,'file') |  exist(fname_tmp,'file'));
disp(sprintf(' %% %s found, not creating',fname_pre));
end;%if ( exist(fname_mat,'file') |  exist(fname_tmp,'file'));
if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
disp(sprintf(' %% %s not found, creating',fname_pre));
save(fname_tmp,'fname_mat');
try; 
test_loader_helper_2(dir_trunk,str_prefix,n_iteration,n_Sample_Label,n_CCOV,n_E_GENE + n_I_GENE,str_Sample_Label_,xeta_.C_rank_,[E_rank_ , I_rank_],xeta_.zeta_EI_un_,xeta_.zeta_EI_vn_,n_zeta_EI,xeta_.beta_EI_un_,xeta_.beta_EI_vn_,n_beta_EI);
catch; disp(sprintf(' %% WARNING: error creating %s',fname_mat)); end;%try;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file') & ~exist(fname_tmp,'file'));
%%%%%%%%;

end;%for nrank=1:n_rank;
