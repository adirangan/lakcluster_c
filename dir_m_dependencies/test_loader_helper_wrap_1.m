function test_loader_helper_wrap_1(dir_trunk,n_iteration,n_CLabel,str_CLabel_,C_VariableName_LC_,E_rank_,I_rank_,xeta_,n_rank);
u_CLabel_ = unique(str_CLabel_);
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

str_prefix = sprintf('E_%s_r%d',xeta_.infix,nrank);
%%%%%%%%;
[AB_E_] = test_loader_helper_1(n_iteration,n_CLabel,n_CCOV,n_E_GENE,str_CLabel_,xeta_.C_rank_,E_rank_,xeta_.zeta_E_un_,xeta_.zeta_E_vn_,n_zeta_E,xeta_.beta_E_un_,xeta_.beta_E_vn_,n_beta_E);
test_loader_plot_absZ_avg__0(AB_E_,C_VariableName_LC_,dir_trunk,str_prefix);
test_loader_plot_C_absZ_0(AB_E_,C_VariableName_LC_,dir_trunk,str_prefix);
test_loader_plot_Z_absZ_0(AB_E_,dir_trunk,str_prefix);
test_loader_helper_tsv_0(AB_E_,dir_trunk,str_prefix);
%%%%%%%%;
for nl=1:3;
if nl==1; nA=find(strcmp(u_CLabel_,'28')); nB=find(strcmp(u_CLabel_,'33')); end;
if nl==2; nA=find(strcmp(u_CLabel_,'22')); nB=find(strcmp(u_CLabel_,'11')); end;
if nl==3; nA=find(strcmp(u_CLabel_,'19')); nB=find(strcmp(u_CLabel_,'15')); end;
disp(sprintf(' %% Comparing Label %d (%s) and %d (%s):',nA,u_CLabel_{nA},nB,u_CLabel_{nB}));
test_loader_plot_sort_list_0(AB_E_,nA,nB,u_CLabel_,C_VariableName_LC_,dir_trunk,str_prefix);
end;%for nl=1:3;
clear AB_E_;
%%%%%%%%;

str_prefix = sprintf('I_%s_r%d',xeta_.infix,nrank);
%%%%%%%%;
[AB_I_] = test_loader_helper_1(n_iteration,n_CLabel,n_CCOV,n_I_GENE,str_CLabel_,xeta_.C_rank_,I_rank_,xeta_.zeta_I_un_,xeta_.zeta_I_vn_,n_zeta_I,xeta_.beta_I_un_,xeta_.beta_I_vn_,n_beta_I);
test_loader_plot_absZ_avg__0(AB_I_,C_VariableName_LC_,dir_trunk,str_prefix);
test_loader_plot_C_absZ_0(AB_I_,C_VariableName_LC_,dir_trunk,str_prefix);
test_loader_plot_Z_absZ_0(AB_I_,dir_trunk,str_prefix);
test_loader_helper_tsv_0(AB_I_,dir_trunk,str_prefix);
%%%%%%%%;
for nl=1:3;
if nl==1; nA=find(strcmp(u_CLabel_,'28')); nB=find(strcmp(u_CLabel_,'33')); end;
if nl==2; nA=find(strcmp(u_CLabel_,'22')); nB=find(strcmp(u_CLabel_,'11')); end;
if nl==3; nA=find(strcmp(u_CLabel_,'19')); nB=find(strcmp(u_CLabel_,'15')); end;
disp(sprintf(' %% Comparing Label %d (%s) and %d (%s):',nA,u_CLabel_{nA},nB,u_CLabel_{nB}));
test_loader_plot_sort_list_0(AB_I_,nA,nB,u_CLabel_,C_VariableName_LC_,dir_trunk,str_prefix);
end;%for nl=1:3;
clear AB_I_;
%%%%%%%%;

str_prefix = sprintf('EI_%s_r%d',xeta_.infix,nrank);
%%%%%%%%;
[AB_EI_] = test_loader_helper_1(n_iteration,n_CLabel,n_CCOV,n_E_GENE + n_I_GENE,str_CLabel_,xeta_.C_rank_,[E_rank_ , I_rank_],xeta_.zeta_EI_un_,xeta_.zeta_EI_vn_,n_zeta_EI,xeta_.beta_EI_un_,xeta_.beta_EI_vn_,n_beta_EI);
test_loader_plot_absZ_avg__0(AB_EI_,C_VariableName_LC_,dir_trunk,str_prefix);
test_loader_plot_C_absZ_0(AB_EI_,C_VariableName_LC_,dir_trunk,str_prefix);
test_loader_plot_Z_absZ_0(AB_EI_,dir_trunk,str_prefix);
test_loader_helper_tsv_0(AB_EI_,dir_trunk,str_prefix);
%%%%%%%%;
for nl=1:3;
if nl==1; nA=find(strcmp(u_CLabel_,'28')); nB=find(strcmp(u_CLabel_,'33')); end;
if nl==2; nA=find(strcmp(u_CLabel_,'22')); nB=find(strcmp(u_CLabel_,'11')); end;
if nl==3; nA=find(strcmp(u_CLabel_,'19')); nB=find(strcmp(u_CLabel_,'15')); end;
disp(sprintf(' %% Comparing Label %d (%s) and %d (%s):',nA,u_CLabel_{nA},nB,u_CLabel_{nB}));
test_loader_plot_sort_list_0(AB_EI_,nA,nB,u_CLabel_,C_VariableName_LC_,dir_trunk,str_prefix);
end;%for nl=1:3;
clear AB_EI_;
%%%%%%%%;

end;%for nrank=1:n_rank;
