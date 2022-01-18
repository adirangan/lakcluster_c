function [mr_A_prm_,mr_Z_prm_] = dcgxpander_uADZSZDA_shuffle_ver2(shuffle_num,M_n_rows_,M_n_cols,A_n_rind_,A_n_cind,Z_n_rind_,T_n_cols,T_n_,T_n_cind,J_c_rind_,K_c_cind);
% creates shuffled row masks respecting categorical-covariates and continuous-covariate sectors ;

n_study=length(M_n_rows_);
M2_n_rows_ = zeros(n_study,1);
A2_n_rind_ = cell(n_study,1);
for nstudy=1:n_study;
M2_n_rows_(nstudy) = M_n_rows_(nstudy) - length(J_c_rind_{nstudy});
A2_n_rind_{nstudy} = setdiff(A_n_rind_{nstudy},J_c_rind_{nstudy});
end;%for nstudy=1:n_study;
A2_n_cind = setdiff(A_n_cind,K_c_cind);

[mr_A_prm_,mr_Z_prm_] = xxxcluster_uADZSZDA_shuffle_ver2(shuffle_num,M_n_rows_,M_n_cols,A2_n_rind_,A2_n_cind,Z_n_rind_,T_n_cols,T_n_,T_n_cind);

for nstudy=1:n_study;
mr_A_prm_{nstudy}(J_c_rind_{nstudy}) = 1;
end;%for nstudy=1:n_study;
