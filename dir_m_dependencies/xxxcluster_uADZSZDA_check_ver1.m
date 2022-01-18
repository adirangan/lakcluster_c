function flag_error = xxxcluster_uADZSZDA_check_ver1(M_n_rows_,M_n_cols,A_n_rind_,A_n_cind,Z_n_rind_,T_n_cols,T_n_,T_n_cind);

flag_error = 0;

nbins = length(M_n_rows_);
if (length(intersect(A_n_cind,1:M_n_cols))<length(A_n_cind));
disp(sprintf(' %% Warning! A_n_cind indices out of bounds in xxxcluster_uADZSZDA_check_ver1.m'));
flag_error = flag_error + 1;
end;%if (length(intersect(A_n_cind,1:M_n_cols))<length(A_n_cind));
if (length(intersect(T_n_cind,1:T_n_cols))<length(T_n_cind));
disp(sprintf(' %% Warning! T_n_cind indices out of bounds in xxxcluster_uADZSZDA_check_ver1.m'));
flag_error = flag_error + 1;
end;%if (length(intersect(T_n_cind,1:T_n_cols))<length(T_n_cind));

for nb1=0:nbins-1;
if (length(intersect(A_n_rind_{1+nb1},Z_n_rind_{1+nb1}))>0); 
disp(sprintf(' %% Warning! nonzero intersection(A_n_rind_,Z_n_rind_) for 1+nb1=%d in xxxcluster_uADZSZDA_check_ver1.m',1+nb1));
flag_error = flag_error + 1;
end; %if (length(intersect(A_n_rind_{1+nb1},Z_n_rind_{1+nb1}))>0); 
if (length(intersect(A_n_rind_{1+nb1},1:M_n_rows_(1+nb1)))<length(A_n_rind_{1+nb1}));
disp(sprintf(' %% Warning! A_n_rind_ indices out of bounds for for 1+nb1=%d in xxxcluster_uADZSZDA_check_ver1.m',1+nb1));
flag_error = flag_error + 1;
end;%if (length(intersect(A_n_rind_{1+nb1},1:M_n_rows_(1+nb1)))<length(A_n_rind_{1+nb1}));
if (length(intersect(Z_n_rind_{1+nb1},1:M_n_rows_(1+nb1)))<length(Z_n_rind_{1+nb1}));
disp(sprintf(' %% Warning! Z_n_rind_ indices out of bounds for for 1+nb1=%d in xxxcluster_uADZSZDA_check_ver1.m',1+nb1));
flag_error = flag_error + 1;
end;%if (length(intersect(Z_n_rind_{1+nb1},1:M_n_rows_(1+nb1)))<length(Z_n_rind_{1+nb1}));
end;%for nb1=0:nbins-1;
