function [A_n_cols,Y_n_cols,T_n_cols,A_n_rows_,Z_n_rows_] = lakcluster_uADZSZDA_check_0(shuffle_num,M_n_,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind);

nbins = length(M_n_);
A_n_cols_ = zeros(nbins,1); for nb1=0:nbins-1; A_n_cols_(1+nb1) = size(M_n_{1+nb1},2); end; if (length(unique(A_n_cols_))>1); disp(sprintf(' %% Warning! A_n_cols_ not unique in lakcluster_uADZSZDA_check_0.m')); end; A_n_cols = unique(A_n_cols_);
Y_n_cols_ = []; Y_n_cols = 0;
T_n_cols_ = zeros(nbins,1); for nb1=0:nbins-1; T_n_cols_(1+nb1) = size(T_n_{1+nb1},2); end; if (length(unique(T_n_cols_))>1); disp(sprintf(' %% Warning! T_n_cols_ not unique in lakcluster_uADZSZDA_check_0.m')); end; T_n_cols = unique(T_n_cols_);
if (length(intersect(A_n_cind,1:A_n_cols))<length(A_n_cind));
disp(sprintf(' %% Warning! A_n_cind indices out of bounds in lakcluster_uADZSZDA_check_0.m'));
end;%if (length(intersect(A_n_cind,1:A_n_cols))<length(A_n_cind));
if (length(intersect(T_n_cind,1:T_n_cols))<length(T_n_cind));
disp(sprintf(' %% Warning! T_n_cind indices out of bounds in lakcluster_uADZSZDA_check_0.m'));
end;%if (length(intersect(T_n_cind,1:T_n_cols))<length(T_n_cind));

A_n_rows_ = zeros(nbins,1);
Z_n_rows_ = zeros(nbins,1);

for nb1=0:nbins-1;
if (length(intersect(A_n_rind_{1+nb1},Z_n_rind_{1+nb1}))>0); 
disp(sprintf(' %% Warning! nonzero intersection(A_n_rind_,Z_n_rind_) for 1+nb1=%d in lakcluster_uADZSZDA_check_0.m',1+nb1));
end; %if (length(intersect(A_n_rind_{1+nb1},Z_n_rind_{1+nb1}))>0); 
if (length(intersect(A_n_rind_{1+nb1},1:size(M_n_{1+nb1},1)))<length(A_n_rind_{1+nb1}));
disp(sprintf(' %% Warning! A_n_rind_ indices out of bounds for for 1+nb1=%d in lakcluster_uADZSZDA_check_0.m',1+nb1));
end;%if (length(intersect(A_n_rind_{1+nb1},1:size(M_n_{1+nb1},1)))<length(A_n_rind_{1+nb1}));
if (length(intersect(Z_n_rind_{1+nb1},1:size(M_n_{1+nb1},1)))<length(Z_n_rind_{1+nb1}));
disp(sprintf(' %% Warning! Z_n_rind_ indices out of bounds for for 1+nb1=%d in lakcluster_uADZSZDA_check_0.m',1+nb1));
end;%if (length(intersect(Z_n_rind_{1+nb1},1:size(M_n_{1+nb1},1)))<length(Z_n_rind_{1+nb1}));
A_n_rows_(1+nb1) = size(M_n_{1+nb1},1);
Z_n_rows_(1+nb1) = size(M_n_{1+nb1},1);
end;%for nb1=0:nbins-1;
