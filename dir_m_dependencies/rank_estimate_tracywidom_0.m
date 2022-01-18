function rank_estimate_tracywidom = rank_estimate_tracywidom_0(A__,p_val_0in);
% iteratively estimates the rank of a matrix A__. ;
% compare with rank_estimate_onecut_0. ;
if (nargin<1);
rng(0);
A__ = randn(128,256);
n_rank = 5;
[U__,S__,V__] = svds(A__,min(size(A__))); S_ = diag(S__); S_(1+n_rank:end) = S_(1+n_rank:end)*0.925;
A__ = U__*diag(S_)*transpose(V__);
rank_estimate_tracywidom = rank_estimate_tracywidom_0(A__);
disp('returning'); return;
end;%if (nargin<1);

if (nargin<2); p_val_0in = 0.01; end;

verbose=1;
[n_row,n_col] = size(A__); n_dim = min(n_row,n_col);
rank_estimate_tracywidom = 0;
B__ = (A__ - mean(A__,'all'))/std(A__,1,'all');
p_val_use = p_val_0in/(1+p_val_0in); %<-- now sum(p_val_use.^[1:Inf]) = p_val_0in. ;
continue_flag=1;
while (continue_flag);
[U_B__,S_B__,V_B__] = svds(B__,n_dim); s_B = S_B__(1,1);
svd_tracy_widom = svd_tracy_widom_0(n_dim,max(n_row,n_col),p_val_use);
if (s_B>=svd_tracy_widom);
rank_estimate_tracywidom = rank_estimate_tracywidom+1;
B__ = transpose(U_B__(:,2:end))*B__;
B__ = (B__ - mean(B__,'all'))/std(B__,1,'all');
[n_row,n_col] = size(B__); n_dim = min(n_row,n_col);
continue_flag=1;
end;%if (s_B>=svd_tracy_widom);
if (s_B< svd_tracy_widom);
continue_flag=0;
end;%if (s_B< svd_tracy_widom);
end;%while (continue_flag);

if (verbose); 
S_A_ = svds(A__,min(size(A__))); n_dim = min(size(A__));
figure(2); clf;
plot(1:n_dim,S_A_,'ko-',1:rank_estimate_tracywidom,S_A_(1:rank_estimate_tracywidom),'rx','MarkerSize',15);
xlim([1,n_dim]); xlabel('rank');ylabel('sigma'); 
title(sprintf('rank %d',rank_estimate_tracywidom));
end;%if (verbose); 

