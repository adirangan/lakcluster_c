function rank_estimate_onecut = rank_estimate_onecut_0(A__,p_val);
% estimates the rank of a matrix A__. ;
if (nargin<1);
rng(0);A__ = randn(128);
n_rank = 5;
[U__,S__,V__] = svd(A__); S_ = diag(S__); S_(1+n_rank:end) = S_(1+n_rank:end)*0.9;
A__ = U__*diag(S_)*transpose(V__);
rank_estimate_onecut = rank_estimate_onecut_0(A__);
disp('returning'); return;
end;%if (nargin<1);

if (nargin<2); p_val = 0.01; end;

verbose=1;
[n_row,n_col] = size(A__);
n_dim = min(n_row,n_col);
B__ = (A__ - mean(A__,'all'))/std(A__,1,'all');
S_B_ = svds(B__,n_dim);
svd_tracy_widom = svd_tracy_widom_0(n_dim,max(n_row,n_col),p_val);
rank_estimate_onecut = max(find(S_B_>=svd_tracy_widom));
if (isempty(rank_estimate_onecut)); rank_estimate_onecut=0; end;
if(verbose); figure(2);clf;plot(1:n_dim,S_B_,'ko-',1:n_dim,ones(1,n_dim)*svd_tracy_widom,'r-'); xlim([1,n_dim]); xlabel('rank');ylabel('sigma'); end;

