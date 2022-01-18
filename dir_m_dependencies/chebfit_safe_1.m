function ...
[ ...
 c_obj ...
,x_use_ ...
,y_use_ ...
,bic_ ...      
,c_ ...
] = ...
chebfit_safe_1( ...
 x_ ...
,y_ ...
,n_degree ...
);

na=0;
if (nargin<1+na); x_ = []; end; na=na+1;
if (nargin<1+na); y_ = []; end; na=na+1;
if (nargin<1+na); n_degree = []; end; na=na+1;

if isempty(n_degree); n_degree = 32; end;

index_use_ = intersect(efind(isfinite(x_)),efind(isfinite(y_)));
n_data_use = numel(index_use_);
n_degree_use = min(n_degree,n_data_use-1);

x_use_ = x_(1+index_use_);
y_use_ = y_(1+index_use_);
x_lim_ = [min(x_use_),max(x_use_)];
C__ = chebpoly([0:n_degree_use],x_lim_);
A__ = C__(x_use_(:)); %<-- n_data_use by n_degree_use. ;

y_use_std = std(y_use_,1);
c__ = zeros(1+n_degree_use,1+n_degree_use);
c_obj_ = cell(1+n_degree_use,1);
bic_ = zeros(1+n_degree_use,1);
for ndegree_use=0:n_degree_use;
tmp_A__ = A__(:,1:1+ndegree_use);
c_ = tmp_A__\y_use_;
c_obj = C__(:,1:1+ndegree_use)*c_;
residual_ = y_use_ - c_obj(x_use_);
loglikelihood = sum(-residual_.^2/max(1e-12,2*y_use_std^2));
bic = (1+ndegree_use)*log(n_data_use) - 2*loglikelihood;
c__(1:1+ndegree_use,1+ndegree_use) = c_;
c_obj_{1+ndegree_use} = c_obj;
bic_(1+ndegree_use) = bic;
end;%for ndegree_use=0:n_degree_use;

[~,tmp_index] = min(bic_); tmp_index=tmp_index-1;
if ( (tmp_index==n_degree_use) & (n_degree_use+1> n_data_use) );
disp(sprintf(' %% Warning, min bic at max degree %d in chebfig_safe_1',n_degree_use));
end;%if ( (tmp_index==n_degree_use) & (n_degree_use+1> n_data_use) );
c_ = c__(:,1+tmp_index);
c_obj = c_obj_{1+tmp_index};

