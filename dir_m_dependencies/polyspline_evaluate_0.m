function ...
[ ...
 A__ ...
,x_value_use_ ...
,y_value_use_ ...
] = ...
polyspline_evaluate_0( ...
 n_patch ...
,x_edge_ ...
,n_degree ...
,x_value_ ...
,y_value_ ...
);
% Builds evaluation matrix for spline with n_patch components. ;
% The degree of each component is n_degree (watch out for conditioning). ;
% x_edge_ of size 1+n_patch should contain the (sorted) patch boundaries. ;
% Constraints encode derivative matching on interior boundaries, ;
% and zero concavity at exterior boundaries. ;
% Note that the polynomial coefficients are stored in matlab order ;
% (i.e., highest power comes first). ;

assert(numel(x_edge_)==1+n_patch);

n_value = numel(x_value_);
assert(numel(y_value_)==n_value);

index_use_ = intersect(efind(isfinite(x_value_)),efind(isfinite(y_value_)));
n_value_use = numel(index_use_);
x_value_use_ = x_value_(1+index_use_);
y_value_use_ = y_value_(1+index_use_);

A__ = zeros(n_value_use,n_patch*(1+n_degree));
for npatch=0:n_patch-1;
%%%%;
if (npatch==0 && npatch==n_patch-1); %<-- only a single patch. ;
tmp_index_ = 0:n_value_use-1; %<-- use everything;
end;%if (npatch==0 && npatch==n_patch-1); %<-- only a single patch. ;
%%%%;
if (npatch==0 && npatch< n_patch-1); %<-- leftmost patch. ;
tmp_index_ = efind(x_value_use_<=x_edge_(1+npatch+1));
end;%if (npatch==0 && npatch< n_patch-1); %<-- leftmost patch. ;
%%%%;
if (npatch> 0 && npatch< n_patch-1); %<-- interior patch. ;
tmp_index_ = efind( (x_value_use_> x_edge_(1+npatch+0)) & (x_value_use_<=x_edge_(1+npatch+1)) );
end;%if (npatch> 0 && npatch< n_patch-1); %<-- interior patch. ;
%%%%;
if (npatch> 0 && npatch==n_patch-1); %<-- rightmost patch. ;
tmp_index_ = efind(x_value_use_> x_edge_(1+npatch+0));
end;%if (npatch> 0 && npatch==n_patch-1); %<-- rightmost patch. ;
%%%%;
n_value_use_local = numel(tmp_index_);
if (n_value_use_local< 1+n_degree);
disp(sprintf(' %% only %d/%d points at npatch %d/%d',n_value_use_local,1+n_degree,npatch,n_patch));
end;%if (n_value_use_local< 1+n_degree);
assert(n_value_use_local>=1+n_degree);
x_eval__ = reshape(x_value_use_(1+tmp_index_),[n_value_use_local,1]).^[n_degree:-1:0];
A__(1+tmp_index_,1+npatch*(1+n_degree)+[0:n_degree]) = x_eval__;
end;%for npatch=0:n_patch-1;
