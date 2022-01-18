function [c_obj,x_use_,y_use_,c_] = chebfit_safe_0(x_,y_,n);
index_use_ = intersect(efind(isfinite(x_)),efind(isfinite(y_)));
n_use = min(n,numel(index_use_)-1);
x_use_ = x_(1+index_use_);
y_use_ = y_(1+index_use_);
x_lim_ = [min(x_use_),max(x_use_)];
C__ = chebpoly([0:n_use],x_lim_);
A__ = C__(x_use_(:));
c_ = A__\y_use_;
c_obj = C__*c_;

