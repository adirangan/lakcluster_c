function [DDinv__,detDDinv] = dolphin_estimate_DD_1(dt,BB__,CC__,eig_tolerance);
if nargin<4; eig_tolerance = 1e-6; end;
DD__ = dt*BB__+CC__;
[Q__,S__] = eig(DD__);
S_ = diag(S__); S_ = max(eig_tolerance,max(eig_tolerance*max(S_),S_));
Sinv_ = 1./S_;
Sinv__ = diag(Sinv_);
DDinv__ = Q__*Sinv__*transpose(Q__);
detDDinv = prod(Sinv_);
