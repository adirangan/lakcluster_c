function DDinv__ = dolphin_estimate_DD_0(dt,BB__,CC__,eig_tolerance);
if nargin<4; eig_tolerance = 1e-6; end;
DD__ = dt*BB__+CC__;
[Q__,S__] = eig(DD__);
S_ = diag(S__); S_ = max(eig_tolerance*max(S_),S_);
Sinv__ = diag(1./S_);
DDinv__ = Q__*Sinv__*transpose(Q__);
