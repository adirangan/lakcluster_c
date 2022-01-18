function output_lpdf_ = tw_chiani_lpdf(x_0in_,beta);
if (nargin<2); beta=1; end;

k_appx_ = [46.44604884387787, 79.6594870666346, 0, 146.0206131050228];   %  K, THETA, ALPHA ;
theta_appx_ = [0.18605402228279347, 0.10103655775856243, 0, 0.05954454047933292];
alpha_appx_ = [9.848007781128567, 9.819607173436484, 0, 11.00161520109004];

output_lpdf_ = -Inf*ones(size(x_0in_));
x_ = x_0in_ + alpha_appx_(beta);
index_ = efind(x_>0);
t = theta_appx_(beta);
k = k_appx_(beta);
output_lpdf_(1+index_) = (k-1).*log(x_(1+index_)) - x_(1+index_)/t - k*log(t) - gammaln(k);

