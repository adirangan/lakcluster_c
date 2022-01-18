function output_pdf_ = tw_chiani_pdf(x_0in_,beta,kta_);
if (nargin<2); beta=1; end;
if (nargin<3); kta_=[]; end;

k_appx_ = [46.44604884387787, 79.6594870666346, 0, 146.0206131050228];   %  K, THETA, ALPHA ;
theta_appx_ = [0.18605402228279347, 0.10103655775856243, 0, 0.05954454047933292];
alpha_appx_ = [9.848007781128567, 9.819607173436484, 0, 11.00161520109004];

if ( isempty(kta_)); alpha = alpha_appx_(beta); end;
if (~isempty(kta_)); alpha = kta_(1+2); end;
if ( isempty(kta_)); theta = theta_appx_(beta); end;
if (~isempty(kta_)); theta = kta_(1+1); end;
if ( isempty(kta_)); k = k_appx_(beta); end;
if (~isempty(kta_)); k = kta_(1+0); end;

output_pdf_ = zeros(size(x_0in_));
x_ = x_0in_ + alpha;
index_ = efind(x_>0);
t = theta;
k = k;
output_pdf_(1+index_) = 1/(gamma(k)*t^k) * x_(1+index_).^(k-1) .* exp(-x_(1+index_)/t);

