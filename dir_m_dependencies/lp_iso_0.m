function [lp] = lp_iso_0(A_,B_);
% calculates the log of the probability that samples in A_ would be drawn from the isotropic gaussian formed from B_. ;
[n_A,d_A] = size(A_);
[n_B,d_B] = size(B_);
assert(d_A==d_B); d = d_B;
[mu_B_,sg_B] = gaussian_iso_0(B_);
C_ = A_ - repmat(mu_B_,n_A,1);
lp = -d*n_A*log(sg_B) - 0.5*d*n_A*log(2*pi) - sum(abs(C_).^2,'all')/(2*sg_B.^2);
