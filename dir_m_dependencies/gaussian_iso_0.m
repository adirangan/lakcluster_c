function [mu_A_,sg_A] = gaussian_iso_0(A_);
% calculates the parameters for the isotropic gaussian associated with A_. ;
[n_A,d] = size(A_);
mu_A_ = mean(A_,1);
B_ = A_ - repmat(mu_A_,n_A,1);
sg_A = sqrt(sum(abs(B_).^2,'all')/(d*n_A));
