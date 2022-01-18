function [kld] = kld_iso_0(A_,B_);
% calculates a version of the (symmetrized) kullbak-leibler divergence between A and B. ;
% briefly: symmetric kld = lp(A|B) + lp(B|A) - lp(A|A) - lp(B|B). ;
% More specifically the asymmetric kld between A_ and B_ is defined as: ;
% kld(A,B) = \int \rho_{A} * log(\rho_{A}/\rho_{B}) d\Omega. ;
% Now in our case we approximate \rho_{A} and \rho_{B} by isotropic gaussians ;
% built using maximum-likelihood-estimation from the sampled data in A_ and B_. ;
% An important modification we make is to replace the area element \rho_{A} d\Omega ;
% with the sampled estimate (i.e., monte-carlo integration). ;
% More specifically: we calculate ;
% kld(A,B) = (1/n_A) * \Sum_{x sampled from A} log( \rho_{A}(x) / \rho_{B}(x) ). ;
% Similarly, ;
% kld(B,A) = (1/n_B) * \Sum_{x sampled from B} log( \rho_{B}(x) / \rho_{A}(x) ). ;
[n_A,d_A] = size(A_);
[n_B,d_B] = size(B_);
assert(d_A==d_B); d = d_B;
[mu_A_,sg_A] = gaussian_iso_0(A_);
[mu_B_,sg_B] = gaussian_iso_0(B_);
A_mu_A_ = A_ - repmat(mu_A_,n_A,1);
A_mu_B_ = A_ - repmat(mu_B_,n_A,1);
B_mu_A_ = B_ - repmat(mu_A_,n_B,1);
B_mu_B_ = B_ - repmat(mu_B_,n_B,1);
kld_A_A = -d*log(sg_A) - 0.5*d*log(2*pi) - sum(abs(A_mu_A_).^2,'all')/(2*sg_A.^2)/n_A;
kld_A_B = -d*log(sg_B) - 0.5*d*log(2*pi) - sum(abs(A_mu_B_).^2,'all')/(2*sg_B.^2)/n_A;
kld_B_A = -d*log(sg_A) - 0.5*d*log(2*pi) - sum(abs(B_mu_A_).^2,'all')/(2*sg_A.^2)/n_B;
kld_B_B = -d*log(sg_B) - 0.5*d*log(2*pi) - sum(abs(B_mu_B_).^2,'all')/(2*sg_B.^2)/n_B;
kld = + kld_A_B + kld_B_A - kld_A_A - kld_B_B;