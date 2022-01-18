function C_ = impute_fill0(B_,ij_missed_);

ij_missed_ = sort(ij_missed_);
n_missed = numel(ij_missed_);
ij_filled_ = setdiff(1:numel(B_),ij_missed_);
n_filled = numel(ij_filled_);
[n_r,n_c] = size(B_);
C_ = B_;

%%%%%%%%;
% fill in missing data with zeros. ;
%%%%%%%%;
C_(ij_missed_) = 0;
