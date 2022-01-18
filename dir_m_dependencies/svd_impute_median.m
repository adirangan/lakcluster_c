function C_ = svd_impute_median(B_,ij_missed_);

ij_missed_ = sort(ij_missed_);
n_missed = numel(ij_missed_);
ij_filled_ = setdiff(1:numel(B_),ij_missed_);
n_filled = numel(ij_filled_);
[n_r,n_c] = size(B_);
C_ = B_;

%%%%%%%%;
% fill in missing data with median of each column. ;
%%%%%%%%;
ij_missed_pre=1; ij_missed_pos=1;
ij_filled_pre=1; ij_filled_pos=1;
for nc=1:n_c;
ij_sub_min = (nc-1)*n_r + 1;
ij_sub_max = (nc-1)*n_r + n_r;
while (ij_missed_pre<n_missed & ij_missed_(ij_missed_pre)<ij_sub_min); ij_missed_pre = ij_missed_pre+1; end;
while (ij_missed_pos<n_missed & ij_missed_(ij_missed_pos)<ij_sub_max); ij_missed_pos = ij_missed_pos+1; end;
ij_missed_sub_ = ij_missed_(ij_missed_pre:ij_missed_pos);
while (ij_filled_pre<n_filled & ij_filled_(ij_filled_pre)<ij_sub_min); ij_filled_pre = ij_filled_pre+1; end;
while (ij_filled_pos<n_filled & ij_filled_(ij_filled_pos)<ij_sub_max); ij_filled_pos = ij_filled_pos+1; end;
ij_filled_sub_ = ij_filled_(ij_filled_pre:ij_filled_pos);
if (length(ij_missed_sub_)>0 & length(ij_filled_sub_)>0);
tmp_m = median(B_(ij_filled_sub_));
C_(ij_missed_sub_) = tmp_m;
end;%if (length(ij_missed_sub_)>0 & length(ij_filled_sub_)>0);
end;%for nc=1:n_c;
