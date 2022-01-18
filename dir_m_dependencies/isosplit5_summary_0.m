function [isosplit5_n,isosplit5_q] = isosplit5_summary_0(data_,label_);
% data_ = n_sample by n_dimension. ;
% label_ = output of isosplit5. ;
% isosplit5_n = number of distinct labels. ;
% isosplit5_q = inverse quality measure: (median intra-cluster distance) / (median inter-cluster distance). ;
n_s = size(data_,1);
aa_ = sum(data_.^2,2);
ab_ = data_*transpose(data_);
d_ = repmat(aa_,1,n_s) + repmat(transpose(aa_),n_s,1) - 2*ab_ ;
d_ = sqrt(d_);
label_u_ = unique(label_);
isosplit5_n = length(label_u_);
isosplit5_q = 1.0;
if isosplit5_n>1;
l_ = repmat(label_(:),1,n_s) - repmat(transpose(label_(:)),n_s,1);
mask_intra_ = find(l_==0);
mask_inter_ = find(l_~=0);
if median(d_(mask_inter_))<=0; disp(sprintf(' %% Warning! zero inter-cluster distances in isosplit5_summary_0')); end;
isosplit5_q = median(d_(mask_intra_)) / median(d_(mask_inter_));
end;%if isosplit5_n>1;
