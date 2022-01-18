function [label_enum_,n_label_] = label_num_to_enum_0(label_num_);
% converts numeric labels to ordered numeric labels, ;
% Note that the numeric labels are 1-based. ;
% Several other functions depend on this feature (do not change). ;
n_u = numel(label_num_);
u_label_ = unique(label_num_);
n_label = length(u_label_);
n_label_ = zeros(n_label,1);
for nlabel = 1:n_label;
n_label_(nlabel) = length(find(label_num_ == u_label_(nlabel)));
end;%for nlabel = 1:n_label;
label_enum_ = zeros(n_u,1);
for nlabel=1:n_label;
tmp_ij_ = find(label_num_ == u_label_(nlabel));
if (~isempty(tmp_ij_));
label_enum_(tmp_ij_) = nlabel;
end;%if (~isempty(tmp_ij_));
end;%for nlabel=1:n_label;

