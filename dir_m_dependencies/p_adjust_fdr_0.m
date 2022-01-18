function output_ = p_adjust_fdr_0(input_);
n = length(input_);
[~,o_] = sort(input_,'descend');
[~,r_] = sort(o_);
i_ = transpose(n:-1:1);
output_ = min(1,cummin((n./i_) .* input_(o_)));
output_ = output_(r_);
