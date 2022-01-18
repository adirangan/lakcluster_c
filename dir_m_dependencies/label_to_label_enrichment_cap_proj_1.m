function [cap_out_] = label_to_label_enrichment_cap_proj_1(cap_ori_rsum_,cap_ori_csum_,cap_0in_,constraint_tolerance);
% project cap_0in_ onto admissible set defined by cap_ori_rsum_ and cap_ori_csum_. ;
% cap_ori_rsum_ is the vector of row-sums (i.e., sum(cap_ori_,2)). ;
% cap_ori_csum_ is the vector of col-sums (i.e., sum(cap_ori_,1)). ;
% constraint_tolerance is a tolerance. ;
if isempty(constraint_tolerance); constraint_tolerance = 1e-9; end;

verbose=0;
if (verbose); disp(sprintf(' %% [entering label_to_label_enrichment_cap_proj_1]')); end;
if (verbose>1); disp(cap_0in_); end;

n_row = numel(cap_ori_rsum_);
n_col = numel(cap_ori_csum_);

cap_out_ = cap_0in_;
n_iteration = 1024*8; niteration=0;
continue_flag=1;
while (continue_flag);
cap_out_rsum_ = sum(cap_out_,2); cap_out_csum_ = sum(cap_out_,1);
constraint_error_rsum_ = cap_out_rsum_ - cap_ori_rsum_;
constraint_error_csum_ = cap_out_csum_ - cap_ori_csum_;
constraint_error_zero_ = min(0,cap_out_);
constraint_error = ...
max([ ...
 max(abs(constraint_error_rsum_),[],'all') ... 
;max(abs(constraint_error_csum_),[],'all') ... 
;max(abs(constraint_error_zero_),[],'all') ...
]);
constraint_error = full(constraint_error);
if (constraint_error> constraint_tolerance);
if (mod(niteration,3)==0);
cap_out_rproj_ = cap_out_ - repmat(constraint_error_rsum_,[1,n_col])/n_col;
cap_out_ = cap_out_rproj_;
end;%if (mod(niteration,3)==0);
if (mod(niteration,3)==1);
cap_out_cproj_ = cap_out_ - repmat(constraint_error_csum_,[n_row,1])/n_row;
cap_out_ = cap_out_cproj_;
end;%if (mod(niteration,3)==1);
if (mod(niteration,3)==2);
%cap_out_zproj_ = cap_out_ - min(constraint_error_zero_,[],'all');
cap_out_zproj_ = max(0,cap_out_);
cap_out_ = cap_out_zproj_;
end;%if (mod(niteration,3)==2);
continue_flag = 1;
end;%if (constraint_error> constraint_tolerance);
if (constraint_error<=constraint_tolerance);
continue_flag = 0;
end;%if (constraint_error<=constraint_tolerance);
niteration=niteration+1; if (niteration>=n_iteration); continue_flag=0; end;
if (verbose); disp(sprintf(' %% niteration %d/%d, constraint_error %0.16f',niteration,n_iteration,constraint_error)); end;
if (verbose>1); disp(cap_out_); end;
end;%while (continue_flag);

if (verbose); disp(sprintf(' %% [finished label_to_label_enrichment_cap_proj_1] constraint_error %0.16f',constraint_error)); end;
