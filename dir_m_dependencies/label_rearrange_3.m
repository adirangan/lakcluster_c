function ...
[ ...
 parameter ...
 ,label_rearrange_niteration ...
 ,output_label_ri__ ...
] = ...
label_rearrange_3( ...
 parameter ...
,input_label_ ...
,data_ ...
);
% rearranges numeric labels; 
% 0. calculates the intra-cluster mean (avg) and std.
% 1. calculates the distance between each point and the other cluster averages. ;
% 2. reassign each point that is *closer* to another cluster average (measured in terms of z-score) than to its own. ;
% if parameter.label_rearrange_flag_style==0, then ;
% std = average-distance-from-the-mean. ;
% if parameter.label_rearrange_flag_style==1, then ;
% std = maximum-likelihood std for isotropic-gaussian ( std^2 = 1/d * mean(squared-distance) ). ;
%%%%%%%%;

if (isempty(parameter)); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'label_rearrange_flag_style')); parameter.label_rearrange_flag_style = 1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'label_rearrange_n_iteration')); parameter.label_rearrange_n_iteration = 1; end; %<-- parameter_bookmark. ;
label_rearrange_flag_style = parameter.label_rearrange_flag_style;
label_rearrange_n_iteration = parameter.label_rearrange_n_iteration;

verbose=0;
if (verbose); disp(sprintf(' %% [entering label_rearrange_3] label_rearrange_flag_style %d',label_rearrange_flag_style)); end;
%%%%%%%%;
[n_row,n_col] = size(data_);
assert(numel(input_label_)==n_row);
output_label_ri__ = zeros(n_row,label_rearrange_n_iteration);
%%%%%%%%;
output_label_old_ = input_label_ ;
label_rearrange_niteration=0;
flag_continue=1;
while flag_continue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
enum_label_ = label_num_to_enum_0(output_label_old_); n_u = length(unique(enum_label_));
if (verbose); disp(sprintf(' %% n_u %d',n_u)); end;
n_u_ = zeros(n_u,1);
u_ij_ = cell(n_u,1);
for nu=1:n_u;
u_ij_{nu} = find(enum_label_==nu);
n_u_(nu) = length(u_ij_{nu});
end;%for nu=1:n_u;
if (verbose); disp(sprintf(' %% n_u_')); disp(num2str(transpose(n_u_))); end;
%%%%%%%%;
avg_ = zeros(n_u,n_col);
std_ = zeros(n_u,1);
for nu1=1:n_u;
tmp_data_ = data_(u_ij_{nu1},:);
avg_(nu1,:) = mean(tmp_data_,1);
if label_rearrange_flag_style==0;
std_(nu1) = mean(sqrt(sum(bsxfun(@minus,tmp_data_,reshape(avg_(nu1,:),[1,n_col])).^2,2)),1);
end;%if label_rearrange_flag_style==0;
if label_rearrange_flag_style==1;
std_(nu1) = sqrt( mean(sum(bsxfun(@minus,tmp_data_,reshape(avg_(nu1,:),[1,n_col])).^2,2),1) / max(1,n_col) );
end;%if label_rearrange_flag_style==1;
clear tmp_data_ ;
end;%for nu1=1:n_u;
if (verbose); disp(sprintf(' %% std_')); disp(num2str(transpose(std_))); end;
%%%%%%%%;
z__ = cell(n_u,1);
for nu1=1:n_u;
z__{nu1} = zeros(n_u_(nu1),n_u);
for nu2=1:n_u_(nu1);
tmp_ij = u_ij_{nu1}(nu2);
tmp_p_ = data_(tmp_ij,:);
for nu3=1:n_u;
%z__{nu1}(nu2,nu3) = fnorm( tmp_p_ - avg_(nu3,:) ) / max(1e-12,std_(nu3)) ; %<-- use inter-cluster std. ; 
z__{nu1}(nu2,nu3) = fnorm( tmp_p_ - avg_(nu3,:) ) / max(1e-12,std_(nu1)) ; %<-- use intra-cluster std. ;
end;%for nu3=1:n_u;
end;%for nu2=1:n_u_(nu1);
end;%for nu1=1:n_u;
%%%%%%%%;
na=0;
f__ = cell(n_u,1);
for nu1=1:n_u;
f__{nu1} = zeros(n_u_(nu1),1);
for nu2=1:n_u_(nu1);
tmp_ij = u_ij_{nu1}(nu2);
[~,f__{nu1}(nu2)] = min(z__{nu1}(nu2,:));
if (verbose>1);
if (f__{nu1}(nu2) ~= nu1);
disp(sprintf(' %% na %d: exchanging point nu2 %d(%d) : cluster %d(%d) --> %d',na,nu2,tmp_ij,nu1,enum_label_(tmp_ij),f__{nu1}(nu2)));
na=na+1;
end;%if (f__{nu1}(nu2) ~= nu1);
end;%if (verbose);
end;%for nu2=1:n_u_(nu1);
end;%for nu1=1:n_u;
%%%%%%%%;
output_label_ = output_label_old_;
for nu1=1:n_u;
for nu2=1:n_u_(nu1);
tmp_ij = u_ij_{nu1}(nu2);
output_label_(tmp_ij) = f__{nu1}(nu2);
end;%for nu2=1:n_u_(nu1);
end;%for nu1=1:n_u;
%%%%%%%%;
flag_continue = (label_rearrange_niteration < label_rearrange_n_iteration-1) & (fnorm(output_label_ - output_label_old_)>0);
output_label_ri__(:,1+label_rearrange_niteration) = output_label_;
output_label_old_ = output_label_;
label_rearrange_niteration = label_rearrange_niteration+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%while flag_continue;

tmp_n = label_rearrange_n_iteration - label_rearrange_niteration;
output_label_ri__(:,1+label_rearrange_niteration:label_rearrange_n_iteration) = repmat(output_label_,[1,tmp_n]); %<-- fix final entries. ;

if (verbose); disp(sprintf(' %% [finished label_rearrange_3] label_rearrange_flag_style %d',label_rearrange_flag_style)); end;
