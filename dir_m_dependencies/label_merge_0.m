function output_label_ = label_merge_0(input_label_,data_);
% merges numeric labels; 
% 0. calculates the intra-cluster mean (avg) and average-distance-from-the-mean (std). ;
% 1. calculates the inter-cluster distance between means. ;
% 2. merges clusters with an inter-cluster distance that is smaller than *both* cluster stds. ; 
%%%%%%%%;
verbose=1;
if (verbose); disp(sprintf(' %% [entering label_merge_0]')); end;
%%%%%%%%;
enum_label_ = label_num_to_enum_0(input_label_); n_u = length(unique(enum_label_));
if (verbose); disp(sprintf(' %% n_u %d',n_u)); end;
n_u_ = zeros(n_u,1);
u_ij_ = cell(n_u,1);
for nu=1:n_u;
u_ij_{nu} = find(enum_label_==nu);
n_u_(nu) = length(u_ij_{nu});
end;%for nu=1:n_u;
if (verbose); disp(sprintf(' %% n_u_')); disp(num2str(transpose(n_u_))); end;
%%%%%%%%;
avg_ = zeros(n_u,size(data_,2));
std_ = zeros(n_u,1);
for nu1=1:n_u;
tmp_data_ = data_(u_ij_{nu1},:);
avg_(nu1,:) = mean(tmp_data_,1);
tmp_d_ = zeros(n_u_(nu1),1);
for nu2=1:n_u_(nu1);
tmp_d_ = fnorm(tmp_data_(nu2,:) - avg_(nu1,:));
end;%for nu2=1:n_u_(nu1);
std_(nu1) = mean(tmp_d_);
if (std_(nu1)<=1e-12); std_(nu1) = 1.0d0 ; end;
clear tmp_data_ tmp_d_ ;
end;%for nu1=1:n_u;
if (verbose); disp(sprintf(' %% std_')); disp(num2str(transpose(std_))); end;
%%%%%%%%;
d__ = zeros(n_u,n_u);
f__ = zeros(n_u,n_u);
for nu1=1:n_u;
for nu2=1:n_u;
d__(nu1,nu2) = fnorm(avg_(nu1,:)-avg_(nu2,:));
if (d__(nu1,nu2) < min(std_(nu1),std_(nu2))); f__(nu1,nu2) = 1; end;
end;%for nu2=1:n_u;
end;%for nu1=1:n_u;
if (verbose); disp(sprintf(' %% d__')); disp(num2str(d__)); end;
if (verbose); disp(sprintf(' %% f__')); disp(num2str(f__)); end;
%%%%%%%%;
disp(sprintf(' fix later: cluster based on f__^length(f__)>0 ')); 
%%%%%%%%;
if (verbose); disp(sprintf(' %% [finished label_merge_0]')); end;
