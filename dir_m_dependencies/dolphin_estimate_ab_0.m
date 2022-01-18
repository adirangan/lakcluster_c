function [a_,b_] = dolphin_estimate_ab_0(aid_,age_,dat__);
% Estimates simple linear drift a_ and offset b_. ;
% Both are assumed constant across all dolphins. ;
% Ignores missing data. ;
  
verbose_flag=1;

u_aid_ = unique(aid_);
n_aid = numel(u_aid_);
n_aid_ = zeros(n_aid,1);
index_aid__ = cell(n_aid,1);
for naid=0:n_aid-1;
index_aid__{1+naid} = efind(aid_==u_aid_(1+naid));
n_aid_(1+naid) = numel(index_aid__{1+naid});
end;%for naid=0:n_aid-1;

[n_smp,n_var] = size(dat__);

sum_1_ = zeros(n_var,1);;
sum_t_ = zeros(n_var,1);;
sum_tt_ = zeros(n_var,1);;
sum_Y_ = zeros(n_var,1);
sum_Yt_ = zeros(n_var,1);

for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nage=0:n_age-1;
tmp_t = tmp_age_sort_(1+nage);
Y_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage),:),[n_var,1]);
tmp_index_ = efind(isfinite(Y_));
if (numel(tmp_index_)>0);
sum_1_(1+tmp_index_) = sum_1_(1+tmp_index_) + 1;
sum_t_(1+tmp_index_) = sum_t_(1+tmp_index_) + tmp_t;
sum_tt_(1+tmp_index_) = sum_tt_(1+tmp_index_) + tmp_t^2;
sum_Y_(1+tmp_index_) = sum_Y_(1+tmp_index_) + Y_(1+tmp_index_);
sum_Yt_(1+tmp_index_) = sum_Yt_(1+tmp_index_) + Y_(1+tmp_index_)*tmp_t;
end;%if (numel(tmp_index_)>0);
end;%for nage=0:n_age-1;
end;%for naid=0:n_aid-1;

a_ = zeros(n_var,1);
b_ = zeros(n_var,1);
for nvar=0:n_var-1;
tmp_AA__ = [sum_tt_(1+nvar),sum_t_(1+nvar);sum_t_(1+nvar),sum_1_(1+nvar)];
tmp_rhs_ = [sum_Yt_(1+nvar);sum_Y_(1+nvar)];
tmp_ab_ = tmp_AA__\tmp_rhs_;
a_(1+nvar) = tmp_ab_(1+0);
b_(1+nvar) = tmp_ab_(1+1);
end;%for nvar=0:n_var-1;
