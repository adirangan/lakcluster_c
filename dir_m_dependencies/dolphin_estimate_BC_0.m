function [BB__,CC__,l2_R__,sum_1,sum_dt,sum_dtdt,sum_DDj__,sum_DDjdt__] = dolphin_estimate_BC_0(aid_,age_,dat_imp__,A__);
% Estimates BB__ and CC__, given A__. ;

verbose_flag=0;

u_aid_ = unique(aid_);
n_aid = numel(u_aid_);
n_aid_ = zeros(n_aid,1);
index_aid__ = cell(n_aid,1);
for naid=0:n_aid-1;
index_aid__{1+naid} = efind(aid_==u_aid_(1+naid));
n_aid_(1+naid) = numel(index_aid__{1+naid});
end;%for naid=0:n_aid-1;

[n_smp,n_var_aug] = size(dat_imp__);

sum_1 = 0;
sum_dt = 0;
sum_dtdt = 0;
sum_DDj__ = zeros(n_var_aug,n_var_aug);
sum_DDjdt__ = zeros(n_var_aug,n_var_aug);
l2_DD__ = zeros(n_var_aug,n_var_aug);
l2_BC__ = zeros(n_var_aug,n_var_aug);

for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat_imp__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
%%%%%%%%;
% single step. ;
%%%%%%%%;
tmp_dt_sort_ = diff(tmp_age_sort_);
for nage=1:n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-1));
assert(fnorm(tmp_dt_sort_(1+nage-1)-tmp_dt)<1e-6);
Y_pre_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-1),:),[n_var_aug,1]);
Y_pos_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage+0),:),[n_var_aug,1]);
dY_ = (Y_pos_ - Y_pre_) - tmp_dt*A__*Y_pre_;
DDj__ = dY_*transpose(dY_);
sum_1 = sum_1 + 1;
sum_dt = sum_dt + tmp_dt;
sum_dtdt = sum_dtdt + tmp_dt.^2;
sum_DDj__ = sum_DDj__ + DDj__;
sum_DDjdt__ = sum_DDjdt__ + DDj__*tmp_dt;
l2_DD__ = l2_DD__ + DDj__.^2;
end;%for nage=1:n_age-1;
%%%%%%%%;
% double step. ;
%%%%%%%%;
tmp_dt_sort_ = tmp_age_sort_(1+2:end)-tmp_age_sort_(1:end-2);
for nage=2:n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-2));
assert(fnorm(tmp_dt_sort_(1+nage-2)-tmp_dt)<1e-6);
Y_pre_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-2),:),[n_var_aug,1]);
Y_pos_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage+0),:),[n_var_aug,1]);
dY_ = (Y_pos_ - Y_pre_) - tmp_dt*A__*Y_pre_;
DDj__ = dY_*transpose(dY_);
sum_1 = sum_1 + 1;
sum_dt = sum_dt + tmp_dt;
sum_dtdt = sum_dtdt + tmp_dt.^2;
sum_DDj__ = sum_DDj__ + DDj__;
sum_DDjdt__ = sum_DDjdt__ + DDj__*tmp_dt;
l2_DD__ = l2_DD__ + DDj__.^2;
end;%for nage=1:n_age-1;
end;%for naid=0:n_aid-1;

BB__ = zeros(n_var_aug,n_var_aug);
CC__ = zeros(n_var_aug,n_var_aug);
for nvar=0:n_var_aug*n_var_aug-1;
tmp_AA__ = [ sum_dtdt , sum_dt ; sum_dt , sum_1 ];
tmp_DDj_ = [ sum_DDjdt__(1+nvar) ; sum_DDj__(1+nvar) ];
tmp_BC_ = tmp_AA__ \ tmp_DDj_;
BB__(1+nvar) = tmp_BC_(1+0);
CC__(1+nvar) = tmp_BC_(1+1);
end;%for nvar=0:n_var_aug*n_var_aug-1;

l2_BC__ = BB__.^2*sum_dtdt + CC__.^2*sum_1 + BB__.*CC__.*2*sum_dt;

l2_R__ = l2_DD__ - l2_BC__;
