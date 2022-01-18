function [A__] = dolphin_estimate_A_0(aid_,age_,dat_imp__,BB__,CC__,dt_lim_);
% Estimates A__, given BB__ and CC__;

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

chebfun_DD__ = chebfun(@(dt) chebout_dolphin_estimate_DD_0(dt,BB__,CC__) , reshape(dt_lim_,[1,2]) , 'splitting' , 'on' );

E_nnmm__ = zeros(n_var_aug^2,n_var_aug^2);
F_nm__ = zeros(n_var_aug,n_var_aug);

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
dY_ = (Y_pos_ - Y_pre_);
YY__ = Y_pre_*transpose(Y_pre_);
DD__ = reshape(chebfun_DD__(tmp_dt),[n_var_aug,n_var_aug]);
E_nnmm__ = E_nnmm__  + tmp_dt^2 * reshape(DD__,[n_var_aug^2,1])*reshape(YY__,[1,n_var_aug^2]) ;
F_nm__ = F_nm__ + tmp_dt * ((DD__ + transpose(DD__))*dY_)*transpose(Y_pre_) ;
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
dY_ = (Y_pos_ - Y_pre_);
YY__ = Y_pre_*transpose(Y_pre_);
DD__ = reshape(chebfun_DD__(tmp_dt),[n_var_aug,n_var_aug]);
E_nnmm__ = E_nnmm__  + tmp_dt^2 * reshape(DD__,[n_var_aug^2,1])*reshape(YY__,[1,n_var_aug^2]) ;
F_nm__ = F_nm__ + tmp_dt * ((DD__ + transpose(DD__))*dY_)*transpose(Y_pre_) ;
end;%for nage=1:n_age-1;
end;%for naid=0:n_aid-1;

E_nmnm__ = reshape(permute(reshape(E_nnmm__,[n_var_aug,n_var_aug,n_var_aug,n_var_aug]),[1,3,2,4]),[n_var_aug^2,n_var_aug^2]);
F_nm_ = reshape(F_nm__,[n_var_aug^2,1]);
A_ = pinv(E_nmnm__,1e-6) * F_nm_ * 0.5;
A__ = reshape(A_,[n_var_aug,n_var_aug]);
