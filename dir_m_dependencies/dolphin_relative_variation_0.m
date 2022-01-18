function [R,R_] = dolphin_relative_variation_0(aid_,age_,dat__,n_step);
% calculates relative variation. ;

if (nargin<4); n_step = 1; end;
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

dt_all_ = [];
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_] = sort(tmp_age_,'ascend'); tmp_index_ = tmp_index_-1;
for nstep=0;n_step-1;
dt_ = tmp_age_sort_(1+(1+nstep):end)-tmp_age_sort_(1:end-(1+nstep));
dt_all_ = [dt_all_;dt_];
end;%for nstep=0;n_step-1;
end;%for naid=0:n_aid-1;

n_dt = numel(dt_all_);
R_ = zeros(n_dt,1);
ndt=0;
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nstep=0:n_step-1;
if (n_age>n_step+1);
tmp_dt_sort_ = tmp_age_sort_(1+(1+nstep):end)-tmp_age_sort_(1:end-(1+nstep));
for nage=(1+nstep):n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-(1+nstep)));
assert(fnorm(tmp_dt_sort_(1+nage-(1+nstep))-tmp_dt)<1e-6);
Y_pre_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-(1+nstep)),:),[n_var,1]);
Y_pos_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage+0),:),[n_var,1]);
dY_ = (Y_pos_ - Y_pre_);
tmp_index_ = efind(isfinite(dY_));
R_(1+ndt) = fnorm(dY_(1+tmp_index_))/fnorm(Y_pre_(1+tmp_index_));
ndt=ndt+1;
end;%for nage=1:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;

R = mean(R_);
