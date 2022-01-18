function [BB_avg__,CC_avg__] = dolphin_estimate_BC_from_aA_crude_0(aid_,age_,dat__,a_,A__,n_step);
% Estimates BB__ and CC__, given a_ and A__. ;
% Uses a crude estimator designed to ensure that BB__ and CC_ are both positive semidefinite. ;

if (nargin<6); n_step = 1; end;
verbose_flag=0;

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
n_prctile = 32; flag_continue=1;
dt_threshold_ = prctile(dt_all_,linspace(0,100,1+n_prctile));
while (flag_continue);
n_count_ = zeros(n_prctile,1);
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
nprctile = max(0,min(n_prctile-1,sum(dt_threshold_<=tmp_dt)-1));
n_count_(1+nprctile) = n_count_(1+nprctile) + 1;
end;%for nage=1:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
flag_continue=0;
if (n_prctile>1 & min(n_count_)<16*n_var);
n_prctile = n_prctile-1;
dt_threshold_ = prctile(dt_all_,linspace(0,100,1+n_prctile));
flag_continue=1;
end;%if (n_prctile>1 & min(n_count_)<16*n_var);
end;%while (flag_continue);
if (verbose_flag); disp(sprintf(' %% n_dt %d; n_prctile %d',n_dt,n_prctile)); end;

n_count_ = zeros(n_prctile,1);
dt_avg_ = zeros(n_prctile,1);
DD_avg___ = cell(n_prctile,1); for nprctile=0:n_prctile-1; DD_avg___{1+nprctile} = zeros(n_var,n_var); end;
BB_avg___ = cell(n_prctile,1); for nprctile=0:n_prctile-1; BB_avg___{1+nprctile} = zeros(n_var,n_var); end;
CC_avg__ = zeros(n_var,n_var);

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
dY_ = (Y_pos_ - Y_pre_) - tmp_dt*(a_ + A__*Y_pre_);
DDj__ = dY_*transpose(dY_);
nprctile = max(0,min(n_prctile-1,sum(dt_threshold_<=tmp_dt)-1));
n_count_(1+nprctile) = n_count_(1+nprctile) + 1;
dt_avg_(1+nprctile) = dt_avg_(1+nprctile) + tmp_dt;
DD_avg___{1+nprctile} = DD_avg___{1+nprctile} + DDj__;
end;%for nage=1:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;

for nprctile=0:n_prctile-1;
dt_avg_(1+nprctile) = dt_avg_(1+nprctile)/max(1,n_count_(1+nprctile));
DD_avg___{1+nprctile} = DD_avg___{1+nprctile}/max(1,n_count_(1+nprctile));
end;%for nprctile=0:n_prctile-1;

CC_avg__ = DD_avg___{1+0};
for nprctile=1:n_prctile-1;
BB_avg___{1+nprctile} = (DD_avg___{1+nprctile} - CC_avg__)/dt_avg_(1+nprctile);
[V_BB__,S_BB__] = eig(BB_avg___{1+nprctile});
BB_avg___{1+nprctile} = real(V_BB__)*max(real(S_BB__),0)*transpose(real(V_BB__));
end;%for nprctile=1:n_prctile-1;
BB_avg__ = zeros(n_var,n_var);
for nprctile=1:n_prctile-1;
BB_avg__ = BB_avg__ + BB_avg___{1+nprctile};
end;%for nprctile=1:n_prctile-1;
BB_avg__ = BB_avg__ / n_prctile;

BB_avg__ = real(BB_avg__);
CC_avg__ = real(CC_avg__);
