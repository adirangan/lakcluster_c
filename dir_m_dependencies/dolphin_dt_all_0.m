function dt_all_ = dolphin_dt_all_0(aid_,age_,n_step);

u_aid_ = unique(aid_);
n_aid = numel(u_aid_);
n_aid_ = zeros(n_aid,1);
index_aid__ = cell(n_aid,1);
for naid=0:n_aid-1;
tmp_index_aid_ = efind(aid_==u_aid_(1+naid));
index_aid__{1+naid} = tmp_index_aid_;
n_aid_(1+naid) = numel(index_aid__{1+naid});
end;%for naid=0:n_aid-1;

n_step = 1;
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
