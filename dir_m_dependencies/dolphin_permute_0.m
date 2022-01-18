function [dat_prm__] = dolphin_permute_0(aid_,age_,dat__,rseed);

u_aid_ = unique(aid_);
n_aid = numel(u_aid_);
n_aid_ = zeros(n_aid,1);
index_aid__ = cell(n_aid,1);
for naid=0:n_aid-1;
tmp_index_aid_ = efind(aid_==u_aid_(1+naid));
index_aid__{1+naid} = tmp_index_aid_;
n_aid_(1+naid) = numel(index_aid__{1+naid});
end;%for naid=0:n_aid-1;

rng(rseed);
dat_prm__ = zeros(size(dat__));
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_p_ = randperm(n_aid_(1+naid));
tmp_dat_prm__ = tmp_dat__(tmp_p_,:);
dat_prm__(1+tmp_index_aid_,:) = tmp_dat_prm__;
end;%for naid=0:n_aid-1;
