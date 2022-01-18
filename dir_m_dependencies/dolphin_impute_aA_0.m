function dat_out__ = dolphin_impute_aA_0(aid_,age_,dat_0in__,a_,A__,index_mss_);
% impute missing data given a_, A__, B__ and C__. ;

na=0;
if (nargin<1+na); aid_=[]; end; na=na+1;
if (nargin<1+na); age_=[]; end; na=na+1;
if (nargin<1+na); dat_0in__=[]; end; na=na+1;
if (nargin<1+na); a_=[]; end; na=na+1;
if (nargin<1+na); A__=[]; end; na=na+1;
if (nargin<1+na); index_mss_=[]; end; na=na+1;

verbose_flag=0;

u_aid_ = unique(aid_);
n_aid = numel(u_aid_);
n_aid_ = zeros(n_aid,1);
index_aid__ = cell(n_aid,1);
for naid=0:n_aid-1;
index_aid__{1+naid} = efind(aid_==u_aid_(1+naid));
n_aid_(1+naid) = numel(index_aid__{1+naid});
end;%for naid=0:n_aid-1;

[n_smp,n_var] = size(dat_0in__);
if isempty(a_); a_ = zeros(n_var,1); end;
if isempty(A__); A__ = zeros(n_var,n_var); end;
if isempty(index_mss_); index_mss_ = efind(~isfinite(dat_0in__)); end;

dat_out__ = dat_0in__;

flag_mss__ = zeros(size(dat_0in__));
flag_mss__(1+index_mss_) = 1;
flag_ret__ = ones(size(dat_0in__));
flag_ret__(1+index_mss_) = 0;

% for each column, fill missing values with median. ;
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat_0in__ = dat_0in__(1+tmp_index_aid_,:);
tmp_dat_out__ = tmp_dat_0in__;
tmp_flag_mss__ = flag_mss__(1+tmp_index_aid_,:);
tmp_flag_ret__ = flag_ret__(1+tmp_index_aid_,:);
for nvar=0:n_var-1;
tmp_index_mss_ = efind(tmp_flag_mss__(:,1+nvar));
tmp_index_ret_ = efind(tmp_flag_ret__(:,1+nvar));
tmp_dat_out__(1+tmp_index_mss_,1+nvar) = 0;
if (numel(tmp_index_mss_)>0 & numel(tmp_index_ret_)>0);
tmp_dat_out__(1+tmp_index_mss_,1+nvar) = median(tmp_dat_out__(1+tmp_index_ret_,1+nvar));
end;%if (numel(tmp_index_mss_)>0 & numel(tmp_index_ret_)>0);
end;%for nvar=0:n_var-1;
dat_out__(1+tmp_index_aid_,:) = tmp_dat_out__;
end;%for naid=0:n_aid-1;

flag_continue=0;
if (~isempty(A__) & ~(isempty(a_))); flag_continue=1; end;
niteration=0;
while(flag_continue);
dat_out_old__ = dat_out__;
% for each column, fill missing values with most likely estimate. ;
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat_out__ = dat_out__(1+tmp_index_aid_,:);
tmp_flag_mss__ = flag_mss__(1+tmp_index_aid_,:);
tmp_flag_ret__ = flag_ret__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[~,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nage=1:n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-1));
Y_pre_ = tmp_dat_out__(1+tmp_index_age_(1+nage-1),:);
Y_pos_ = tmp_dat_out__(1+tmp_index_age_(1+nage+0),:);
Y_est_ = reshape(Y_pre_,[n_var,1]) + tmp_dt*(a_ + A__*reshape(Y_pre_,[n_var,1]));
tmp_index_mss_ = efind(tmp_flag_mss__(1+tmp_index_age_(1+nage+0),:));
tmp_dat_out__(1+tmp_index_age_(1+nage+0),1+tmp_index_mss_) = reshape(Y_est_(1+tmp_index_mss_),[1,numel(tmp_index_mss_)]);
end;%for nage=1:n_age-1;
dat_out__(1+tmp_index_aid_,:) = tmp_dat_out__;
end;%for naid=0:n_aid-1;
dat_out_new__ = dat_out__;
I_diff = fnorm(dat_out_new__-dat_out_old__);
if (verbose_flag); disp(sprintf(' %% niteration %d; I_diff %0.16f',niteration,I_diff)); end;
flag_continue=0;
if (I_diff>0 & niteration<32); flag_continue=1; end;
niteration=niteration+1;
end;%while(flag_continue);
