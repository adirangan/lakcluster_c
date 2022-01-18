function [BB__,CC__,l2_R__,sum_1,sum_dt,sum_dtdt,sum_DDj__,sum_DDjdt__] = dolphin_estimate_BC_from_aA_0(aid_,age_,dat__,a_,A__,n_step);
% Estimates BB__ and CC__, given a_ and A__. ;

if (nargin<6); n_step = 1; end;
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

sum_1 = 0;
sum_dt = 0;
sum_dtdt = 0;
sum_DDj__ = zeros(n_var,n_var);
sum_DDjdt__ = zeros(n_var,n_var);
l2_DD__ = zeros(n_var,n_var);
l2_BC__ = zeros(n_var,n_var);

tmp_t = tic();
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
sum_1 = sum_1 + 1;
sum_dt = sum_dt + tmp_dt;
sum_dtdt = sum_dtdt + tmp_dt.^2;
sum_DDj__ = sum_DDj__ + DDj__;
sum_DDjdt__ = sum_DDjdt__ + DDj__*tmp_dt;
l2_DD__ = l2_DD__ + DDj__.^2;
end;%for nage=1:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% sum_DDjdt__: time %0.2fs',tmp_t)); end;

tmp_t = tic();
BB__ = zeros(n_var,n_var);
CC__ = zeros(n_var,n_var);
for nvar=0:n_var*n_var-1;
tmp_AA__ = [ sum_dtdt , sum_dt ; sum_dt , sum_1 ];
tmp_DDj_ = [ sum_DDjdt__(1+nvar) ; sum_DDj__(1+nvar) ];
tmp_BC_ = tmp_AA__ \ tmp_DDj_;
BB__(1+nvar) = tmp_BC_(1+0);
CC__(1+nvar) = tmp_BC_(1+1);
end;%for nvar=0:n_var*n_var-1;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% BB__ and C__: time %0.2fs',tmp_t)); end;

l2_BC__ = BB__.^2*sum_dtdt + CC__.^2*sum_1 + BB__.*CC__.*2*sum_dt;

l2_R__ = l2_DD__ - l2_BC__;

flag_rectify = 1;
if flag_rectify;
tmp_t = tic();
[V_BB__,S_BB__] = eig(BB__);
BB__ = real(V_BB__)*max(real(S_BB__),0)*transpose(real(V_BB__));
[V_CC__,S_CC__] = eig(CC__);
CC__ = real(V_CC__)*max(real(S_CC__),0)*transpose(real(V_CC__));
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% rectification: time %0.2fs',tmp_t)); end;
end;%if flag_rectify;

BB__ = real(BB__);
CC__ = real(CC__);
