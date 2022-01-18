function [a_,A__,L] = dolphin_estimate_aA_from_BC_0(aid_,age_,dat__,BB__,CC__,dt_lim_,n_step);
% Estimates a_ and A__, given BB__ and CC__;

if (nargin<7); n_step = 1; end;
  
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

tmp_t = tic();
chebfun_DD__ = chebfun(@(dt) chebout_dolphin_estimate_DD_0(dt,BB__,CC__) , reshape(dt_lim_,[1,2]) , 'splitting' , 'on' );
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% chebfun: time %0.2fs',tmp_t)); end;

E00_nn_ = zeros(n_var^2,1);
E01_nnm__ = zeros(n_var^2,n_var^1);
E10_mnn__ = zeros(n_var^1,n_var^2);
E11_nnmm__ = zeros(n_var^2,n_var^2);
F0_n_ = zeros(n_var,1);
F1_nm__ = zeros(n_var,n_var);

tmp_t = tic();
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nstep=0:n_step-1;
if (n_age>n_step+1);
tmp_dt_sort_ = tmp_age_sort_(1+nstep+1:end)-tmp_age_sort_(1:end-(nstep+1));
for nage=1+nstep:n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-(1+nstep)));
assert(fnorm(tmp_dt_sort_(1+nage-(1+nstep))-tmp_dt)<1e-6);
Y_pre_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-(1+nstep)),:),[n_var,1]);
Y_pos_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage+0),:),[n_var,1]);
dY_ = (Y_pos_ - Y_pre_);
YY__ = Y_pre_*transpose(Y_pre_);
DD__ = reshape(chebfun_DD__(tmp_dt),[n_var,n_var]);
E00_nn_ = E00_nn_ + tmp_dt^2 * reshape(DD__,[n_var^2,1]);
E01_nnm__ = E01_nnm__ + tmp_dt^2 * reshape(DD__,[n_var^2,1])*reshape(Y_pre_,[1,n_var]);
E10_mnn__ = E10_mnn__ + tmp_dt^2 * reshape(Y_pre_,[n_var,1])*reshape(DD__,[1,n_var^2]);
E11_nnmm__ = E11_nnmm__  + tmp_dt^2 * reshape(DD__,[n_var^2,1])*reshape(YY__,[1,n_var^2]) ;
F0_n_ = F0_n_ + tmp_dt * (DD__*dY_) ;
F1_nm__ = F1_nm__ + tmp_dt * (DD__*dY_)*transpose(Y_pre_) ;
end;%for nage=1+nstep:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
E00_nn__ = reshape(E00_nn_,[n_var^1,n_var^1]);
E01_nnm__ = reshape(E01_nnm__,[n_var^1,n_var^2]);
E10_nmn__ = reshape(permute(reshape(E10_mnn__,[n_var,n_var,n_var]),[2,1,3]),[n_var^2,n_var^1]);
E11_nmnm__ = reshape(permute(reshape(E11_nnmm__,[n_var,n_var,n_var,n_var]),[1,3,2,4]),[n_var^2,n_var^2]);
F1_nm_ = reshape(F1_nm__,[n_var^2,1]);
E_nmnm__ = [ E00_nn__ , E01_nnm__ ; E10_nmn__ , E11_nmnm__ ];
F_nm_ = [ F0_n_ ; F1_nm_ ];
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% E__ and F__: time %0.2fs',tmp_t)); end;

tmp_t = tic();
aA_ = pinv(E_nmnm__,1e-6) * F_nm_;
a_ = aA_(1:n_var);
A__ = reshape(aA_(n_var + [1:n_var^2]),[n_var,n_var]);
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% a_ and A__: time %0.2fs',tmp_t)); end;

a_ = real(a_);
A__ = real(A__);

L = 0;
if (nargout>2);
tmp_t = tic();
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nstep=0:n_step-1;
if (n_age>n_step+1);
tmp_dt_sort_ = tmp_age_sort_(1+nstep+1:end)-tmp_age_sort_(1:end-(nstep+1));
for nage=1+nstep:n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-(1+nstep)));
assert(fnorm(tmp_dt_sort_(1+nage-(1+nstep))-tmp_dt)<1e-6);
Y_pre_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-(1+nstep)),:),[n_var,1]);
Y_pos_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage+0),:),[n_var,1]);
dY_ = (Y_pos_ - Y_pre_);
YY__ = Y_pre_*transpose(Y_pre_);
DD__ = reshape(chebfun_DD__(tmp_dt),[n_var,n_var]);
tmp_L = transpose(dY_ - tmp_dt*(a_ + A__*Y_pre_)) * DD__ * (dY_ - tmp_dt*(a_ + A__*Y_pre_));
L = L + 0.5*n_var*log(2*pi) - 0.5*log(det(DD__)) + 0.5*tmp_L;
end;%for nage=1+nstep:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% L: time %0.2fs',tmp_t)); end;
end;%if (nargout>2);

