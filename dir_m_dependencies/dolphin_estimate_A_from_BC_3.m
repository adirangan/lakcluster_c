function ...
[ ...
 parameter ...
,A__ ...
,L_bkp ...
] = ...
dolphin_estimate_A_from_BC_3( ...
 parameter ...
,aid_ ...
,age_ ...
,dat__ ...
,BB_inv__ ...
,CC_inv__ ...
);
% Assuming a_=0, estimates A__, given BB_inv__ and CC_inv__;
% Uses centered difference Y_{j+1}-Y_{j-1} rather than Y_{j+1}-Y_{j}. ;

verbose_flag=0;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); aid_=[]; end; na=na+1;
if (nargin<1+na); age_=[]; end; na=na+1;
if (nargin<1+na); dat__=[]; end; na=na+1;
if (nargin<1+na); BB_inv__=[]; end; na=na+1;
if (nargin<1+na); CC_inv__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end;
if (~isfield(parameter,'n_step')); parameter.n_step = 1; end;
tolerance_master = parameter.tolerance_master;
n_step = parameter.n_step;
  
flag_check=0;

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
n_total = 0;
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nstep=0:n_step-1;
for nage=2*(1+nstep):n_age-1;
n_total = n_total+1;
end;%for nage=2*(1+nstep):n_age-1;
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% n_total: time %0.2fs',tmp_t)); end;

tmp_t = tic();
dt_all_ = zeros(n_total,1);
ntotal=0;
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nstep=0:n_step-1;
tmp_dt_sort_ = tmp_age_sort_(1+2*(1+nstep):end)-tmp_age_sort_(1:end-2*(1+nstep));
for nage=2*(1+nstep):n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-2*(1+nstep)));
assert(fnorm(tmp_dt_sort_(1+nage-2*(1+nstep))-tmp_dt)<1e-6);
dt_all_(1+ntotal) = tmp_dt;
ntotal = ntotal+1;
end;%for nage=2*(1+nstep):n_age-1;
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
assert(ntotal==n_total);
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% dt_all_: time %0.2fs',tmp_t)); end;
dt2_all_ = dt_all_.^2;
dt_lim_ = [0,max(dt_all_)];

tmp_t = tic();
Y_mid_all_vt__ = zeros(n_var,n_total);
YY_mid_all_vvt__ = zeros(n_var^2,n_total);
dY_all_vt__ = zeros(n_var,n_total);
ntotal=0;
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nstep=0:n_step-1;
tmp_dt_sort_ = tmp_age_sort_(1+2*(1+nstep):end)-tmp_age_sort_(1:end-2*(1+nstep));
for nage=2*(1+nstep):n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-2*(1+nstep)));
assert(fnorm(tmp_dt_sort_(1+nage-2*(1+nstep))-tmp_dt)<1e-6);
Y_pre_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-2*(1+nstep)),:),[n_var,1]);
Y_mid_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-1*(1+nstep)),:),[n_var,1]);
Y_pos_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-0*(1+nstep)),:),[n_var,1]);
Y_mid_all_vt__(:,1+ntotal) = Y_mid_(:);
dY_ = (Y_pos_ - Y_pre_);
dY_all_vt__(:,1+ntotal) = dY_;
YY_mid__ = Y_mid_*transpose(Y_mid_);
YY_mid_all_vvt__(:,1+ntotal) = YY_mid__(:);
ntotal = ntotal+1;
end;%for nage=2*(1+nstep):n_age-1;
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
assert(ntotal==n_total);
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% Y_mid_all_vt__ and YY_mid_all_vvt__: time %0.2fs',tmp_t)); end;
Y_mid_all_tv__ = transpose(Y_mid_all_vt__);
YY_mid_all_tvv__ = transpose(YY_mid_all_vvt__);

tmp_t = tic();
[ ...
 parameter ...
,DD_inv_vvt__ ...
,det_DD_inv_t_ ...
] = ...
chebout_dolphin_estimate_DD_inv_3( ...
 parameter ...
,dt_all_ ...
,BB_inv__ ...
,CC_inv__ ...
);
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% DD_inv_vvt__: time %0.2fs',tmp_t)); end;

E00_bkp_nn_ = zeros(n_var^2,1);
tmp_t = tic();
E00_bkp_nn_ = DD_inv_vvt__*dt2_all_;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% E00_bkp_nn_: time %0.2fs',tmp_t)); end;

E01_bkp_nnm__ = zeros(n_var^2,n_var^1);
tmp_t = tic();
E01_bkp_nnm__ = DD_inv_vvt__*(sparse(1:n_total,1:n_total,dt2_all_,n_total,n_total)*Y_mid_all_tv__);
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% E01_bkp_nnm__: time %0.2fs',tmp_t)); end;

E10_bkp_mnn__ = transpose(E01_bkp_nnm__);

E11_bkp_nnmm__ = zeros(n_var^2,n_var^2);
tmp_t = tic();
E11_bkp_nnmm__ = DD_inv_vvt__*(sparse(1:n_total,1:n_total,dt2_all_,n_total,n_total)*YY_mid_all_tvv__);
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% E11_bkp_nnmm__: time %0.2fs',tmp_t)); end;

tmp_t = tic();
F0_bkp_n_ = zeros(n_var,1);
F1_bkp_nm__ = zeros(n_var,n_var);
for nt=0:n_total-1;
dY_ = dY_all_vt__(:,1+nt);
Y_mid_ = Y_mid_all_vt__(:,1+nt);
dt_DD_inv_dY_ = dt_all_(1+nt)*reshape(DD_inv_vvt__(:,1+nt),[n_var,n_var])*dY_;
dt_DD_inv_dYY__ = dt_DD_inv_dY_*transpose(Y_mid_);
F0_bkp_n_ = F0_bkp_n_ + dt_DD_inv_dY_;
F1_bkp_nm__ = F1_bkp_nm__ + dt_DD_inv_dYY__;
end;%for nt=0:n_total-1;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% F0_bkp_n_ and F1_bkp_nm__: time %0.2fs',tmp_t)); end;

if flag_check;
%%%%%%%%;
% check against brute force. ;
%%%%%%%%;
E00_nn_ = zeros(n_var^2,1);
E01_nnm__ = zeros(n_var^2,n_var^1);
E10_mnn__ = zeros(n_var^1,n_var^2);
E11_nnmm__ = zeros(n_var^2,n_var^2);
F0_n_ = zeros(n_var,1);
F1_nm__ = zeros(n_var,n_var);
%%%%%%%%;
tmp_t = tic();
chebfun_DD_inv__ = chebfun(@(dt) chebout_dolphin_estimate_DD_0(dt,BB_inv__,CC_inv__) , reshape(dt_lim_,[1,2]) , 'splitting' , 'on' );
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% chebfun: time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nstep=0:n_step-1;
tmp_dt_sort_ = tmp_age_sort_(1+2*(1+nstep):end)-tmp_age_sort_(1:end-2*(1+nstep));
for nage=2*(1+nstep):n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-2*(1+nstep)));
assert(fnorm(tmp_dt_sort_(1+nage-2*(1+nstep))-tmp_dt)<1e-6);
Y_pre_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-2*(1+nstep)),:),[n_var,1]);
Y_mid_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-1*(1+nstep)),:),[n_var,1]);
Y_pos_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-0*(1+nstep)),:),[n_var,1]);
dY_ = (Y_pos_ - Y_pre_);
YY_mid__ = Y_mid_*transpose(Y_mid_);
DD_inv__ = reshape(chebfun_DD_inv__(tmp_dt),[n_var,n_var]);
E00_nn_ = E00_nn_ + tmp_dt^2 * reshape(DD_inv__,[n_var^2,1]);
E01_nnm__ = E01_nnm__ + tmp_dt^2 * reshape(DD_inv__,[n_var^2,1])*reshape(Y_mid_,[1,n_var]);
E10_mnn__ = E10_mnn__ + tmp_dt^2 * reshape(Y_mid_,[n_var,1])*reshape(DD_inv__,[1,n_var^2]);
E11_nnmm__ = E11_nnmm__  + tmp_dt^2 * reshape(DD_inv__,[n_var^2,1])*reshape(YY_mid__,[1,n_var^2]) ;
F0_n_ = F0_n_ + tmp_dt * (DD_inv__*dY_) ;
F1_nm__ = F1_nm__ + tmp_dt * (DD_inv__*dY_)*transpose(Y_mid_) ;
end;%for nage=2*(1+nstep):n_age-1;
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% E__ and F__: time %0.2fs',tmp_t)); end;
%%%%%%%%;
disp(sprintf(' %% E00_nn_ vs E00_bkp_nn_: %0.16f',fnorm(E00_nn_ - E00_bkp_nn_)/fnorm(E00_nn_)));
disp(sprintf(' %% E01_nnm__ vs E01_bkp_nnm__: %0.16f',fnorm(E01_nnm__ - E01_bkp_nnm__)/fnorm(E01_nnm__)));
disp(sprintf(' %% E10_mnn__ vs E01_bkp_mnn__: %0.16f',fnorm(E10_mnn__ - E10_bkp_mnn__)/fnorm(E10_mnn__)));
disp(sprintf(' %% E11_nnmm__ vs E11_bkp_nnmm__: %0.16f',fnorm(E11_nnmm__ - E11_bkp_nnmm__)/fnorm(E11_nnmm__)));
disp(sprintf(' %% F0_n_ vs F0_bkp_n_: %0.16f',fnorm(F0_n_ - F0_bkp_n_)/fnorm(F0_n_)));
disp(sprintf(' %% F1_nm__ vs F1_bkp_nm__: %0.16f',fnorm(F1_nm__ - F1_bkp_nm__)/fnorm(F1_nm__)));
end;%if flag_check;

tmp_t = tic();
E00_bkp_nn__ = reshape(E00_bkp_nn_,[n_var^1,n_var^1]);
E01_bkp_nnm__ = reshape(E01_bkp_nnm__,[n_var^1,n_var^2]);
E10_bkp_nmn__ = reshape(permute(reshape(E10_bkp_mnn__,[n_var,n_var,n_var]),[2,1,3]),[n_var^2,n_var^1]);
E11_bkp_nmnm__ = reshape(permute(reshape(E11_bkp_nnmm__,[n_var,n_var,n_var,n_var]),[1,3,2,4]),[n_var^2,n_var^2]);
F1_bkp_nm_ = reshape(F1_bkp_nm__,[n_var^2,1]);
%E_nmnm__ = [ E00_bkp_nn__ , E01_bkp_nnm__ ; E10_bkp_nmn__ , E11_bkp_nmnm__ ];
%F_nm_ = [ F0_bkp_n_ ; F1_bkp_nm_ ];
E_nmnm__ = E11_bkp_nmnm__;
F_nm_ = F1_bkp_nm_;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% E__ and F__: time %0.2fs',tmp_t)); end;

tmp_t = tic();
%aA_ = pinv(E_nmnm__,tolerance_master) * F_nm_;
%a_ = aA_(1:n_var);
%A__ = reshape(aA_(n_var + [1:n_var^2]),[n_var,n_var]);
A__ = reshape(pinv(E_nmnm__,tolerance_master) * F_nm_,[n_var,n_var]);
%tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% a_ and A__: time %0.2fs',tmp_t)); end;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% A__: time %0.2fs',tmp_t)); end;

A__ = real(A__);

if (nargout>2);

tmp_t = tic();
L_bkp = 0;
for nt=0:n_total-1;
tmp_dt = dt_all_(1+nt);
dY_ = dY_all_vt__(:,1+nt);
Y_mid_ = Y_mid_all_vt__(:,1+nt);
YY__ = reshape(YY_mid_all_vvt__(:,1+nt),[n_var,n_var]);
DD_inv__ = reshape(DD_inv_vvt__(:,1+nt),[n_var,n_var]);
%tmp_L = transpose(dY_ - tmp_dt*(a_ + A__*Y_mid_)) * DD_inv__ * (dY_ - tmp_dt*(a_ + A__*Y_mid_));
tmp_L = transpose(dY_ - tmp_dt*(A__*Y_mid_)) * DD_inv__ * (dY_ - tmp_dt*(A__*Y_mid_));
det_DD_inv = det_DD_inv_t_(1+nt);
L_bkp = L_bkp + 0.5*n_var*log(2*pi) - 0.5*log(det_DD_inv) + 0.5*tmp_L;
end;%for nt=0:n_total-1;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% L_bkp: time %0.2fs',tmp_t)); end;

if flag_check;
%%%%%%%%;
% check against brute force. ;
%%%%%%%%;
L = 0;
tmp_t = tic();
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nstep=0:n_step-1;
tmp_dt_sort_ = tmp_age_sort_(1+2*(1+nstep):end)-tmp_age_sort_(1:end-2*(1+nstep));
for nage=2*(1+nstep):n_age-1;
tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-2*(1+nstep)));
assert(fnorm(tmp_dt_sort_(1+nage-2*(1+nstep))-tmp_dt)<1e-6);
Y_pre_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-2*(1+nstep)),:),[n_var,1]);
Y_mid_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-1*(1+nstep)),:),[n_var,1]);
Y_pos_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-0*(1+nstep)),:),[n_var,1]);
dY_ = (Y_pos_ - Y_pre_);
YY_mid__ = Y_mid_*transpose(Y_mid_);
DD_inv__ = reshape(chebfun_DD_inv__(tmp_dt),[n_var,n_var]);
%tmp_L = transpose(dY_ - tmp_dt*(a_ + A__*Y_mid_)) * DD_inv__ * (dY_ - tmp_dt*(a_ + A__*Y_mid_));
tmp_L = transpose(dY_ - tmp_dt*(A__*Y_mid_)) * DD_inv__ * (dY_ - tmp_dt*(A__*Y_mid_));
L = L + 0.5*n_var*log(2*pi) - 0.5*log(det(DD_inv__)) + 0.5*tmp_L;
end;%for nage=1+nstep:n_age-1;
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
tmp_t = toc(tmp_t); if (verbose_flag); disp(sprintf(' %% L: time %0.2fs',tmp_t)); end;
disp(sprintf(' %% L vs L_bkp: %0.16f',fnorm(L - L_bkp)/fnorm(L)));
end;%if flag_check;

end;%if (nargout>2);

