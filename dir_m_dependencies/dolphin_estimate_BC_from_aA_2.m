function ...
[ ...
 parameter ...
,BB_inv__ ...
,CC_inv__ ...
,l2_R__ ...
,sum_1 ...
,sum_dt ...
,sum_dtdt ...
,sum_DD_j__ ...
,sum_DD_j_dt__ ...
] = ...
dolphin_estimate_BC_from_aA_2( ...
 parameter ...
,aid_ ...
,age_ ...
,dat__ ...
,a_ ...
,A__ ...
);
% Estimates BB_inv__ and CC_inv__, given a_ and A__. ;
% if flag_constrain_CC then we constrain CC_inv to be diagonal. ;

verbose=0;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); aid_=[]; end; na=na+1;
if (nargin<1+na); age_=[]; end; na=na+1;
if (nargin<1+na); dat__=[]; end; na=na+1;
if (nargin<1+na); a_=[]; end; na=na+1;
if (nargin<1+na); A__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end;
if (~isfield(parameter,'n_step')); parameter.n_step = 1; end;
if (~isfield(parameter,'flag_constrain_CC')); parameter.flag_constrain_CC = 0; end;
tolerance_master = parameter.tolerance_master;
n_step = parameter.n_step;
flag_constrain_CC = parameter.flag_constrain_CC;

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
sum_DD_inv_j__ = zeros(n_var,n_var);
sum_DD_inv_j_dt__ = zeros(n_var,n_var);
l2_DD_inv__ = zeros(n_var,n_var);
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
DD_inv_j__ = dY_*transpose(dY_);
sum_1 = sum_1 + 1;
sum_dt = sum_dt + tmp_dt;
sum_dtdt = sum_dtdt + tmp_dt.^2;
sum_DD_inv_j__ = sum_DD_inv_j__ + DD_inv_j__;
sum_DD_inv_j_dt__ = sum_DD_inv_j_dt__ + DD_inv_j__*tmp_dt;
l2_DD_inv__ = l2_DD_inv__ + DD_inv_j__.^2;
end;%for nage=1:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% sum_DD_inv_j_dt__: time %0.2fs',tmp_t)); end;

tmp_t = tic();
BB_inv__ = zeros(n_var,n_var);
CC_inv__ = zeros(n_var,n_var);
tmp_AA__ = [ sum_dtdt , sum_dt ; sum_dt , sum_1 ];
for nvar0=0:n_var-1;
for nvar1=nvar0:n_var-1;
if (flag_constrain_CC & nvar1~=nvar0);
BB_inv__(1+nvar0,1+nvar1) = sum_DD_inv_j_dt__(1+nvar0,1+nvar1)/sum_dtdt;
CC_inv__(1+nvar0,1+nvar1) = 0;
else;
tmp_DD_inv_j_ = [ sum_DD_inv_j_dt__(1+nvar0,1+nvar1) ; sum_DD_inv_j__(1+nvar0,1+nvar1) ];
tmp_BC_ = tmp_AA__ \ tmp_DD_inv_j_;
BB_inv__(1+nvar0,1+nvar1) = tmp_BC_(1+0);
CC_inv__(1+nvar0,1+nvar1) = tmp_BC_(1+1);
end;%if;
BB_inv__(1+nvar1,1+nvar0) = BB_inv__(1+nvar0,1+nvar1);
CC_inv__(1+nvar1,1+nvar0) = CC_inv__(1+nvar0,1+nvar1);
end;%for nvar1=nvar0:n_var-1;
end;%for nvar0=0:n_var-1;

tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% BB_inv__ and CC_inv__: time %0.2fs',tmp_t)); end;

l2_BC__ = BB_inv__.^2*sum_dtdt + CC_inv__.^2*sum_1 + BB_inv__.*CC_inv__.*2*sum_dt;

l2_R__ = l2_DD_inv__ - l2_BC__;

flag_rectify = 1;
if flag_rectify;
tmp_t = tic();
[V_BB_inv__,S_BB_inv__] = eig(BB_inv__);
BB_inv__ = real(V_BB_inv__)*max(0,real(S_BB_inv__))*transpose(real(V_BB_inv__));
[V_CC_inv__,S_CC_inv__] = eig(CC_inv__);
CC_inv__ = real(V_CC_inv__)*max(0,real(S_CC_inv__))*transpose(real(V_CC_inv__));
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% rectification: time %0.2fs',tmp_t)); end;
end;%if flag_rectify;

BB_inv__ = real(BB_inv__);
CC_inv__ = real(CC_inv__);
