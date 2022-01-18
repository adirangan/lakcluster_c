function ...
[ ...
 parameter ...
,BB_inv_avg_bkp__ ...
,CC_inv_avg_bkp__ ...
] = ...
dolphin_estimate_BC_from_aA_crude_2( ...
 parameter ...
,aid_ ...
,age_ ...
,dat__ ...
,a_ ...
,A__ ...
);
% Estimates BB_inv__ and CC_inv__, given a_ and A__. ;
% Uses a crude estimator designed to ensure that BB_inv__ and CC_inv_ are both positive semidefinite. ;
% If flag_constrain_CC then we constrain CC_inv to be diagonal. ;
% Estimates BB_inv__ and CC_inv__, given a_ and A__. ;
% if flag_constrain_CC then we constrain CC_inv to be diagonal. ;

verbose=0;
flag_check=0;

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

tmp_t = tic();
n_total = 0;
dt_all_ = zeros(n_total,1);
for naid=0:n_aid-1;
tmp_index_aid_ = index_aid__{1+naid};
tmp_dat__ = dat__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
for nstep=0:n_step-1;
if (n_age>n_step+1);
tmp_dt_sort_ = tmp_age_sort_(1+nstep+1:end)-tmp_age_sort_(1:end-(nstep+1));
for nage=1+nstep:n_age-1;
n_total = n_total+1;
end;%for nage=1+nstep:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% n_total: time %0.2fs',tmp_t)); end;

tmp_t = tic();
dt_all_ = zeros(n_total,1);
ntotal=0;
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
dt_all_(1+ntotal) = tmp_dt;
ntotal = ntotal+1;
end;%for nage=1+nstep:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
assert(ntotal==n_total);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% dt_all_: time %0.2fs',tmp_t)); end;

tmp_t = tic();
n_prctile = 32; flag_continue=1;
dt_threshold_ = prctile(dt_all_,linspace(0,100,1+n_prctile));
while (flag_continue);
index_percentile_from_dt_ = sum( repmat(dt_all_,[1,1+n_prctile]) >= repmat(reshape(dt_threshold_,[1,1+n_prctile]),[n_total,1]) , 2 ) - 1;
index_percentile_from_dt_ = max(0 , min(n_prctile-1 , index_percentile_from_dt_ ));
n_count_bkp_ = full( sparse( 1+index_percentile_from_dt_ , 1 , 1 , n_prctile , 1 ) );
flag_continue=0;
if (n_prctile>1 & min(n_count_bkp_)<16*n_var);
n_prctile = n_prctile-1;
dt_threshold_ = prctile(dt_all_,linspace(0,100,1+n_prctile));
flag_continue=1;
end;%if (n_prctile>1 & min(n_count_bkp_)<16*n_var);
end;%while (flag_continue);
if (verbose); disp(sprintf(' %% n_prctile %d',n_prctile)); end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% n_count_bkp_: time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
percentile_indicator__ = zeros(n_total,n_prctile);
for nprctile=0:n_prctile-1;
percentile_indicator__(1+efind(index_percentile_from_dt_==nprctile),1+nprctile)=1;
end;%for nprctile=0:n_prctile-1;
n_count_bkp_ = transpose(ones(1,n_total)*percentile_indicator__);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% percentile_indicator__: time %0.2fs',tmp_t)); end;

if flag_check;
tmp_t = tic();
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
if (verbose); disp(sprintf(' %% n_dt %d; n_prctile %d',n_dt,n_prctile)); end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% n_count_: time %0.2fs',tmp_t)); end;
disp(sprintf(' %% n_count_ vs n_count_bkp_: %0.16f',fnorm(n_count_ - n_count_bkp_)/fnorm(n_count_)));
end;%if flag_check;

%tmp_t = tic();
%dY_all_vt__ = zeros(n_var,n_total);
%DD_all_vvt__ = zeros(n_var^2,n_total);
%ntotal=0;
%for naid=0:n_aid-1;
%tmp_index_aid_ = index_aid__{1+naid};
%tmp_dat__ = dat__(1+tmp_index_aid_,:);
%tmp_age_ = age_(1+tmp_index_aid_);
%[tmp_age_sort_,tmp_index_age_] = sort(tmp_age_,'ascend'); tmp_index_age_ = tmp_index_age_-1; n_age = numel(tmp_index_age_);
%for nstep=0:n_step-1;
%if (n_age>n_step+1);
%tmp_dt_sort_ = tmp_age_sort_(1+nstep+1:end)-tmp_age_sort_(1:end-(nstep+1));
%for nage=1+nstep:n_age-1;
%tmp_dt = tmp_age_(1+tmp_index_age_(1+nage+0)) - tmp_age_(1+tmp_index_age_(1+nage-(1+nstep)));
%assert(fnorm(tmp_dt_sort_(1+nage-(1+nstep))-tmp_dt)<1e-6);
%Y_pre_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage-(1+nstep)),:),[n_var,1]);
%Y_pos_ = reshape(tmp_dat__(1+tmp_index_age_(1+nage+0),:),[n_var,1]);
%dY_ = (Y_pos_ - Y_pre_) - tmp_dt*(a_ + A__*Y_pre_);
%%dY_all_vt__(:,1+ntotal) = dY_;
%DD__ = dY_*transpose(dY_);
%DD_all_vt__(:,1+ntotal) = DD__(:);
%ntotal = ntotal+1;
%end;%for nage=1+nstep:n_age-1;
%end;%if (n_age>n_step+1);
%end;%for nstep=0:n_step-1;
%end;%for naid=0:n_aid-1;
%assert(ntotal==n_total);
%tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% dY_all_vt__: time %0.2fs',tmp_t)); end;

%tmp_t = tic();
%dt_avg_bkp_ = transpose(transpose(dt_all_)*percentile_indicator__)./n_count_bkp_;
%DD_avg_bkp___ = reshape(DD_all_vvt__*(percentile_indicator__*sparse(1:n_prctile,1:n_prctile,1./n_count_bkp_,n_prctile,n_prctile)),[n_var,n_var,n_prctile]);
%tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% DD_avg_bkp__: time %0.2fs',tmp_t)); end;

tmp_t = tic();
n_count_bkp_ = zeros(n_prctile,1);
dt_avg_bkp_ = zeros(n_prctile,1);
DD_avg_bkp___ = zeros(n_var,n_var,n_prctile);
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
n_count_bkp_(1+nprctile) = n_count_bkp_(1+nprctile) + 1;
dt_avg_bkp_(1+nprctile) = dt_avg_bkp_(1+nprctile) + tmp_dt;
DD_avg_bkp___(:,:,1+nprctile) = DD_avg_bkp___(:,:,1+nprctile) + DDj__;
end;%for nage=1:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for naid=0:n_aid-1;
for nprctile=0:n_prctile-1;
dt_avg_bkp_(1+nprctile) = dt_avg_bkp_(1+nprctile)/max(1,n_count_bkp_(1+nprctile));
DD_avg_bkp___(:,:,1+nprctile) = DD_avg_bkp___(:,:,1+nprctile)/max(1,n_count_bkp_(1+nprctile));
end;%for nprctile=0:n_prctile-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% DD_avg_bkp___: time %0.2fs',tmp_t)); end;

tmp_t = tic();
CC_inv_avg_bkp__ = DD_avg_bkp___(:,:,1+0);
if (flag_constrain_CC); CC_inv_avg_bkp__ = diag(diag(CC_inv_avg_bkp__)); end;
BB_inv_avg_bkp___ = zeros(n_var,n_var,n_prctile);
for nprctile=1:n_prctile-1;
BB_inv_avg_bkp___(:,:,1+nprctile) = (DD_avg_bkp___(:,:,1+nprctile) - CC_inv_avg_bkp__)/dt_avg_bkp_(1+nprctile);
[V_BB_inv__,S_BB_inv__] = eig(BB_inv_avg_bkp___(:,:,1+nprctile));
BB_inv_avg_bkp___(:,:,1+nprctile) = real(V_BB_inv__)*max(0,real(S_BB_inv__))*transpose(real(V_BB_inv__));
end;%for nprctile=1:n_prctile-1;
BB_inv_avg_bkp__ = sum(BB_inv_avg_bkp___,3) / n_prctile;
BB_inv_avg_bkp__ = real(BB_inv_avg_bkp__);
CC_inv_avg_bkp__ = real(CC_inv_avg_bkp__);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% BB_inv_avg_bkp__: time %0.2fs',tmp_t)); end;

if (flag_check);

BB_inv_avg___ = cell(n_prctile,1); for nprctile=0:n_prctile-1; BB_inv_avg___{1+nprctile} = zeros(n_var,n_var); end;
CC_inv_avg__ = zeros(n_var,n_var);

tmp_t = tic();
n_count_ = zeros(n_prctile,1);
dt_avg_ = zeros(n_prctile,1);
DD_avg___ = cell(n_prctile,1); for nprctile=0:n_prctile-1; DD_avg___{1+nprctile} = zeros(n_var,n_var); end;
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
%%%%%%%%;
for nprctile=0:n_prctile-1;
dt_avg_(1+nprctile) = dt_avg_(1+nprctile)/max(1,n_count_(1+nprctile));
DD_avg___{1+nprctile} = DD_avg___{1+nprctile}/max(1,n_count_(1+nprctile));
end;%for nprctile=0:n_prctile-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% DD_avg___: time %0.2fs',tmp_t)); end;

tmp_t = tic();
CC_inv_avg__ = DD_avg___{1+0};
if (flag_constrain_CC); CC_inv_avg__ = diag(diag(CC_inv_avg__)); end;
for nprctile=1:n_prctile-1;
BB_inv_avg___{1+nprctile} = (DD_avg___{1+nprctile} - CC_inv_avg__)/dt_avg_(1+nprctile);
[V_BB_inv__,S_BB_inv__] = eig(BB_inv_avg___{1+nprctile});
BB_inv_avg___{1+nprctile} = real(V_BB_inv__)*max(real(S_BB_inv__),0)*transpose(real(V_BB_inv__));
end;%for nprctile=1:n_prctile-1;
BB_inv_avg__ = zeros(n_var,n_var);
for nprctile=1:n_prctile-1;
BB_inv_avg__ = BB_inv_avg__ + BB_inv_avg___{1+nprctile};
end;%for nprctile=1:n_prctile-1;
BB_inv_avg__ = BB_inv_avg__ / n_prctile;
BB_inv_avg__ = real(BB_inv_avg__);
CC_inv_avg__ = real(CC_inv_avg__);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% BB_inv_avg__: time %0.2fs',tmp_t)); end;

disp(sprintf(' %% dt_avg_ vs dt_avg_bkp_: %0.16f',fnorm(dt_avg_ - dt_avg_bkp_)/fnorm(dt_avg_)));
disp(sprintf(' %% BB_inv_avg__ vs BB_inv_avg_bkp__: %0.16f',fnorm(BB_inv_avg__ - BB_inv_avg_bkp__)/fnorm(BB_inv_avg__)));
disp(sprintf(' %% CC_inv_avg__ vs CC_inv_avg_bkp__: %0.16f',fnorm(CC_inv_avg__ - CC_inv_avg_bkp__)/fnorm(CC_inv_avg__)));

end;%if (flag_check);
