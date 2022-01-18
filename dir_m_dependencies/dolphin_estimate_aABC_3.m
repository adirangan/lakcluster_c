function ...
[ ...
 parameter ...
,a_ ...
,A__ ...
,BB__ ...
,CC__ ...
,L ...
,niteration ...
,a__ ...
,A___ ...
,BB___ ...
,CC___ ...
,L_ ...
,index_min ...
] = ...
dolphin_estimate_aABC_3( ...
 parameter ...
,aid_ ...
,age_ ...
,dat__ ...
);

if nargin<3; parameter = []; end;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'n_step'); parameter.n_step = 1; end;
if ~isfield(parameter,'n_iteration'); parameter.n_iteration = 4; end;
if ~isfield(parameter,'flag_constrain_CC'); parameter.flag_constrain_CC = 0; end;
if ~isfield(parameter,'flag_crude'); parameter.flag_crude = 0; end;

n_step = parameter.n_step;
n_iteration = parameter.n_iteration;
flag_verbose = parameter.flag_verbose;
flag_constrain_CC = parameter.flag_constrain_CC;
flag_crude = parameter.flag_crude;

[n_smp,n_var] = size(dat__);

u_aid_ = unique(aid_);
n_aid = numel(u_aid_);
n_aid_ = zeros(n_aid,1);
index_aid__ = cell(n_aid,1);
for naid=0:n_aid-1;
tmp_index_aid_ = efind(aid_==u_aid_(1+naid));
index_aid__{1+naid} = tmp_index_aid_;
n_aid_(1+naid) = numel(index_aid__{1+naid});
end;%for naid=0:n_aid-1;

tmp_t = tic();
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
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% dt_all_: time %0.3fs',tmp_t)); end;
dt_lim_ = [0;max(dt_all_)];

a_ = zeros(n_var,1);
A__ = zeros(n_var,n_var); 

tmp_t = tic();
index_mss_ = efind(~isfinite(dat__));
BB__ = []; CC__ = [];
%%%%%%%%;
% Now impute data using A. ;
%%%%%%%%;
dat_imp__ = dolphin_impute_aA_0(aid_,age_,dat__,a_,A__,index_mss_);
%%%%%%%%;
% Now estimate initial B,C. ;
%%%%%%%%;
if (~flag_crude); [BB__,CC__] = dolphin_estimate_BC_from_aA_1(aid_,age_,dat_imp__,a_,A__,n_step,flag_constrain_CC); end;
if ( flag_crude); [BB__,CC__] = dolphin_estimate_BC_from_aA_crude_1(aid_,age_,dat_imp__,a_,A__,n_step,flag_constrain_CC); end;
%%%%%%%%;
% Now estimate a_ and A__ from BB__ and CC__. ;
%%%%%%%%;
[a_,A__,L] = dolphin_estimate_aA_from_BC_1(aid_,age_,dat_imp__,BB__,CC__,dt_lim_,n_step);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% initial estimate: time %0.3fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% initial: negative-log-likelihood %0.16f',L)); end;
%%%%%%%%;
% Now iterate a few times. ;
%%%%%%%%;
L_ = zeros(n_iteration+1,1);
a__ = cell(1+n_iteration,1);
A___ = cell(1+n_iteration,1);
BB___ = cell(1+n_iteration,1);
CC___ = cell(1+n_iteration,1);
niteration=0;
L_old = L;
L_(1+niteration) = L_old;
a__{1+niteration} = a_;
A___{1+niteration} = A__;
BB___{1+niteration} = BB__;
CC___{1+niteration} = CC__;
flag_continue=1;
while (flag_continue);
tmp_t = tic();
% Re-impute data using A. ;
dat_imp__ = dolphin_impute_aA_0(aid_,age_,dat__,a_,A__,index_mss_);
% Re-estimate B,C. ;
if (~flag_crude); [BB__,CC__] = dolphin_estimate_BC_from_aA_1(aid_,age_,dat_imp__,a_,A__,n_step,flag_constrain_CC); end;
if ( flag_crude); [BB__,CC__] = dolphin_estimate_BC_from_aA_crude_1(aid_,age_,dat_imp__,a_,A__,n_step,flag_constrain_CC); end;
% Re-estimate a_ and A__ from BB__ and CC__. ;
[a_,A__,L] = dolphin_estimate_aA_from_BC_1(aid_,age_,dat_imp__,BB__,CC__,dt_lim_,n_step);
L_new = L; L_(1+niteration+1)=L;
a__{1+niteration+1} = a_;
A___{1+niteration+1} = A__;
BB___{1+niteration+1} = BB__;
CC___{1+niteration+1} = CC__;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% iteration %d: time %0.3fs',tmp_t,niteration)); end;
if (flag_verbose); disp(sprintf(' %% iteration %d: negative-log-likelihood %0.16f',niteration,L)); end;
flag_continue=0;
niteration=niteration+1;
if (niteration<n_iteration & fnorm(L_old-L_new)/fnorm(L_old)>1e-3); flag_continue=1; end;
L_old = L_new;
end;%while (flag_continue);

%%%%%%%%;
% pick lowest loss as output. ;
%%%%%%%%;
[~,index_min] = min(L_(1:1+niteration)); index_min = index_min-1;
a_ = a__{1+index_min};
A__ = A___{1+index_min};
BB__ = BB___{1+index_min};
CC__ = CC___{1+index_min};
L = L_(1+index_min);
