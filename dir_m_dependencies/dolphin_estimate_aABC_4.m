function ...
[ ...
 parameter ...
,a_ ...
,A__ ...
,BB_inv__ ...
,CC_inv__ ...
,L ...
,niteration ...
,a__ ...
,A___ ...
,BB_inv___ ...
,CC_inv___ ...
,L_ ...
,index_min ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ ...
,age_ ...
,dat__ ...
);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); aid_=[]; end; na=na+1;
if (nargin<1+na); age_=[]; end; na=na+1;
if (nargin<1+na); dat__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'verbose'); parameter.verbose = 0; end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
if ~isfield(parameter,'n_iteration'); parameter.n_iteration = 16; end;
if ~isfield(parameter,'flag_filter'); parameter.flag_filter = 0; end;
verbose = parameter.verbose;
tolerance_master = parameter.tolerance_master;
n_iteration = parameter.n_iteration;
flag_filter = parameter.flag_filter;

[n_smp,n_var] = size(dat__);

%%%%%%%%;
% Initialize a and A. ;
%%%%%%%%;
a_ = zeros(n_var,1);
A__ = zeros(n_var,n_var); 
%%%%%%%%;
index_mss_ = efind(~isfinite(dat__));
BB_inv__ = []; CC_inv__ = [];
%%%%%%%%;
% Now impute data using A. ;
%%%%%%%%;
dat_imp__ = dolphin_impute_aA_0(aid_,age_,dat__,a_,A__,index_mss_);
%%%%%%%%;
% Now estimate initial B,C. ;
%%%%%%%%;
[ ...
 parameter ...
,BB_inv__ ...
,CC_inv__ ...
] = ....
dolphin_estimate_BB_inv_CC_inv_from_aA_2( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_imp__ ...
,a_ ...
,A__ ...
);
%%%%%%%%;
% Now estimate a_ and A__ using unfiltered data. ;
%%%%%%%%;
[ ...
 parameter ...
,a_ ...
,A__ ...
,L ...
] = ...
dolphin_estimate_aA_from_BC_3( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_imp__ ...
,1*BB_inv__ ...
,1*CC_inv__ ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now iterate a few times. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
P_nvv___ = {};
L_ = zeros(n_iteration+1,1);
a__ = cell(1+n_iteration,1);
A___ = cell(1+n_iteration,1);
BB_inv___ = cell(1+n_iteration,1);
CC_inv___ = cell(1+n_iteration,1);
L_old = +Inf;
niteration=0;
flag_continue=1;
while (flag_continue);
tmp_t = tic();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  flag_filter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% First filter. ;
%%%%%%%%;
[ ...
 parameter ...
,dat_imp_fil__ ...
,P_nvv___ ...
] = ...
dolphin_estimate_kalman_filter_0( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_imp__ ...
,a_ ...
,A__ ...
,1*BB_inv__ ...
,1*CC_inv__ ...
,P_nvv___ ...
);
%%%%%%%%;
% Estimate modified B,C using filtered data. ;
%%%%%%%%;
[ ...
 parameter ...
,BB_inv_fil__ ...
,CC_inv_fil__ ...
] = ....
dolphin_estimate_BB_inv_CC_inv_from_aA_2( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_imp_fil__ ...
,a_ ...
,A__ ...
);
%%%%%%%%;
% Now estimate a_ and A__ using filtered data and modifed B,C. ;
%%%%%%%%;
[ ...
 parameter ...
,a_ ...
,A__ ...
,L ...
] = ...
dolphin_estimate_aA_from_BC_3( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_imp_fil__ ...
,1*BB_inv_fil__ ...
,1*CC_inv_fil__ ...
);
%{
%%%%%%%%;
% Now estimate B using filtered data. ;
%%%%%%%%;
[ ...
 parameter ...
,BB_inv__ ...
,~ ...
] = ....
dolphin_estimate_BB_inv_CC_inv_from_aA_2( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_imp_fil__ ...
,a_ ...
,A__ ...
);
 %}
%%%%%%%%;
% Now estimate B,C using unfiltered data. ;
%%%%%%%%;
[ ...
parameter ...
,BB_inv__ ...
,CC_inv__ ...
] = ....
dolphin_estimate_BB_inv_CC_inv_from_aA_2( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_imp__ ...
,a_ ...
,A__ ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if  flag_filter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~flag_filter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now estimate a_ and A__ using unfiltered data. ;
%%%%%%%%;
if (niteration> 0);
[ ...
 parameter ...
,a_ ...
,A__ ...
,L ...
] = ...
dolphin_estimate_aA_from_BC_3( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_imp__ ...
,1*BB_inv__ ...
,1*CC_inv__ ...
);
end;%if (niteration> 0);
%%%%%%%%;
% Estimate B,C using unfiltered data. ;
%%%%%%%%;
[ ...
parameter ...
,BB_inv__ ...
,CC_inv__ ...
] = ....
dolphin_estimate_BB_inv_CC_inv_from_aA_2( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_imp__ ...
,a_ ...
,A__ ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~flag_filter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now re-impute data using a and A. ;
%%%%%%%%;
dat_imp__ = dolphin_impute_aA_0(aid_,age_,dat__,a_,A__,index_mss_);
%%%%%%%%;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% niteration %d/%d: time %0.3fs',niteration,n_iteration,tmp_t)); end;
if (verbose); disp(sprintf(' %% negative-log-likelihood %0.2f',L)); end;
%%%%%%%%;
L_(1+niteration) = L;
a__{1+niteration} = a_;
A___{1+niteration} = A__;
BB_inv___{1+niteration} = BB_inv__;
CC_inv___{1+niteration} = CC_inv__;
L_err = fnorm(L - L_old)/max(1e-12,fnorm(L));
flag_continue = ( (niteration<n_iteration) & (L_err> tolerance_master) );
L_old = L;
niteration=niteration+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%while (flag_continue);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_iteration = niteration;
L_ = L_(1:n_iteration);
a__ = a__(1:n_iteration);
A___ = A___(1:n_iteration);
BB_inv___ = BB_inv___(1:n_iteration);
CC_inv___ = CC_inv___(1:n_iteration);

%%%%%%%%;
% pick lowest loss as output. ;
%%%%%%%%;
[~,index_min] = min(L_); index_min = index_min-1;
a_ = a__{1+index_min};
A__ = A___{1+index_min};
BB_inv__ = BB_inv___{1+index_min};
CC_inv__ = CC_inv___{1+index_min};
L = L_(1+index_min);