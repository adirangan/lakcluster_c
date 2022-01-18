function ...
[ ...
 parameter ...
,P_vvt___ ...
] = ...
dolphin_estimate_kalman_asymptotic_covariance_0( ...
 parameter ...
,A__ ...
,BB_inv__ ...
,CC_inv__ ...
,dt_ ...
);

if nargin<1;
%%%%%%%%;
tolerance_master = 1e-2;
n_var = 2; T_max = 64;
if n_var==1; A_tru__ = -15/T_max*eye(n_var); end;
if n_var==2; A_tru__ = [-1/T_max , +1 ; -1 -2/T_max]; end;
a_tru_ = A_tru__*ones(n_var,1)*-0.3;
BB_tru__ = [];
BB_inv_tru__ = 0.01*eye(n_var);
CC_tru__ = [];
CC_inv_tru__ = 1.25*eye(n_var);
dt_ = linspace(0,0.25,1024);
parameter = struct('type','parameter');
parameter.tolerance_master = 1e-6;
parameter.n_iteration = 128;
[ ...
 parameter ...
,P_vvt___ ...
] = ...
dolphin_estimate_kalman_asymptotic_covariance_0( ...
 parameter ...
,A_tru__ ...
,BB_inv_tru__ ...
,CC_inv_tru__ ...
,dt_ ...
);
%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;

verbose=0;

na=0;
if (nargin<(1+na)); parameter=[]; end; na=na+1;
if (nargin<(1+na)); A__=[]; end; na=na+1;
if (nargin<(1+na)); BB_inv__=[]; end; na=na+1;
if (nargin<(1+na)); CC_inv__=[]; end; na=na+1;
if (nargin<(1+na)); dt_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end;
if (~isfield(parameter,'n_iteration')); parameter.n_iteration = 16; end;
tolerance_master = parameter.tolerance_master;
n_iteration = parameter.n_iteration;

if isempty(dt_); dt_ = linspace(0,0.1,11); end;
n_dt = numel(dt_);

n_var = size(A__,1);
P_vvt___ = zeros(n_var,n_var,n_dt);
P_old__ = ones(n_var,n_var);
for ndt=n_dt-1:-1:0;
dt = dt_(1+ndt);
F__ = eye(n_var,n_var) + A__*dt;
R__ = CC_inv__;
Q__ = BB_inv__*dt;
niteration=0;
flag_continue=1;
while flag_continue;
P_new__ = F__*(P_old__-P_old__*pinv(P_old__+R__,tolerance_master)*P_old__)*transpose(F__) + Q__;
P_err = fnorm(P_new__-P_old__)/max(1e-12,fnorm(P_new__));
flag_continue = (P_err> tolerance_master) & niteration<n_iteration;
P_old__ = P_new__;
niteration = niteration + 1;
end;%while flag_continue;
if (verbose);
disp(sprintf(' %% ndt %d/%d dt %0.2f: niteration %d P_err %0.6f',ndt,n_dt,dt,niteration,P_err));
disp(sprintf(' %% P_new__: ')); disp(P_new__);
end;%if (verbose);
P_vvt___(:,:,1+ndt) = P_new__;
end;%for ndt=n_dt-1:-1:0;


