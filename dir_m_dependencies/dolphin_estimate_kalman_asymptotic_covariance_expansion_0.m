function ...
[ ...
 parameter ...
,P1__ ...
,P2__ ...
,P3__ ...
] = ...
dolphin_estimate_kalman_asymptotic_covariance_0( ...
 parameter ...
,A__ ...
,BB_inv__ ...
,CC_inv__ ...
,CC__ ...
);

if nargin<1;
%%%%%%%%;
tolerance_master = 1e-2;
n_var = 2; T_max = 64;
if n_var==1; A_tru__ = -15/T_max*eye(n_var); end;
if n_var==2; A_tru__ = [-1/T_max , +2 ; -1 -2/T_max]; end;
a_tru_ = A_tru__*ones(n_var,1)*-0.3;
BB_tru__ = [];
BB_inv_tru__ = 1.50*eye(n_var);
if n_var==2; BB_inv_tru__ = 1.50*0.1*[1 0.3 ; 0.3 1]; end;
CC_tru__ = [];
CC_inv_tru__ = 1.25*eye(n_var);
if n_var==2; CC_inv_tru__ = 1.25*[1.0 0.9 ; 0.9 1.0]; end;
dt_ = linspace(0,0.25,1024); n_dt = numel(dt_);
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
parameter = struct('type','parameter');
parameter.tolerance_master = 1e-6;
[ ...
 parameter ...
,P1__ ...
,P2__ ...
,P3__ ...
] = ...
dolphin_estimate_kalman_asymptotic_covariance_expansion_0( ...
 parameter ...
,A_tru__ ...
,BB_inv_tru__ ...
,CC_inv_tru__ ...
,CC_tru__ ...
);
P1__,;
P2__,;
P3__,;
%%%%%%%%;
P1_vvt___ = zeros(n_var,n_var,n_dt);
P2_vvt___ = zeros(n_var,n_var,n_dt);
P3_vvt___ = zeros(n_var,n_var,n_dt);
for ndt=0:n_dt-1;
dt = dt_(1+ndt);
Pe__ = P1__*sqrt(dt);
P1_vvt___(:,:,1+ndt) = Pe__;
Pe__ = P1__*sqrt(dt) + P2__*dt;
P2_vvt___(:,:,1+ndt) = Pe__;
Pe__ = P1__*sqrt(dt) + P2__*dt + P3__*dt*sqrt(dt);
P3_vvt___(:,:,1+ndt) = Pe__;
end;%for ndt=0:n_dt-1;
figure(1);clf;figbig;
for nvar0=0:n_var-1;
for nvar1=0:n_var-1;
subplot(n_var,n_var,1+nvar1+nvar0*n_var);
hold on;
plot(sqrt(dt_),squeeze(P_vvt___(1+nvar0,1+nvar1,:)),'k.');
plot(sqrt(dt_),squeeze(P1_vvt___(1+nvar0,1+nvar1,:)),'b-');
plot(sqrt(dt_),squeeze(P2_vvt___(1+nvar0,1+nvar1,:)),'g-');
plot(sqrt(dt_),squeeze(P3_vvt___(1+nvar0,1+nvar1,:)),'r-');
hold off;
xlim([0,max(sqrt(dt_))]); xlabel('sqrt(dt)');
ylim(prctile(P_vvt___(1+nvar0,1+nvar1,:),[0,100],'all'));
title(sprintf('%d,%d',nvar0,nvar1),'Interpreter','none');
end;%for nvar1=0:n_var-1;
end;%for nvar0=0:n_var-1;
%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;

verbose=0;

na=0;
if (nargin<(1+na)); parameter=[]; end; na=na+1;
if (nargin<(1+na)); A__=[]; end; na=na+1;
if (nargin<(1+na)); BB_inv__=[]; end; na=na+1;
if (nargin<(1+na)); CC_inv__=[]; end; na=na+1;
if (nargin<(1+na)); CC__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;

if isempty(CC__); CC__ = pinv(CC_inv__,tolerance_master); end;

n_var = size(A__,1);

S0__ = CC__;

Z1__ = [ zeros(n_var,n_var) , -S0__ ; -BB_inv__ , zeros(n_var,n_var) ];
[U__,D__] = eigs(Z1__,n_var,'smallestreal');
U0__ = U__(1 + 0*n_var + [0:n_var-1],:);
U1__ = U__(1 + 1*n_var + [0:n_var-1],:);
P1__ = U1__*pinv(U0__,tolerance_master);

S1__ = -CC__*P1__*CC__;

Z2A__ = CC__*P1__;
Z2Q__ = (A__*P1__ + P1__*transpose(A__) + P1__*CC__*P1__*CC__*P1__);
Z2__ = [ -Z2A__ , zeros(n_var,n_var) ; -Z2Q__ , +transpose(Z2A__) ];
[U__,D__] = eigs(Z2__,n_var,'smallestreal');
U0__ = U__(1 + 0*n_var + [0:n_var-1],:);
U1__ = U__(1 + 1*n_var + [0:n_var-1],:);
P2__ = U1__*pinv(U0__,tolerance_master);

S2__ = CC__*P1__*CC__*P1__*CC__ - CC__*P2__*CC__;
PSP2__ = P1__*CC__*P1__;
Z3A__ = CC__*P1__;
Z3Q__ = P1__*S2__*P1__ + P1__*S1__*P2__ + P2__*S1__*P1__ + P2__*CC__*P2__ + ...
  A__*PSP2__ + PSP2__*transpose(A__) + ...
  -A__*P2__ - P2__*transpose(A__) ...
  ;
Z3__ = [ -Z3A__ , zeros(n_var,n_var) ; +Z3Q__ , +transpose(Z3A__) ];
[U__,D__] = eigs(Z3__,n_var,'smallestreal');
U0__ = U__(1 + 0*n_var + [0:n_var-1],:);
U1__ = U__(1 + 1*n_var + [0:n_var-1],:);
P3__ = U1__*pinv(U0__,tolerance_master);


