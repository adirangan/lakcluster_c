function ...
[ ...
 parameter ...
,age_ ...
,X__ ...
,Y__ ...
,R_avg ...
] = ...
SDE_generate_data_0( ...
 parameter ...
,a_tru_ ...
,A_tru__ ...
,BB_tru__ ...
,BB_inv_tru__ ...
,CC_tru__ ...
,CC_inv_tru__ ...
,X_0_ ...
);

if nargin<1;
tolerance_master = 1e-6;
n_var = 2; T_max = 64;
if n_var==1; A_tru__ = -1*eye(n_var); end;
if n_var==2; A_tru__ = [-1/T_max , +1 ; -1 -2/T_max]; end;
a_tru_ = A_tru__*ones(n_var,1)*0.3;
A_inv_a_ = pinv(A_tru__,tolerance_master)*a_tru_;
disp(sprintf(' %% A_inv_a_: '));
disp(A_inv_a_);
if n_var==1; BB_tru__ = 1.5*eye(n_var); end;
if n_var==2;
B_inv__ = [ +1 , 0.1 ; -1 , 0.2]; BB_inv_tru__ = B_inv__*transpose(B_inv__);
BB_tru__ = pinv(BB_inv_tru__,tolerance_master);
end;%if n_var==2;
if n_var==1; CC_tru__ = 1*eye(n_var); end;
if n_var==2;
C_inv__ = [ +0.5 , +0.3 ; -0.5 , +0.7]; CC_inv_tru__ = C_inv__*transpose(C_inv__);
CC_tru__ = pinv(CC_inv_tru__,tolerance_master);
end;%if n_var==2;
parameter = struct('type','parameter');
parameter.dt_avg = 0.15;
parameter.T_max = T_max;
parameter.rseed = 1;
[ ...
 parameter ...
,age_ ...
,X__ ...
,Y__ ...
,R_avg ...
] = ...
SDE_generate_data_0( ...
 parameter ...
,a_tru_ ...
,A_tru__ ...
,BB_tru__ ...
,BB_inv_tru__ ...
,CC_tru__ ...
,CC_inv_tru__ ...
);
%%%%;
figure;clf;figbig;figbeach();
subplot(2,2,1);plot(age_,X__(1,:),'r-',age_,X__(2,:),'b-');
xlabel('time');ylabel('X'); xlim([0,T_max]);
subplot(2,2,2);s = surfline_0(X__(1,:),X__(2,:),age_); set(s,'LineWidth',3);
xlabel('x0');ylabel('x1');
subplot(2,2,3);plot(age_,Y__(1,:),'r-',age_,Y__(2,:),'b-');
xlabel('time');ylabel('Y'); xlim([0,T_max]);
subplot(2,2,4);s = surfline_0(Y__(1,:),Y__(2,:),age_); set(s,'LineWidth',3);
xlabel('y0');ylabel('y1');
%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); a_tru_=[]; end; na=na+1;
if (nargin<1+na); A_tru__=[]; end; na=na+1;
if (nargin<1+na); BB_tru__=[]; end; na=na+1;
if (nargin<1+na); BB_inv_tru__=[]; end; na=na+1;
if (nargin<1+na); CC_tru__=[]; end; na=na+1;
if (nargin<1+na); CC_inv_tru__=[]; end; na=na+1;
if (nargin<1+na); X_0_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end;
if (~isfield(parameter,'dt_avg')); parameter.dt_avg = 0.25; end;
if (~isfield(parameter,'T_max')); parameter.T_max = 64; end;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end;
tolerance_master = parameter.tolerance_master;
dt_avg = parameter.dt_avg;
T_max = parameter.T_max;
rseed = parameter.rseed;
n_T = floor(T_max/parameter.dt_avg);

n_var = numel(a_tru_);
[Psi__,Lambda__] = eig(A_tru__); Psi_inv__ = inv(Psi__);
Lambda_ = diag(Lambda__);
A_inv_a_ = pinv(A_tru__,tolerance_master)*a_tru_;
%%%%%%%%;
if ( isempty(BB_inv_tru__));
BB_inv_tru__ = pinv(BB_tru__,tolerance_master);
end;%if ( isempty(BB_inv_tru__));
B_pinv_tru__ = sqrtm(BB_inv_tru__);
if ( isempty(CC_inv_tru__));
CC_inv_tru__ = pinv(CC_tru__,tolerance_master);
end;%if ( isempty(CC_inv_tru__));
C_pinv_tru__ = sqrtm(CC_inv_tru__);
%%%%%%%%;

rng(rseed);
X__ = zeros(n_var,n_T);
if  isempty(X_0_); X__(:,1+0) = randn(n_var,1); end;
if ~isempty(X_0_); X__(:,1+0) = reshape(X_0_,[n_var,1]); end;
Y__ = zeros(n_var,n_T);
Y__(:,1+0) = X__(:,1+0) + C_pinv_tru__*randn(n_var,1);
age_ = zeros(n_T,1);
dt_ = zeros(n_T-1,1);
R_ = zeros(n_T-1,1);
for nT=1:n_T-1;
dt=-dt_avg * log(1-rand());
dt_(1+nT-1) = dt;
age_(1+nT) = age_(1+nT-1) + dt;
X_pre_ = X__(:,1+nT-1);
tmp_exp_Adt__ = Psi__ * diag(exp(+Lambda_*dt)) * Psi_inv__;
X_pos_ = real( -A_inv_a_ + tmp_exp_Adt__ * ( (X_pre_ + A_inv_a_) + SDE_sample_int_expnAsBdW_1(dt,Psi__,Lambda_,Psi_inv__,B_pinv_tru__,1e-2)) );
X__(:,1+nT+0) = X_pos_;
Y__(:,1+nT+0) = X_pos_ + C_pinv_tru__*randn(n_var,1);
R_(1+nT-1) = fnorm(Y__(:,1+nT+0) - Y__(:,1+nT-1))/fnorm(Y__(:,1+nT-1));
end;%for nT=1:n_T-1;
R_avg = mean(R_);
