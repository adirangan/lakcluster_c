function output = SDE_sample_int_expnAsBdW_1(t,Psi_A__,Lambda_A_,Psi_A_inv__,B_inv__,tolerance_SDE);
% samples the function: ;
% \int_{s=0}^{s=t} exp(-A__ * s) * B_inv__ * dW(s) ;
% = \int_{s=0}^{s=t} Psi_A__ * exp(-diag(Lambda_A_) * s) * Psi_A_inv__ * B_inv__ * dW(s) ;
% where Psi_A__ * diag(Lambda_A_) * Psi_A_inv__ is an eigendecomposition of A__. ;
% Discretizes the total time t into increments dt that are sufficiently small that: ;
%  | exp(-Lambda_A_*dt) - (1-Lambda_A_*dt) | ;
% is (uniformly) less than tolerance_SDE. ;
% This is done simply by setting dt = sqrt(2*tolerance_SDE)/max(abs(Lambda_A_)). ;
% Consequently: ;
% |exp(-Lambda_A_ * dt) - (1-Lambda_A_*dt)| ;
% \approx | 1 - Lambda_A_*dt + 0.5 * (Lambda_A_*dt)^.2 - (1-Lambda_A_*dt) | ;
% == 0.5 * (Lambda_A_*dt).^2 ;
% <= 0.5 * (max(|Lambda_A_|)*dt).^2 ;
% == 0.5 * 2 * tolerance_SDE ;

n_var = numel(Lambda_A_);

Lambda_A_max = max(abs(Lambda_A_));
if (Lambda_A_max<=0); dt_0in = t; end;
if (Lambda_A_max> 0); dt_0in = sqrt(2*tolerance_SDE)/Lambda_A_max; end;
n_dt = max(2,ceil(t/dt_0in)); %<-- take at least 2 steps. ;
dt = t/n_dt;

S__ = Psi_A_inv__*B_inv__ ;
E__ = diag(exp(-Lambda_A_*dt));

output = zeros(n_var,1);
for ndt=0:n_dt-1;
output = output + S__ * randn(n_var,1) * sqrt(dt);
S__ = E__ * S__;
end;%for ndt=0:n_dt-1;
output = Psi_A__ * output;


