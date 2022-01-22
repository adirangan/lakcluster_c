function [output,error_final,iteration_final] = adi_median(input_,tolerance_p,verbose);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% function output = adi_median(input_); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% This function reads in an array 'input_' of size d-by-N. ;
% This input_ is interpreted as a collection of N points in dimension d. ;
% The output is the median (in dimension d). ;
% The median is calculated by minimizing the L2-distance: ;
% E1 = \sum_{n} \| m - x_{n} \|, ;
% where m is the median, and x_{n} is the n^{th} point. ;
% ;
% Inputs: ;
% input_ (double) d-by-N input array. ;
% tolerance_p (double) error tolerance used for gradient (default 1e-9). ;
% verbose (integer) verbosity level (default 0). ;
% ;
% test by running with no arguments: ;
% adi_median();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

if (nargin<1);
rng(1);
verbose = 1;
tolerance_p = 1e-9;
input_ = rand(2,1024*8);
t_start = tic;
[output,error_final,iteration_final] = adi_median(input_,tolerance_p,verbose);
t_elapsed = toc(t_start);
disp(sprintf(' %% final error %0.16f; iteration %d; time elapsed %0.16f',error_final,iteration_final,t_elapsed));
disp('returning'); return;
end;%if (nargin<1);

if (nargin<3); verbose=0; end;
if (nargin<2); tolerance_p=1e-9; end;

m = mean(input_,2);
p = 1.0;
iteration_outer_max = 1024;
iteration_outer_cur = 1;
continue_outer_flag = 1;
while (continue_outer_flag);
iteration_inner_max = 1024;
iteration_inner_cur = 0;
continue_inner_flag = 1;
while (continue_inner_flag);
D1_ = D1(input_,m);
E1_ = E1(D1_);
Ep_ = Ep(p,E1_);
dEp_ = dEp(p,E1_,D1_);
ddEp_ = ddEp(p,E1_,D1_);
dm = ddEp_ \ dEp_ ;
m = m+dm;
iteration_inner_cur = iteration_inner_cur + 1;
continue_inner_flag = (iteration_inner_cur<iteration_inner_max & norm(dm)>tolerance_p & norm(dEp_)>tolerance_p);
end;%while (continue_inner_flag);
D1_ = D1(input_,m);
E1_ = E1(D1_);
dEp_target_ = dEp(0.5,E1_,D1_);
if (verbose); disp(sprintf(' %% p %0.2f, iteration_inner_cur %d, norm(dm) %0.16f, norm(dEp_) %0.16f norm(dEp_target_) %0.16f',p,iteration_inner_cur,norm(dm),norm(dEp_),norm(dEp_target_))); end;
p = 0.5 + (p-0.5)*0.5;
iteration_outer_cur = iteration_outer_cur + 1;
continue_outer_flag = (iteration_outer_cur<iteration_outer_max & norm(dEp_target_)>tolerance_p);
end;%while (continue_outer_flag);
error_final = norm(dEp_target_);
iteration_final = iteration_outer_cur;

output = m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function ddEp_ = ddEp(p,E1_,D1_);
[d,N] = size(D1_);
E1p1_ = repmat(E1_.^(p-1),d*d,1);
E1p2_ = repmat(E1_.^(p-2),d*d,1);
Id = repmat(reshape(eye(d),d*d,1),1,N);
DD = zeros(d*d,N);
for nd1=0:d-1; for nd2=0:d-1;
DD(1+nd1+nd2*d,:) = D1_(1+nd1,:).*D1_(1+nd2,:);
end;end;%for nd1=0:d-1; for nd2=0:d-1;
ddEp_ = reshape(sum(p.*(p-1).*E1p2_.*4.*DD + p.*E1p1_.*2.*Id,2),d,d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function dEp_ = dEp(p,E1_,D1_);
[d,N] = size(D1_);
dEp_ = sum(p.*repmat(E1_.^(p-1),d,1).*2.*D1_,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function Ep_ = Ep(p,E1_);
Ep_ = sum(E1_.^p,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function E1_ = E1(D1_);
E1_ = sum(D1_.^2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function D1_ = D1(input_,m);
[d,N] = size(input_);
D1_ = input_ - repmat(m,1,N);

