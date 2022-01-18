function [DDL_all_,USL_n_,VL_n_,SL_n_] = SVD_discat_0(A_n_,label_A_,n_rank,alpha_over_beta);
% returns supervised SVD results for categorical covariates listed in label_A_. ;
%%%%%%%%;
% Input: ;
% A_n_ : double array of size (M,N). ; Data array. ;
% label_A_: double vector of size M. ; categorical covariates, assumed numerical. ;
% n_rank: number of principal-components requested. ;
% alpha_over_beta: double. ; ratio of weighting factors intra-cluster/inter-cluster in loss-function. ;
%%%%%%%%;
% Output: ;
% DDL_all_: double array of size (N,N). ; Hessian of loss-function. ;
% USL_n_: double array of size (M,n_rank). ; Projected vectors associated with data-array: A_n_*V_n_. ;
% VL_n_: double array of size (N,n_rank). ; Projection vectors associated with data-array. ;
% SL_n_: double array of size (n_rank,n_rank). ; Diagonal matrix of eigenvalues. ;
%       Each eigenvalue is equal to twice the loss-function associated with the corresponding column of V_n_. ;
%%%%%%%%;

if (nargin<1);
disp(sprintf(' %% testing SVD_discat_0'));

verbose = 1;
k_p_r_max = 1.0;
k_eq_d = 1/(2*pi)*k_p_r_max;
[ ...
 n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,J_node_ ...
,J_weight_ ...
,J_chebfun_ ...
,J_polyval_ ...
] = ...
sample_sphere_7( ...
 verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
) ;
A_n_ = [ 3*k_c_0_all_ , 2*k_c_1_all_ , 1*k_c_2_all_ ];
label_A_ = zeros(size(A_n_,1),1);
label_A_(1+efind(k_c_2_all_>0.5*k_p_r_max)) = 1;
%%%%%%%%;
u_label_A_ = unique(label_A_);
n_label_A = numel(u_label_A_);
n_label_A_ = zeros(n_label_A,1);
index_label_A_ = cell(n_label_A,1);
for nlabel_A=0:n_label_A-1;
index_label_A_{1+nlabel_A} = efind(label_A_==u_label_A_(1+nlabel_A));
n_label_A_(1+nlabel_A) = numel(index_label_A_{1+nlabel_A});
end;%for nlabel_A=0:n_label_A-1;
%%%%%%%%;

[ ...
 ~ ...
,shell_azimu_b_all_ ...
,shell_polar_a_all_ ...
,~ ...
,shell_k_c_0_all_ ...
,shell_k_c_1_all_ ...
,shell_k_c_2_all_ ...
,~ ...
,~ ...
,~ ...
] = ...
sample_shell_5( ...
 k_p_r_max ...
,k_eq_d/4 ...
,'L' ...
) ;

v_ = randn(3,1);
for nd=0:3-1;
dv_ = zeros(3,1); dv_(1+nd)=1*1e-4;
L0_ = get_L_01(A_n_,index_label_A_{1+0},index_label_A_{1+1},v_ - dv_);
L0_bkp_ = get_L_01_bkp(A_n_,index_label_A_{1+0},index_label_A_{1+1},v_ - dv_);
disp(sprintf(' %% testing get_L_01: error %0.16f',fnorm(L0_ - L0_bkp_)/fnorm(L0_)));
L1_ = get_L_01(A_n_,index_label_A_{1+0},index_label_A_{1+1},v_ + dv_);
L1_bkp_ = get_L_01_bkp(A_n_,index_label_A_{1+0},index_label_A_{1+1},v_ + dv_);
disp(sprintf(' %% testing get_L_01: error %0.16f',fnorm(L1_ - L1_bkp_)/fnorm(L1_)));
DL_diff_(:,1+nd) = (L1_ - L0_)/2;
DL_form_(:,1+nd) = transpose(dv_)*get_DDL_01(A_n_,index_label_A_{1+0},index_label_A_{1+1})*v_;
end;%for nd=0:3-1;
disp(sprintf(' %% testing get_DDL_01: error %0.16f',fnorm(DL_diff_ - DL_form_)/fnorm(DL_diff_)));

flag_test=0;
if flag_test;
%%%%%%%%;
% visualize loss. ;
%%%%%%%%;
v__ = transpose([ shell_k_c_0_all_ , shell_k_c_1_all_ , shell_k_c_2_all_ ]);
L_all_ = get_L_all(A_n_,index_label_A_,v__);
L_lim_ = [ min(L_all_,[],'all') , max(L_all_,[],'all') ];
disp(sprintf(' %% L_lim_ %+0.2f %+0.2f',L_lim_));
figure(1);clf;figbig; c_ = colormap_beach();
imagesc_polar_a_azimu_b_0(shell_polar_a_all_,shell_azimu_b_all_,L_all_,L_lim_,c_,0);
%xlim([0,2*pi]);ylim([0,pi]);
xlim([-1,+1]); ylim([-1,+1]); zlim([-1,+1]); axis equal; axis vis3d; xlabel('x'); ylabel('y'); zlabel('z');
end;%if flag_test;

[~,USL_,VL_,SL_] = SVD_discat_0(A_n_,label_A_);

[UA_,SA_,VA_] = svds(A_n_,3);
USA_ = UA_*SA_;
figure(1);clf;figbig;figbeach();
subplot(1,2,1);
scatter(USA_(:,1),USA_(:,2),15,label_A_,'filled');
axis equal;
subplot(1,2,2);
scatter(USL_(:,1),USL_(:,2),15,label_A_,'filled');
axis equal;

disp('returning'); return;
end;%if (nargin<1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (nargin<3); n_rank = []; end;
if (nargin<4); alpha_over_beta = []; end;
if isempty(n_rank); n_rank = min(min(size(A_n_)),3); end;
if isempty(alpha_over_beta); alpha_over_beta = +1; end;

u_label_A_ = unique(label_A_);
n_label_A = numel(u_label_A_);
n_label_A_ = zeros(n_label_A,1);
index_label_A_ = cell(n_label_A,1);
for nlabel_A=0:n_label_A-1;
index_label_A_{1+nlabel_A} = efind(label_A_==u_label_A_(1+nlabel_A));
n_label_A_(1+nlabel_A) = numel(index_label_A_{1+nlabel_A});
end;%for nlabel_A=0:n_label_A-1;

DDL_all_ = get_DDL_all(A_n_,index_label_A_,alpha_over_beta);
if (nargout>1);
%%%%%%%%;
[VL_n_,SL_n_] = eigs(DDL_all_,n_rank,'largestreal');
USL_n_ = A_n_*VL_n_;
%%%%%%%%;
flag_check=0;
if flag_check;
for nrank=0:n_rank-1;
L = get_L_all(A_n_,index_label_A_,VL_n_(:,1+nrank),alpha_over_beta);
disp(sprintf(' %% nrank %d/%d eig vs 2*L %0.16f',nrank,n_rank,fnorm(SL_n_(1+nrank,1+nrank)-2*L)/fnorm(2*L)));
end;%for nrank=0:n_rank-1;
end;%if flag_check;
%%%%%%%%;
end;%if (nargout>1);

function output_ = get_DDL_all(A_n_,index_label_A_,alpha_over_beta);
% returns the Hessian of the loss-function: ;
% L = \sum_{index_0_} \sum_{index_1_} L_01 / Z01 * delta;
% where L_01 is defined below (for each pair of sets index_0_ and index_1_), ;
% and Z01 is a normalizing-factor, ;
% and delta:=-alpha if the sets are the same, and delta:=beta if the sets are different.; 
if (nargin<3); alpha_over_beta = +1; end;
[M,N] = size(A_n_);
output_ = zeros(N,N);
n_label_A = numel(index_label_A_);
for nlabel_A_0=0:n_label_A-1;
index_0_ = index_label_A_{1+nlabel_A_0};
for nlabel_A_1=0:n_label_A-1;
index_1_ = index_label_A_{1+nlabel_A_1};
if (nlabel_A_0==nlabel_A_1); delta = -alpha_over_beta; end;
if (nlabel_A_0~=nlabel_A_1); delta = 1.0; end;
output_ = output_ + delta * get_DDL_01(A_n_,index_0_,index_1_) / get_Z01(index_0_,index_1_);
end;%for nlabel_A_1=0:n_label_A-1;
end;%for nlabel_A_0=0:n_label_A-1;

function output = get_Z01(index_0_,index_1_);
% normalizing-factor for L_01. ;
n_0 = numel(index_0_);
n_1 = numel(index_1_);
%output = 1;;
output = 1/(n_0*n_1);

function output_ = get_DDL_01(A_n_,index_0_,index_1_);
% returns the Hessian of the loss-function: ;
% L_01 = \sum_{j0\in index_0_} \sum_{j1\in index_1_} \| A_n_(j0,:)*v - A_n_(j1,:)*v \|^{2} ;
% with respect to the vector v_. ;
% Note that the full gradient requires multiplying the output_ by v_ on the right. ;
n_0 = numel(index_0_);
n_1 = numel(index_1_);
A_t_ = transpose(A_n_);
e0_n_ = ones(n_0,1); e0_t_ = transpose(e0_n_);
e1_n_ = ones(n_1,1); e1_t_ = transpose(e1_n_);
output_ = 2*( ...
 + A_t_(:,1+index_0_) * ( n_1*A_n_(1+index_0_,:) - e0_n_*e1_t_*A_n_(1+index_1_,:) ) ...
 + A_t_(:,1+index_1_) * ( n_0*A_n_(1+index_1_,:) - e1_n_*e0_t_*A_n_(1+index_0_,:) ) ...
	      ) ;

function output_ = get_L_all(A_n_,index_label_A_,v_,alpha_over_beta);
% returns the loss-function: ;
% L = \sum_{index_0_} \sum_{index_1_} L_01 / Z01 * delta;
% where L_01 is defined below (for each pair of sets index_0_ and index_1_), ;
% and Z01 is a normalizing-factor, ;
% and delta:=-alpha if the sets are the same, and delta:=beta if the sets are different.; 
if (nargin<4); alpha_over_beta = +1; end;
[M,N] = size(A_n_);
output_ = zeros(1,size(v_,2));
n_label_A = numel(index_label_A_);
for nlabel_A_0=0:n_label_A-1;
index_0_ = index_label_A_{1+nlabel_A_0};
for nlabel_A_1=0:n_label_A-1;
index_1_ = index_label_A_{1+nlabel_A_1};
if (nlabel_A_0==nlabel_A_1); delta = -alpha_over_beta; end;
if (nlabel_A_0~=nlabel_A_1); delta = 1.0; end;
output_ = output_ + delta * get_L_01(A_n_,index_0_,index_1_,v_) / get_Z01(index_0_,index_1_);
end;%for nlabel_A_1=0:n_label_A-1;
end;%for nlabel_A_0=0:n_label_A-1;

function output = get_L_01(A_n_,index_0_,index_1_,v_);
% returns loss-function L_01 above. ;
n_0 = numel(index_0_);
n_1 = numel(index_1_);
e0_n_ = ones(n_0,1); e0_t_ = transpose(e0_n_);
e1_n_ = ones(n_1,1); e1_t_ = transpose(e1_n_);
A_n_0_v_ = A_n_(1+index_0_,:)*v_;
A_n_1_v_ = A_n_(1+index_1_,:)*v_;
output = n_1*e0_t_*A_n_0_v_.^2 + n_0*e1_t_*A_n_1_v_.^2 - 2*(e0_t_*A_n_0_v_).*(e1_t_*A_n_1_v_);

function output = get_L_01_bkp(A_n_,index_0_,index_1_,v_);
% returns loss-function L_01 above. ;
n_0 = numel(index_0_);
n_1 = numel(index_1_);
e0_n_ = ones(n_0,1); e0_t_ = transpose(e0_n_);
e1_n_ = ones(n_1,1); e1_t_ = transpose(e1_n_);
A_n_0_v_ = A_n_(1+index_0_,:)*v_;
A_n_1_v_ = A_n_(1+index_1_,:)*v_;
output = 0;
for n0=0:n_0-1;
for n1=0:n_1-1;
output = output + (A_n_0_v_(1+n0) - A_n_1_v_(1+n1)).^2;
end;%for n1=0:n_1-1;
end;%for n0=0:n_0-1;


