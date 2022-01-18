function [lpv,lP_0,flag_method,cap_,cup_] = label_to_label_enrichment_4(label_A_,label_B_);
% Calculates enrichment p-value. ;
% Assumes numeric labels, although labels do not need to be sequential. ;
% lpv = log of p-value (estimated). ;
% lP_0 = log of probability of achieving observed intersections. ;
%        Note that, up to a constant, this is the same as the entropy of the observed intersections. ;
% flag_method = 0 for trivial case, 
% flag_method = 1 for monte-carlo, 
% flag_method = 2 for monte-carlo + importance-sampling ;
% flag_method = 3 for simple lP0 (i.e., when importance-sampling fails). ;
% cap_ = set of intersections. ;
% cup_ = set of unions. ;
%%%%%%%%;

if (nargin<1);
disp(sprintf(' %% [testing label_to_label_enrichment_4]'));
setup_OptiPlex;
fname_mat = sprintf('~/dir_bcc/dir_jamison/dir_mat/label_to_label_enrichment_4.mat');
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
n_u = 36;
A0 = floor(1*n_u/9);
A1 = floor(2*n_u/9);
A2 = n_u - A0 - A1;
B0 = floor(1*n_u/6);
B1 = n_u - B0;
nlplim_ = [0,27];
c__ = colormap(colormap_nlpv); n_c = size(c__,1);
lP0_base = sum(gammaln(1+[A0;A1;A2])) + sum(gammaln(1+[B0;B1])) - gammaln(1+n_u);
nlP0__ = -Inf*ones(1+A0,1+A1);
x__ = zeros(4,(1+A0)*(1+A1));
y__ = zeros(4,(1+A0)*(1+A1));
c___ = zeros(1,(1+A0)*(1+A1),3);
na=0;
for nC10=0:A1;
for nC00=0:A0;
nC01 = A0 - nC00;
nC11 = A1 - nC10;
nC20 = B0 - nC00 - nC10;
nC21 = B1 - nC01 - nC11;
C__ = [ nC00 , nC01 ; nC10 , nC11 ; nC20 , nC21 ];
nlP0 = +Inf;
if ( (nC00>=0) & (nC01>=0) & (nC10>=0) & (nC11>=0) & (nC20>=0) & (nC21>=0) );
lP0 = lP0_base - sum(gammaln(1+C__),'all'); nlP0 = -lP0;
end;%if ( (nC00>=0) & (nC01>=0) & (nC10>=0) & (nC11>=0) & (nC20>=0) & (nC21>=0) );
nlP0__(1+nC00,1+nC10) = nlP0;
x__(:,1+na) = -0.5 + nC00 + [0;0;1;1] ;
y__(:,1+na) = -0.5 + nC10 + [0;1;1;0] ;
nc = max(0,min(n_c-1,floor(n_c*(nlP0 - min(nlplim_))/diff(nlplim_))));
c___(1,1+na,:) = c__(1+nc,:);
na=na+1;
end;%for nC00=0:A0;
end;%for nC10=0:A1;
%%%%%%%%;
n_fine = 4;
A0_fine_ = linspace(0,A0,1+n_fine*A0);
A1_fine_ = linspace(0,A1,1+n_fine*A1);
nlP0_fine__ = -Inf*ones(1+n_fine*A0,1+n_fine*A1);
na=0;
for nA0=0:numel(A0_fine_)-1; for nA1=0:numel(A1_fine_)-1;
nC00 = A0_fine_(1+nA0);
nC10 = A1_fine_(1+nA1);
nC01 = A0 - nC00;
nC11 = A1 - nC10;
nC20 = B0 - nC00 - nC10;
nC21 = B1 - nC01 - nC11;
C__ = [ nC00 , nC01 ; nC10 , nC11 ; nC20 , nC21 ];
nlP0 = +Inf;
if ( (nC00>=0) & (nC01>=0) & (nC10>=0) & (nC11>=0) & (nC20>=0) & (nC21>=0) );
lP0 = lP0_base - sum(gammaln(1+C__),'all'); nlP0 = -lP0;
end;%if ( (nC00>=0) & (nC01>=0) & (nC10>=0) & (nC11>=0) & (nC20>=0) & (nC21>=0) );
nlP0_fine__(1+nA0,1+nA1) = nlP0;
na=na+1;
end;end;%for nA0=0:numel(A0_fine_)-1; for nA1=0:numel(A1_fine_)-1;
%%%%%%%%;
tmp_index_ = efind(~isfinite(nlP0__));
c___(1,1+tmp_index_,:) = repmat([0.75,0.90,0.98],[numel(tmp_index_),1]);
%%%%%%%%;
n_label_A = 3; n_label_A_ = [A0;A1;A2]; num_A_ = zeros(n_u,1); num_A_(1:A0)=1+0; num_A_(A0+[1:A1])=1+1; num_A_(A0+A1+[1:A2])=1+2;
n_label_B = 2; n_label_B_ = [B0;B1]; num_B_ = zeros(n_u,1); num_B_(1:B0)=1+0; num_B_(B0+[1:B1])=1+1;
n_iteration = 1024; n_T = 32; ctol = 1e-3; dt=0.1; n_timestep = ceil(n_T/dt);
rng(0);
[lN__,lP__,cap_abst____] = ...
label_to_label_enrichment_sample_path_1( ...
 n_label_A ...
,n_label_A_ ...
,num_A_ ...
,n_label_B ...
,n_label_B_ ...
,num_B_ ...
,lP0_base ...
,n_iteration ...
,n_timestep ...
,dt ...
,ctol ...
);
sample_C00_tn__ = permute(reshape(cap_abst____(1+0,1+0,:,:),[n_iteration,n_timestep]),[2,1]);
sample_C10_tn__ = permute(reshape(cap_abst____(1+1,1+0,:,:),[n_iteration,n_timestep]),[2,1]);
sample_C01_tn__ = A0 - sample_C00_tn__;
sample_C11_tn__ = A1 - sample_C10_tn__;
sample_C20_tn__ = B0 - sample_C00_tn__ - sample_C10_tn__;
sample_C21_tn__ = B1 - sample_C01_tn__ - sample_C11_tn__;
save(fname_mat);
end;%if (~exist(fname_mat,'file'));
load(fname_mat);
%%%%%%%%;
[u_nlP0_,index_u_from_array_,index_array_from_u_] = unique(nlP0__(:));
index_u_from_array_ = index_u_from_array_ - 1; index_array_from_u_ = index_array_from_u_ - 1;
tmp_index_ = efind(isfinite(u_nlP0_));
u_nlP0_ = u_nlP0_(1+tmp_index_);
index_u_from_array_ = index_u_from_array_(1+tmp_index_);
index_array_from_u_ = index_array_from_u_(1+tmp_index_);
n_u_nlP0 = numel(u_nlP0_);
nlpv_true_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
nlP0 = u_nlP0_(1+nu_nlP0);
tmp_index_ = efind(nlP0__>=nlP0);
nlpv_true = -log(sum(exp(-nlP0__(1+tmp_index_))));
nlpv_true_(1+nu_nlP0) = nlpv_true;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
nlpv_mc00_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
nlP0 = u_nlP0_(1+nu_nlP0);
ntimestep=0;
tmp_lP_ = lP__(:,1+ntimestep);
tmp_lP_ij_ = find(tmp_lP_<=-nlP0);
if ( isempty(tmp_lP_ij_)); tmp_lpv = -Inf; end;
if (~isempty(tmp_lP_ij_)); tmp_lpv = log(numel(tmp_lP_ij_)/n_iteration); end;
nlpv_mc00_(1+nu_nlP0) = -tmp_lpv;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
nlpv_mcis_ = zeros(n_u_nlP0,1);
nlpv_mcis__ = zeros(n_timestep,n_u_nlP0);
for nu_nlP0=0:n_u_nlP0-1;
nlP0 = u_nlP0_(1+nu_nlP0);
tmp_lpv_ = zeros(n_timestep,1);
for ntimestep=0:n_timestep-1;
tmp_lP_ = lP__(:,1+ntimestep);
tmp_lP_ij_ = find(tmp_lP_<=-nlP0);
if ( isempty(tmp_lP_ij_)); tmp_lpv_(1+ntimestep) = -Inf; end;
if (~isempty(tmp_lP_ij_)); 
tmp_lw_ = tmp_lP_ - lP__(:,1); %<-- weight should be P(shifted)/P(original) = exp(lP_(ntimestep))/exp(lP_(1)). ;
tmp_lw_max = max(tmp_lw_); %<-- factor out largest exponential. ;
tmp_lw_adj_ = tmp_lw_ - tmp_lw_max; %<-- all adjusted weights are smaller than one. ;
tmp_lpv_(1+ntimestep) = tmp_lw_max + log(sum(exp(tmp_lw_adj_(tmp_lP_ij_)))/n_iteration) ;
end;%if (~isempty(tmp_lP_ij_)); 
end;%for ntimestep=0:n_timestep-1;
tmp_lpv = max(tmp_lpv_);
nlpv_mcis__(:,1+nu_nlP0) = -tmp_lpv_;
nlpv_mcis_(1+nu_nlP0) = -tmp_lpv;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;

%%%%%%%%;
% Now calculating gradient and hessian. ;
%%%%%%%%;
F_nlP0 = @(C__) ...
+ sum(gammaln(1+C__(:))) ...
+ gammaln(1+sum(C__,'all')) ...
- sum(gammaln(1+sum(C__,1))) ...
- sum(gammaln(1+sum(C__,2))) ...
;
D_nlP0 = @(C__) ...
+ dgammaln(1+C__) ...
+ dgammaln(1+sum(C__,'all')) ...
- repmat(dgammaln(1+sum(C__,1)),[size(C__,1),1]) ...
- repmat(dgammaln(1+sum(C__,2)),[1,size(C__,2)]) ...
;
tmp_size_ = size(C__);
tmp_C__ = 13 + rand(tmp_size_);
tmp_dC__ = rand(tmp_size_)*1e-3;
D0_ = (F_nlP0(tmp_C__ + tmp_dC__) - F_nlP0(tmp_C__ - tmp_dC__))/2;
D1_ = sum(D_nlP0(tmp_C__).*tmp_dC__,'all');
disp(sprintf(' %% error D0_ vs D1_ = %0.16f',fnorm(D0_-D1_)/fnorm(D0_)));
%%%%%%%%;
D_nlP0 = @(C__) ...
+ dgammaln(1+C__) ...
+ dgammaln(1+sum(C__,'all')) ...
- repmat(dgammaln(1+sum(C__,1)),[size(C__,1),1]) ...
- repmat(dgammaln(1+sum(C__,2)),[1,size(C__,2)]) ...
;
H_nlP0 = @(C__) ...
+ diag(ddgammaln(1+C__(:))) ...
+ ddgammaln(1+sum(C__,'all')) ...
- reshape(permute(repmat(diag(ddgammaln(1+sum(C__,1))),[1,1,size(C__,1),size(C__,1)]),[3,1,4,2]),[numel(C__),numel(C__)]) ...
- reshape(permute(repmat(diag(ddgammaln(1+sum(C__,2))),[1,1,size(C__,2),size(C__,2)]),[1,3,2,4]),[numel(C__),numel(C__)]) ...
;
tmp_size_ = size(C__);
tmp_C__ = 13 + rand(tmp_size_);
tmp_dC__ = rand(tmp_size_)*1e-3;
DD0_ = (D_nlP0(tmp_C__ + tmp_dC__) - D_nlP0(tmp_C__ - tmp_dC__))/2;
DD1_ = reshape(H_nlP0(tmp_C__)*tmp_dC__(:),tmp_size_);
disp(sprintf(' %% error DD0_ vs DD1_ = %0.16f',fnorm(DD0_-DD1_)/fnorm(DD0_)));
%%%%%%%%;
% Now build basis for admissible space. ;
%%%%%%%%;
n_C = n_label_A*n_label_B;
Q_off__ = zeros(n_C,n_label_A+n_label_B);
nl=0;
for nlabel_A=0:n_label_A-1;
tmp_ = zeros(n_label_A,n_label_B); tmp_(1+nlabel_A,:) = 1;
Q_off__(:,1+nl) = tmp_(:);
nl=nl+1;
end;%for nlabel_A=0:n_label_A-1;
for nlabel_B=0:n_label_B-1;
tmp_ = zeros(n_label_A,n_label_B); tmp_(:,1+nlabel_B) = 1;
Q_off__(:,1+nl) = tmp_(:);
nl=nl+1;
end;%for nlabel_B=0:n_label_B-1;
[tmp_Q_,tmp_R_] = qr([Q_off__,randn(n_C,n_C - n_label_A - n_label_B)]);
Q_off__ = tmp_Q_(:,1:n_label_A+n_label_B-1);
Q_0on__ = tmp_Q_(:,n_label_A+n_label_B:n_C);
n_d = size(Q_0on__,2);
%%%%%%%%;
% Now define critical point and measure concavity and build quadratic-approximation. ;
%%%%%%%%;
C_crit__ = n_u * (n_label_A_/n_u)*(transpose(n_label_B_)/n_u);
F_nlP0_crit = F_nlP0(C_crit__);
D_nlP0_crit_ = D_nlP0(C_crit__);
H_nlP0_crit__ = H_nlP0(C_crit__);
F_nlP0_quad = @(C__) F_nlP0_crit + 0.5*reshape(C__-C_crit__,[1,n_C])*H_nlP0_crit__*reshape(C__-C_crit__,[n_C,1]);
%%%%%%%%;
nlP0_quad_fine__ = -Inf*ones(1+n_fine*A0,1+n_fine*A1);
na=0;
for nA0=0:numel(A0_fine_)-1; for nA1=0:numel(A1_fine_)-1;
nC00 = A0_fine_(1+nA0);
nC10 = A1_fine_(1+nA1);
nC01 = A0 - nC00;
nC11 = A1 - nC10;
nC20 = B0 - nC00 - nC10;
nC21 = B1 - nC01 - nC11;
C__ = [ nC00 , nC01 ; nC10 , nC11 ; nC20 , nC21 ];
nlP0_quad_fine__(1+nA0,1+nA1) = F_nlP0_quad(C__);
na=na+1;
end;end;%for nA0=0:numel(A0_fine_)-1; for nA1=0:numel(A1_fine_)-1;
%%%%%%%%;
% Now define asymptotic p-value. ;
% Using the quadratic-approximation, we expect the nlP0 to act like: ;
% F_nlP0_quad = F_nlP0_crit + 0.5 * transpose(dC_) * H__ * dC_ ;
% with dC limited to the span of Q__ := Q_0on__. ;
% Therefore, we can write: ;
% F_nlP0_quad = F_nlP0_crit + 0.5 * transpose(dC_) * [Q__*transpose(Q__)] * H__ * [Q__*transpose(Q__)] * dC_ ;
% or: ;
% F_nlP0_quad = F_nlP0_crit + 0.5 * transpose(QdC_) * QHQ__ * QdC_ ;
% where QdC_ = transpose(Q__)*dC_ and QHQ__ = transpose(Q__)*H__*Q__. ;
% This is a simple anisotropic-gaussian with principal-variances given by the inverse-eigenvalues of QHQ__. ;
% The normalizing constant (i.e. integral) of this distribution is: ;
% I_quad = exp(-F_nlP0_crit) * \Pi_{dimensions} sqrt(2*pi)*sigma. ;
%%%%%%%%;
QHQ__ = transpose(Q_0on__)*H_nlP0_crit__*Q_0on__;
sigmainv_quad_ = sqrt(eig(QHQ__));
sigma_quad_ = 1./sigmainv_quad_;
nl_I_quad = F_nlP0_crit - log(sqrt(2*pi))*n_d - sum(log(sigma_quad_));
%%%%%%%%;
% Now the quadratic-approximation to the distribution is given by: ;
% rho(QdC_) = exp(-nl_I_quad) * exp(-F_nlP0_crit) * exp(-0.5 * (transpose(QdC_) * QHQ__ * QdC_) ) ;
% which, after aligning QdC_ to the eigenvectors of QHQ__, is equivalent to: ;
% rho(Y_) = exp(+nl_I_quad) * exp(-F_nlP0_crit) * exp(-0.5 * (transpose(Y_) * diag(sigmainv_quad_.^2) * Y_) ). ;
%%%%;
% Now given a particular nlP0_data from the data, we could achieve that nlP0 for values of QdC_ such that: ;
% nlP0_data = F_nlP0_crit + 0.5 * transpose(QdC_) * QHQ__ * QdC_ ;
% or ;
% 2*(nlP0_data - F_nlP0_crit) = transpose(Y_) * diag(sigmainv_quad_.^2) * Y_, ;
% which implies that any particular coordinate Y_(ny) should be at least: ;
% Y_(ny) = R*sigma_quad(ny), ;
% with R = sqrt(2*(nlP0_data-F_nlP0_crit)) ;
%%%%;
% The integral we are interested in is: ;
% p(lP0) = \int_{\omega(C)} rho(Y_) dY_, ;
% where \omega(C) is the region where nlP0>=nlP0_data. ;
% Now making the transformation X_ = Y_./sigma_quad_, we have: ;
% p(lP0) = prod(sigma_quad_)*\int_{diag(sigmainv_quad_)*\omega(C)} rho(X_) dX_, ;
%        = prod(sigma_quad_)*\int_{B(R)} exp(+nl_I_quad) * exp(-F_nlP0_crit) * exp(-0.5 * (transpose(X_) * X_) ) dX ;
%        = prod(sigma_quad_)*exp(+nl_I_quad)*exp(-F_nlP0_crit) * \int_{B(R)}  * exp(-0.5 * (transpose(X_) * X_) ) dX ;
%        = 1/sqrt(2*pi)^n_d * [Surface of sphere in n_d-dimensions] * \int_{R}^{\infty}  * exp(-0.5 * r^2) r^(n_d-1)dr ;
% where B(R) is now the exterior of the sphere of radius R. ;
%%%%;
% The integral I(R;n_d) = [ \int_{R}^{\infty} exp(-0.5*r^2) r^(n_d-1) dr ] satisfies: ;
% I(R;1) = sqrt(pi/2) * erfc(R/sqrt(2));
% I(R;2) = exp(-0.5*R^2) ;
% I(R;n_d>=3) = exp(-0.5*R^2)*R^(n_d-2) + (n_d-2)*I(R;n_d-2). ;
%%%%;
% So putting this all together, we have: ;
% p(lP0) = I(R;n_d) ;
%%%%%%%%;
nlpv_quad_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
nlP0 = u_nlP0_(1+nu_nlP0);
tmp_R = sqrt(2*(nlP0 - F_nlP0_crit));
tmp_I = 1/sqrt(2*pi)^n_d * exp(-0.5*tmp_R^2);
nlpv_quad_(1+nu_nlP0) = 0.5*tmp_R^2;
end;%for nu_nlP0=0:n_u_nlP0-1;

figure(1);clf; figbig;
c__ = colormap(colormap_nlpv); n_c = size(c__,1);
niteration_use_ = [2,3,5,8];
markersize_use = 8;
c_lime__ = colormap_nlpv(64,0.25,0,2.0); %<-- lime green. ;
n_c_lime = size(c_lime__,1);
subplot(3,3,[1,4]);
hold on;
patch(x__,y__,c___,'EdgeColor','none');
for nt=1:n_timestep-1;
nc_lime = max(0,min(n_c_lime-1,floor(n_c_lime*nt/n_timestep)));
line(sample_C00_tn__(1+[nt-1,nt],1+niteration_use_),sample_C10_tn__(1+[nt-1,nt],1+niteration_use_),'LineStyle','-','Color','k','LineWidth',3);
line(sample_C00_tn__(1+[nt-1,nt],1+niteration_use_),sample_C10_tn__(1+[nt-1,nt],1+niteration_use_),'LineStyle','-','Color',c_lime__(1+nc_lime,:),'LineWidth',2);
end;%for nt=1:n_timestep-1;
for nt=round(linspace(0,n_timestep-1,n_timestep/16));
nc_lime = max(0,min(n_c_lime-1,floor(n_c_lime*nt/n_timestep)));
plot(sample_C00_tn__(1+nt,1+niteration_use_),sample_C10_tn__(1+nt,1+niteration_use_),'o','Color',c_lime__(1+nc_lime,:),'MarkerSize',markersize_use,'MarkerEdgeColor','k','MarkerFaceColor',c_lime__(1+nc_lime,:));
end;%for nt=round(linspace(0,n_timestep-1,n_timestep/16));
hold off;
xlabel('C[1,1]'); ylabel('C[2,1]');
set(gca,'XTick',[0:A0],'XTickLabel',[0:A0]);
set(gca,'YTick',[0:A1],'YTickLabel',[0:A1]);
set(gca,'TickLength',[0,0]);
xlim([-0.5+0,A0+0.5]);
ylim([-0.5+0,A1+0.5]);
%axis image;
title('-log(phat) values');
tmp_c_ = colorbar; set(tmp_c_,'YTick',[],'YTickLabel',[]);
%%%%%%%%;
subplot(3,3,[2,5]);
hold on;
patch( -0.5 + [0;0;A0+1;A0+1] ,  -0.5 + [0;A1+1;A1+1;0] , 0.95*[1,1,1] , 'EdgeColor','none' );
for nl0=0:A0+1; line(nl0-0.5+[0;0],-0.5+[0;A1+1],'Color',0.85*[1,1,1],'LineWidth',0.5); end;%for nl0=0:A0+1;
for nl1=0:A1+1; line(-0.5+[0;A0+1],nl1-0.5+[0;0],'Color',0.85*[1,1,1],'LineWidth',0.5); end;%for nl1=0:A1+1;
contour(A0_fine_,A1_fine_,transpose(nlP0_fine__),[0:1:27],'LineWidth',2); set(gca,'Ydir','normal');
for nt=1:n_timestep-1;
nc_lime = max(0,min(n_c_lime-1,floor(n_c_lime*nt/n_timestep)));
line(sample_C00_tn__(1+[nt-1,nt],1+niteration_use_),sample_C10_tn__(1+[nt-1,nt],1+niteration_use_),'LineStyle','-','Color','k','LineWidth',3);
line(sample_C00_tn__(1+[nt-1,nt],1+niteration_use_),sample_C10_tn__(1+[nt-1,nt],1+niteration_use_),'LineStyle','-','Color',c_lime__(1+nc_lime,:),'LineWidth',2);
end;%for nt=1:n_timestep-1;
for nt=round(linspace(0,n_timestep-1,n_timestep/16));
nc_lime = max(0,min(n_c_lime-1,floor(n_c_lime*nt/n_timestep)));
plot(sample_C00_tn__(1+nt,1+niteration_use_),sample_C10_tn__(1+nt,1+niteration_use_),'o','Color',c_lime__(1+nc_lime,:),'MarkerSize',markersize_use,'MarkerEdgeColor','k','MarkerFaceColor',c_lime__(1+nc_lime,:));
end;%for nt=round(linspace(0,n_timestep-1,n_timestep/16));
hold off;
xlabel('C[1,1]'); ylabel('C[2,1]');
set(gca,'XTick',[0:A0],'XTickLabel',[0:A0]);
set(gca,'YTick',[0:A1],'YTickLabel',[0:A1]);
set(gca,'TickLength',[0,0]);
xlim([-0.5+0,A0+0.5]);
ylim([-0.5+0,A1+0.5]);
%axis image;
title('-log(phat) contours');
tmp_c_ = colorbar; set(tmp_c_,'YTick',[],'YTickLabel',[]);
%%%%%%%%;
subplot(3,3,[3,6]);
hold on;
patch( -0.5 + [0;0;A0+1;A0+1] ,  -0.5 + [0;A1+1;A1+1;0] , 0.95*[1,1,1] , 'EdgeColor','none' );
for nl0=0:A0+1; line(nl0-0.5+[0;0],-0.5+[0;A1+1],'Color',0.85*[1,1,1],'LineWidth',0.5); end;%for nl0=0:A0+1;
for nl1=0:A1+1; line(-0.5+[0;A0+1],nl1-0.5+[0;0],'Color',0.85*[1,1,1],'LineWidth',0.5); end;%for nl1=0:A1+1;
contour(A0_fine_,A1_fine_,transpose(nlP0_quad_fine__),[0:1:27],'LineWidth',2); set(gca,'Ydir','normal');
for nt=1:n_timestep-1;
nc_lime = max(0,min(n_c_lime-1,floor(n_c_lime*nt/n_timestep)));
line(sample_C00_tn__(1+[nt-1,nt],1+niteration_use_),sample_C10_tn__(1+[nt-1,nt],1+niteration_use_),'LineStyle','-','Color','k','LineWidth',3);
line(sample_C00_tn__(1+[nt-1,nt],1+niteration_use_),sample_C10_tn__(1+[nt-1,nt],1+niteration_use_),'LineStyle','-','Color',c_lime__(1+nc_lime,:),'LineWidth',2);
end;%for nt=1:n_timestep-1;
for nt=round(linspace(0,n_timestep-1,n_timestep/16));
nc_lime = max(0,min(n_c_lime-1,floor(n_c_lime*nt/n_timestep)));
plot(sample_C00_tn__(1+nt,1+niteration_use_),sample_C10_tn__(1+nt,1+niteration_use_),'o','Color',c_lime__(1+nc_lime,:),'MarkerSize',markersize_use,'MarkerEdgeColor','k','MarkerFaceColor',c_lime__(1+nc_lime,:));
end;%for nt=round(linspace(0,n_timestep-1,n_timestep/16));
hold off;
xlabel('C[1,1]'); ylabel('C[2,1]');
set(gca,'XTick',[0:A0],'XTickLabel',[0:A0]);
set(gca,'YTick',[0:A1],'YTickLabel',[0:A1]);
set(gca,'TickLength',[0,0]);
xlim([-0.5+0,A0+0.5]);
ylim([-0.5+0,A1+0.5]);
%axis image;
title('contours of quadratic approximation');
tmp_c_ = colorbar; set(tmp_c_,'YTick',[],'YTickLabel',[]);
%%%%%%%%;
subplot(3,3,[7]);
nu_nlP0 = 20;
ylim_ = [7.5,10.5];
plot(1:n_timestep,min(max(ylim_),nlpv_mcis__(:,1+nu_nlP0)),'m.-');
xlim([1,n_timestep]);
ylim(ylim_);
xlabel('time');
ylabel('-log(p)');
set(gca,'XTick',0:25:325,'XTickLabel',0:25:325);
set(gca,'YTick',10:13,'YTickLabel',10:13);
title(sprintf('time-sweep for fixed -log(phat(C)) = %0.2f',u_nlP0_(1+nu_nlP0)));
subplot(3,3,[8]);
markersize_dot = 12; linewidth_use = 1;
ylim_ = [0,12];
hold on;
plot(u_nlP0_,nlpv_mcis_,'m.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
plot(u_nlP0_,nlpv_true_,'k.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
plot(u_nlP0_,nlpv_quad_,'r.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
plot(u_nlP0_,min(max(ylim_),nlpv_mc00_),'g.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
plot(u_nlP0_,u_nlP0_,'c.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
hold off;
legend({'mcis','true','quad','mc','phat'},'Location','SouthEast');
xlim([min(u_nlP0_),max(u_nlP0_)]);
ylim(ylim_);
ylabel('-log(p)');
xlabel('-log(phat)');
title(sprintf('-log(p)'));
subplot(3,3,[9]);
markersize_dot = 12; linewidth_use = 1;
ylim_ = [0,12]/2;
hold on;
plot(u_nlP0_,nlpv_mcis_,'m.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
plot(u_nlP0_,nlpv_true_,'k.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
plot(u_nlP0_,nlpv_quad_,'r.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
plot(u_nlP0_,min(max(ylim_),nlpv_mc00_),'g.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
plot(u_nlP0_,u_nlP0_,'c.-','LineWidth',linewidth_use,'MarkerSize',markersize_dot);
hold off;
legend({'mcis','true','quad','mc','phat'},'Location','SouthEast');
xlim([min(u_nlP0_),max(u_nlP0_)/1.5]);
ylim(ylim_);
ylabel('-log(p)');
xlabel('-log(phat)');
title(sprintf('-log(p)'));
%%%%%%%%;
figbig;
disp(sprintf(' %% writing %s',sprintf('~/dir_bcc/dir_jamison/dir_jpg/label_to_label_enrichment_cartoon_4.eps')));
print('-depsc',sprintf('~/dir_bcc/dir_jamison/dir_jpg/label_to_label_enrichment_cartoon_4.eps'));
print('-djpeg',sprintf('~/dir_bcc/dir_jamison/dir_jpg/label_to_label_enrichment_cartoon_4.jpg'));
close(gcf);

disp('returning'); return;
end;%if (nargin<1);

%%%%%%%%;
verbose=0;
if (verbose); disp(sprintf(' %% [entering label_to_label_enrichment_4]')); end;
lpv = 0; lP_0 = 0; flag_method = 0;
%%%%%%%%;
n_A = numel(label_A_);
u_label_A_ = unique(label_A_);
n_label_A = length(u_label_A_);
n_label_A_ = zeros(n_label_A,1);
for nlabel_A = 1:n_label_A;
ij_A_{nlabel_A} = find(label_A_ == u_label_A_(nlabel_A));
n_label_A_(nlabel_A) = numel(ij_A_{nlabel_A});
end;%for nlabel_A=1:n_label_A;
%%%%%%%%;
n_B = numel(label_B_);
u_label_B_ = unique(label_B_);
n_label_B = length(u_label_B_);
n_label_B_ = zeros(n_label_B,1);
for nlabel_B = 1:n_label_B;
ij_B_{nlabel_B} = find(label_B_ == u_label_B_(nlabel_B));
n_label_B_(nlabel_B) = numel(ij_B_{nlabel_B});
end;%for nlabel_B=1:n_label_B;
%%%%%%%%;
if (n_A~=n_B); disp(sprintf(' %% Warning! n_A %d n_B %d in label_to_label_enrichment',n_A,n_B)); end;
assert(n_A==n_B);
n_X = n_A;
%%%%%%%%;
cap_ = zeros(n_label_A,n_label_B);
cup_ = zeros(n_label_A,n_label_B);
for nlabel_A=1:n_label_A;
for nlabel_B=1:n_label_B;
cap_(nlabel_A,nlabel_B) = numel(intersect(ij_A_{nlabel_A},ij_B_{nlabel_B}));
cup_(nlabel_A,nlabel_B) = numel(union(ij_A_{nlabel_A},ij_B_{nlabel_B}));
end;%for nlabel_B=1:n_label_B;
end;%for nlabel_A=1:n_label_A;
%%%%%%%%;
% Given a list of specific sets from A, ;
% as well as the cardinality of the sets in B, ;
% the number of ways of finding the precise intersections listed in cap_ is: ;
% lN = sum(gammaln(1+n_label_A_)) - sum(gammaln(1+cap_(:))). ;
% Similarly, the probability of observing cap_ is: ;
% lP = sum(gammaln(1+n_label_A_)) - sum(gammaln(1+cap_(:))) + sum(gammaln(1+n_label_B_)) - gammaln(1+n_X). ;
%%%%%%%%;
num_A_ = zeros(n_A,1);
for nlabel_A=1:n_label_A;
num_A_(ij_A_{nlabel_A}) = nlabel_A;
end;%for nlabel_A=1:n_label_A;
num_B_ = zeros(n_B,1);
for nlabel_B=1:n_label_B;
num_B_(ij_B_{nlabel_B}) = nlabel_B;
end;%for nlabel_B=1:n_label_B;
lP_base = sum(gammaln(1+n_label_A_)) + sum(gammaln(1+n_label_B_)) - gammaln(1+n_X);
cap_0_ = sparse(num_A_,num_B_,1,n_label_A,n_label_B);
tmp_0 = sum(gammaln(1+cap_0_),'all');
lN_0 = +tmp_0;
lP_0 = lP_base - tmp_0;
flag_method = 0;
%%%%%%%%;
if (lP_0<log(0.50));
%%%%%%%%;
% First we check to see if the p-value can be approximated using monte-carlo. ;
%%%%%%%%;
n_iteration=1024;
lN_ = zeros(n_iteration,1);
lP_ = zeros(n_iteration,1);
cap__ = zeros(n_label_A,n_label_B,n_iteration);
for niteration=1:n_iteration;
tmp_cap_ = sparse(num_A_,num_B_(randperm(n_B)),1,n_label_A,n_label_B);
tmp = sum(gammaln(1+tmp_cap_),'all');
lN_(niteration) = +tmp;
lP_(niteration) = lP_base - tmp;
cap__(:,:,niteration) = tmp_cap_;
end;%for niteration=1:n_iteration;
clear tmp_cap_;
% compare mean(cap__,3) with (n_label_A_*transpose(n_label_B_)/n_X) ;
% try: ;
%{
  cap_avg__ = mean(cap__,3);
  cap_bar__ = (n_label_A_*transpose(n_label_B_)/n_X);
  disp(sprintf(' %% relative deviation: %0.16f',fnorm(cap_avg__-cap_bar__)/fnorm(cap_bar__)));
  %}
[~,u_cap_ij_] = unique(transpose(reshape(cap__,n_label_A*n_label_B,n_iteration)),'row');
clear cap__;
n_u = numel(u_cap_ij_);
lN_u_ = lN_(u_cap_ij_);
lP_u_ = lP_(u_cap_ij_);
%%%%%%%%;
% if lN_0 falls in the range observed, we can ;
% model the probability of observing lN_u (w.r.t. sampling label_B_) ;
% as a gaussian in lN_u. ;
%%%%%%%%;
if (lN_0<prctile(lN_u_,95));
if (verbose); disp(sprintf(' %% using monte-carlo (flag_method == 1)')); end;
n_iteration=1024*8;
lN_ = zeros(n_iteration,1);
lP_ = zeros(n_iteration,1);
cap__ = zeros(n_label_A,n_label_B,n_iteration);
for niteration=1:n_iteration;
tmp_cap_ = sparse(num_A_,num_B_(randperm(n_B)),1,n_label_A,n_label_B);
tmp = sum(gammaln(1+tmp_cap_),'all');
lN_(niteration) = +tmp;
lP_(niteration) = lP_base - tmp;
cap__(:,:,niteration) = tmp_cap_;
end;%for niteration=1:n_iteration;
[~,u_cap_ij_] = unique(transpose(reshape(cap__,n_label_A*n_label_B,n_iteration)),'row');
clear tmp_cap_; clear cap__;
n_u = numel(u_cap_ij_);
lN_u_ = lN_(u_cap_ij_);
lP_u_ = lP_(u_cap_ij_);
mu = mean(lN_u_); sg = max(1e-12,std(lN_u_,1));
flag_plot=0;
if flag_plot;
n_h = 33;
hlN_ = mu + sg*4.5*linspace(-1,1,n_h);
h_ = hist(lN_u_,hlN_); h_ = h_/sum(h_)/mean(diff(hlN_));
g_ = 1/sqrt(2*pi) / sg * exp(-(hlN_-mu).^2/(2*sg^2));
subplot(1,2,1); plot(hlN_,g_,'r-',hlN_,h_,'k.'); xlabel('lN_u'); ylabel('p'); title('p');
subplot(1,2,2); plot(hlN_,log(g_),'r-',hlN_,log(h_),'k.'); ylim([min(log(g_)),max(log(g_))]); xlabel('lN_u'); ylabel('log(p)'); title('log(p)');
end;%if flag_plot;
z = (lN_0-mu)/sg;
lpv = z_to_p_0(z); flag_method = 1;
 else;%end;%if (lN_0<prctile(lN_u_,95));
%%%%%%%%;
% Otherwise, use gradient-flow to shift observations. ;
%%%%%%%%;
if (verbose); disp(sprintf(' %% using gradient-flow to sample rare regions (flag_method == 2)')); end;
flag_method = 2;
%%%%%%%%;
tmp_n_iteration=4; tmp_n_timestep=min(max(n_label_A_),max(n_label_B_)); 
[ ...
 tmp_lN__ ...
,tmp_lP__ ...
,~ ...
] = ...
label_to_label_enrichment_sample_path_1( ...
 n_label_A ...
,n_label_A_ ...
,num_A_ ...
,n_label_B ...
,n_label_B_ ...
,num_B_ ...
,lP_base ...
,tmp_n_iteration ...
,tmp_n_timestep ...
,[] ...
,[] ...
);
lP_avg_ = mean(tmp_lP__,1); lP_ij_ = find(lP_avg_<lP_0);
if  isempty(lP_ij_); 
disp(sprintf(' %% Warning: Preliminary search insufficient to locate rare regions.'));
disp(sprintf(' %% Increasing tmp_n_timestep is unlikely to improve sampling.'));
disp(sprintf(' %% Instead, we approximate integral with lP_0 (flag_method == 3).'));
lpv = lP_0; flag_method = 3;
end;%if  isempty(lP_ij_); 
if ~isempty(lP_ij_);
n_timestep = ceil(1+1.25*min(lP_ij_));
if (verbose); disp(sprintf(' %% found rare region in %d --> %d steps',min(lP_ij_),n_timestep)); end;
n_iteration = 256; %<-- typically enough for 1 digit accuracy. ;
[ ...
 lN__ ...
,lP__ ...
,~ ...
] = ...
label_to_label_enrichment_sample_path_1( ...
 n_label_A ...
,n_label_A_ ...
,num_A_ ...
,n_label_B ...
,n_label_B_ ...
,num_B_ ...
,lP_base ...
,n_iteration ...
,n_timestep ...
,[] ...
,[] ...
);
lpv_ = zeros(n_timestep,1);
for ntimestep=1:n_timestep;
tmp_lN_ = lN__(:,ntimestep);
tmp_lP_ = lP__(:,ntimestep);
tmp_lN_ij_ = find(tmp_lN_>lN_0);
if ( isempty(tmp_lN_ij_)); lpv_(ntimestep) = -Inf; end;
if (~isempty(tmp_lN_ij_)); 
tmp_lw_ = tmp_lP_ - lP__(:,1); %<-- weight should be P(shifted)/P(original) = exp(lP_(ntimestep))/exp(lP_(1)). ;
tmp_lw_max = max(tmp_lw_); %<-- factor out largest exponential. ;
tmp_lw_adj_ = tmp_lw_ - tmp_lw_max; %<-- all adjusted weights are smaller than one. ;
lpv_(ntimestep) = tmp_lw_max + log(sum(exp(tmp_lw_adj_(tmp_lN_ij_)))/n_iteration) ;
if (verbose); disp(sprintf(' %% ntimestep %d/%d; fraction %d/%d; tmp_lw_max %0.2f lpv %0.2f',ntimestep,n_timestep,length(tmp_lN_ij_),n_iteration,tmp_lw_max,lpv_(ntimestep))); end;
end;%if (~isempty(tmp_lN_ij_)); 
end;%for ntimestep=1:n_timestep;
lpv = max(lpv_); flag_method = 2; %<-- Note that here we take the largest value under the presumption that it is better to overestimate the p-value. ;
end;%if ~isempty(lP_ij_);
%%%%%%%%;
end;%if (lN_0<prctile(lN_u_,95));
end;%if (lP_0<log(0.50));
if (verbose); disp(sprintf(' %% [finished label_to_label_enrichment_4]')); end;
