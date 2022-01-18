function [lpv,lP_0,flag_method,cap_,cup_] = label_to_label_enrichment_quad_4(label_A_,label_B_);
% Calculates enrichment p-value using quadratic approximation. ;
% Assumes numeric labels, although labels do not need to be sequential. ;
% lpv = log of p-value (estimated). ;
% lP_0 = log of probability of achieving observed intersections. ;
%        Note that, up to a constant, this is the same as the entropy of the observed intersections. ;
% flag_method = 4 for quadratic approximation. ;
% cap_ = set of intersections. ;
% cup_ = set of unions. ;
%%%%%%%%;

if (nargin<1);

disp(sprintf(' %% [testing label_to_label_enrichment_quad_4]'));
setup_OptiPlex;

F_nlP0 = @(C__) ...
+ sum(gammaln(1+C__(:))) ...
+ gammaln(1+sum(C__,'all')) ...
- sum(gammaln(1+sum(C__,1))) ...
- sum(gammaln(1+sum(C__,2))) ...
;

fname_mat = sprintf('~/dir_bcc/dir_jamison/dir_mat/label_to_label_enrichment_quad_4_2d.mat');
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% 2-d example. ;
%%%%%%%%;
n_u = 64;
n_label_A = 3;
n_label_B = 2;
A0 = floor(n_u/3);
A1 = A0+floor(n_u/9);
A2 = n_u - A0 - A1;
B0 = floor(n_u/2);
B1 = n_u - B0;
nlplim_ = [0,27];
lP0_base = sum(gammaln(1+[A0;A1;A2])) + sum(gammaln(1+[B0;B1])) - gammaln(1+n_u);
nlP0_form__ = -Inf*ones(1+A0,1+A1);
nlpv_quad__ = -Inf*ones(1+A0,1+A1);
na=0;
for nC10=0:A1;
for nC00=0:A0;
nC01 = A0 - nC00;
nC11 = A1 - nC10;
nC20 = B0 - nC00 - nC10;
nC21 = B1 - nC01 - nC11;
C__ = [ nC00 , nC01 ; nC10 , nC11 ; nC20 , nC21 ];
nlP0_form = +Inf;
if ( (nC00>=0) & (nC01>=0) & (nC10>=0) & (nC11>=0) & (nC20>=0) & (nC21>=0) );
lP0 = lP0_base - sum(gammaln(1+C__),'all'); nlP0_form = -lP0;
nlP0_form__(1+nC00,1+nC10) = nlP0_form;
tmp_label_A_ = zeros(n_u,1); tmp_label_B_ = zeros(n_u,1);
nu=0;
for nc=0:numel(C__(:))-1;
C = C__(1+nc);
tmp_label_A_(1+[nu:nu+C-1]) = mod(nc,n_label_A);
tmp_label_B_(1+[nu:nu+C-1]) = floor(nc/n_label_A);
nu=nu+C;
end;%for nc=0:numel(C__(:))-1;
assert(nu==n_u);
[nlpv_quad,nlP0_quad] = label_to_label_enrichment_quad_4(tmp_label_A_,tmp_label_B_); nlpv_quad = -nlpv_quad; nlP0_quad = -nlP0_quad;
nlP0_error = fnorm(nlP0_form-nlP0_quad)/fnorm(nlP0_form);
assert(nlP0_error<1e-3);
nlpv_quad__(1+nC00,1+nC10) = nlpv_quad;
end;%if ( (nC00>=0) & (nC01>=0) & (nC10>=0) & (nC11>=0) & (nC20>=0) & (nC21>=0) );
na=na+1;
end;%for nC00=0:A0;
end;%for nC10=0:A1;
%%%%%%%%;
n_label_A_ = [A0;A1;A2]; num_A_ = zeros(n_u,1); num_A_(1:A0)=1+0; num_A_(A0+[1:A1])=1+1; num_A_(A0+A1+[1:A2])=1+2;
n_label_B_ = [B0;B1]; num_B_ = zeros(n_u,1); num_B_(1:B0)=1+0; num_B_(B0+[1:B1])=1+1;
n_iteration = 512; n_T = 32; ctol = 1e-3; dt=0.1; n_timestep = ceil(n_T/dt);
rng(0);
[lN__,lP__] = ...
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
save(fname_mat);
end;%if (~exist(fname_mat,'file'));
load(fname_mat);
%%%%%%%%;
nlpv_form__ = -Inf*ones(1+A0,1+A1);
for nC10=0:A1; for nC00=0:A0;
nlP0_form = nlP0_form__(1+nC00,1+nC10);
tmp_index_ = efind(isfinite(nlP0_form__) & (nlP0_form__>=nlP0_form));
nlpv_form = -log(sum(exp(-nlP0_form__(1+tmp_index_))));
nlpv_form__(1+nC00,1+nC10) = nlpv_form;
end;end;%for nC10=0:A1; for nC00=0:A0;
%%%%%%%%;
[u_nlP0_,index_u_from_array_,index_array_from_u_] = unique(nlP0_form__(:));
index_u_from_array_ = index_u_from_array_ - 1; index_array_from_u_ = index_array_from_u_ - 1;
tmp_index_ = efind(isfinite(u_nlP0_));
u_nlP0_ = u_nlP0_(1+tmp_index_);
index_u_from_array_ = index_u_from_array_(1+tmp_index_);
index_array_from_u_ = index_array_from_u_(1+tmp_index_);
n_u_nlP0 = numel(u_nlP0_);
%%%%%%%%;
nlpv_form_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
nlP0_form = u_nlP0_(1+nu_nlP0);
tmp_index_ = efind(nlP0_form__>=nlP0_form);
nlpv_form = -log(sum(exp(-nlP0_form__(1+tmp_index_))));
nlpv_form_(1+nu_nlP0) = nlpv_form;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
nlpv_mc00_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
nlP0_form = u_nlP0_(1+nu_nlP0);
ntimestep=0;
tmp_lP_ = lP__(:,1+ntimestep);
tmp_lP_ij_ = find(tmp_lP_<=-nlP0_form);
if ( isempty(tmp_lP_ij_)); tmp_lpv = -Inf; end;
if (~isempty(tmp_lP_ij_)); tmp_lpv = log(numel(tmp_lP_ij_)/n_iteration); end;
nlpv_mc00_(1+nu_nlP0) = -tmp_lpv;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
nlpv_mcis_ = zeros(n_u_nlP0,1);
nlpv_mcis__ = zeros(n_timestep,n_u_nlP0);
for nu_nlP0=0:n_u_nlP0-1;
nlP0_form = u_nlP0_(1+nu_nlP0);
tmp_lpv_ = zeros(n_timestep,1);
for ntimestep=0:n_timestep-1;
tmp_lP_ = lP__(:,1+ntimestep);
tmp_lP_ij_ = find(tmp_lP_<=-nlP0_form);
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
n_d = (n_label_A-1)*(n_label_B-1);
C_crit__ = n_u * (n_label_A_/n_u)*(transpose(n_label_B_)/n_u);
F_nlP0_crit = F_nlP0(C_crit__);
nlpv_quad_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
R = sqrt(2*max(0,u_nlP0_(1+nu_nlP0) - F_nlP0_crit));
d = (n_label_A-1) * (n_label_B-1);
lpv = lpv_erfc_iso_0(R,d);
nlpv_quad_(1+nu_nlP0) = -lpv;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
figure(1);clf; figbig; fignlpv; nlpvlim = [0,27]; markersize_use = 24;
%subplot(2,3,1); imagesc(nlP0_form__,nlpvlim); axisnotick; axis image; title('form');
%subplot(2,3,2); imagesc(nlpv_form__,nlpvlim); axisnotick; axis image; title('form');
%subplot(2,3,3); imagesc(nlpv_quad__,nlpvlim); axisnotick; axis image; title('quad');
%subplot(2,3,4); plot(u_nlP0_,nlpv_form__(1+index_u_from_array_),'ko',u_nlP0_,nlpv_quad__(1+index_u_from_array_),'ro'); title('2d');
subplot(1,3,1); plot(u_nlP0_,nlpv_form_,'k.',u_nlP0_,nlpv_quad_,'r.',u_nlP0_,nlpv_mc00_,'g.',u_nlP0_,nlpv_mcis_,'m.','MarkerSize',markersize_use); title('2d');

fname_mat = sprintf('~/dir_bcc/dir_jamison/dir_mat/label_to_label_enrichment_quad_4_3d.mat');
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% 3-d example. ;
%%%%%%%%;
n_u = 64;
n_label_A = 4;
n_label_B = 2;
A0 = floor(n_u/3);
A1 = floor(n_u/9);
A2 = floor(n_u/6);
A3 = n_u - A0 - A1 - A2;
B0 = floor(n_u/3);
B1 = n_u - B0;
nlplim_ = [0,27];
lP0_base = sum(gammaln(1+[A0;A1;A2;A3])) + sum(gammaln(1+[B0;B1])) - gammaln(1+n_u);
nlP0_form__ = -Inf*ones(1+A0,1+A1,1+A2);
nlpv_quad__ = -Inf*ones(1+A0,1+A1,1+A2);
na=0;
for nC20=0:A2;
for nC10=0:A1;
for nC00=0:A0;
nC01 = A0 - nC00;
nC11 = A1 - nC10;
nC21 = A2 - nC20;
nC30 = B0 - nC00 - nC10 - nC20;
nC31 = B1 - nC01 - nC11 - nC21;
C__ = [ nC00 , nC01 ; nC10 , nC11 ; nC20 , nC21 ; nC30 , nC31 ];
nlP0_form = +Inf;
if ( min(C__,[],'all')>=0 );
lP0 = lP0_base - sum(gammaln(1+C__),'all'); nlP0_form = -lP0;
nlP0_form__(1+nC00,1+nC10,1+nC20) = nlP0_form;
tmp_label_A_ = zeros(n_u,1); tmp_label_B_ = zeros(n_u,1);
nu=0;
for nc=0:numel(C__(:))-1;
C = C__(1+nc);
tmp_label_A_(1+[nu:nu+C-1]) = mod(nc,n_label_A);
tmp_label_B_(1+[nu:nu+C-1]) = floor(nc/n_label_A);
nu=nu+C;
end;%for nc=0:numel(C__(:))-1;
assert(nu==n_u);
[nlpv_quad,nlP0_quad] = label_to_label_enrichment_quad_4(tmp_label_A_,tmp_label_B_); nlpv_quad = -nlpv_quad; nlP0_quad = -nlP0_quad;
nlP0_error = fnorm(nlP0_form-nlP0_quad)/fnorm(nlP0_form);
assert(nlP0_error<1e-3);
nlpv_quad__(1+nC00,1+nC10,1+nC20) = nlpv_quad;
end;%if ( min(C__,[],'all')>=0 );
na=na+1;
end;%for nC00=0:A0;
end;%for nC10=0:A1;
end;%for nC20=0:A2;
%%%%%%%%;
n_label_A_ = [A0;A1;A2;A3]; num_A_ = zeros(n_u,1); num_A_(1:A0)=1+0; num_A_(A0+[1:A1])=1+1; num_A_(A0+A1+[1:A2])=1+2; num_A_(A0+A1+A2+[1:A3])=1+3;
n_label_B_ = [B0;B1]; num_B_ = zeros(n_u,1); num_B_(1:B0)=1+0; num_B_(B0+[1:B1])=1+1;
n_iteration = 256; n_T = 32; ctol = 1e-3; dt=0.1; n_timestep = ceil(n_T/dt);
rng(0);
[lN__,lP__] = ...
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
save(fname_mat);
end;%if (~exist(fname_mat,'file'));
load(fname_mat);
%%%%%%%%;
nlpv_form__ = -Inf*ones(1+A0,1+A1,1+A2);
for nC20=0:A2; for nC10=0:A1; for nC00=0:A0;
nlP0_form = nlP0_form__(1+nC00,1+nC10,1+nC20);
tmp_index_ = efind(isfinite(nlP0_form__) & (nlP0_form__>=nlP0_form));
nlpv_form = -log(sum(exp(-nlP0_form__(1+tmp_index_))));
nlpv_form__(1+nC00,1+nC10,1+nC20) = nlpv_form;
end;end;end;%for nC20=0:A2; for nC10=0:A1; for nC00=0:A0;
%%%%%%%%;
[u_nlP0_,index_u_from_array_,index_array_from_u_] = unique(nlP0_form__(:));
index_u_from_array_ = index_u_from_array_ - 1; index_array_from_u_ = index_array_from_u_ - 1;
tmp_index_ = efind(isfinite(u_nlP0_));
u_nlP0_ = u_nlP0_(1+tmp_index_);
index_u_from_array_ = index_u_from_array_(1+tmp_index_);
index_array_from_u_ = index_array_from_u_(1+tmp_index_);
n_u_nlP0 = numel(u_nlP0_);
%%%%%%%%;
nlpv_form_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
nlP0_form = u_nlP0_(1+nu_nlP0);
tmp_index_ = efind(nlP0_form__>=nlP0_form);
nlpv_form = -log(sum(exp(-nlP0_form__(1+tmp_index_))));
nlpv_form_(1+nu_nlP0) = nlpv_form;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
nlpv_mc00_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
nlP0_form = u_nlP0_(1+nu_nlP0);
ntimestep=0;
tmp_lP_ = lP__(:,1+ntimestep);
tmp_lP_ij_ = find(tmp_lP_<=-nlP0_form);
if ( isempty(tmp_lP_ij_)); tmp_lpv = -Inf; end;
if (~isempty(tmp_lP_ij_)); tmp_lpv = log(numel(tmp_lP_ij_)/n_iteration); end;
nlpv_mc00_(1+nu_nlP0) = -tmp_lpv;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
nlpv_mcis_ = zeros(n_u_nlP0,1);
nlpv_mcis__ = zeros(n_timestep,n_u_nlP0);
for nu_nlP0=0:n_u_nlP0-1;
nlP0_form = u_nlP0_(1+nu_nlP0);
tmp_lpv_ = zeros(n_timestep,1);
for ntimestep=0:n_timestep-1;
tmp_lP_ = lP__(:,1+ntimestep);
tmp_lP_ij_ = find(tmp_lP_<=-nlP0_form);
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
n_d = (n_label_A-1)*(n_label_B-1);
C_crit__ = n_u * (n_label_A_/n_u)*(transpose(n_label_B_)/n_u);
F_nlP0_crit = F_nlP0(C_crit__);
nlpv_quad_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
R = sqrt(2*max(0,u_nlP0_(1+nu_nlP0) - F_nlP0_crit));
d = (n_label_A-1) * (n_label_B-1);
lpv = lpv_erfc_iso_0(R,d);
nlpv_quad_(1+nu_nlP0) = -lpv;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
figure(1); 
subplot(1,3,2); plot(u_nlP0_,nlpv_form_,'k.',u_nlP0_,nlpv_quad_,'r.',u_nlP0_,nlpv_mc00_,'g.',u_nlP0_,nlpv_mcis_,'m.','MarkerSize',markersize_use); title('3d');

fname_mat = sprintf('~/dir_bcc/dir_jamison/dir_mat/label_to_label_enrichment_quad_4_4d.mat');
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% 4-d example. ;
%%%%%%%%;
n_u = 64;
n_label_A = 3;
n_label_B = 3;
A0 = floor(n_u/3);
A1 = floor(n_u/9);
A2 = n_u - A0 - A1 ;
B0 = floor(n_u/3);
B1 = floor(n_u/6);
B2 = n_u - B0 - B1;
nlplim_ = [0,27];
lP0_base = sum(gammaln(1+[A0;A1;A2])) + sum(gammaln(1+[B0;B1;B2])) - gammaln(1+n_u);
nlP0_form__ = -Inf*ones(1+A0,1+A1,1+A0,1+A1);
nlpv_quad__ = -Inf*ones(1+A0,1+A1,1+A0,1+A1);
na=0;
for nC11=0:A1;
for nC01=0:A0;
for nC10=0:A1;
for nC00=0:A0;
nC02 = A0 - nC00 - nC01;
nC12 = A1 - nC10 - nC11;
nC20 = B0 - nC00 - nC10;
nC21 = B1 - nC01 - nC11;
nC22 = B2 - nC02 - nC12;
C__ = [ nC00 , nC01 , nC02 ; nC10 , nC11 , nC12 ; nC20 , nC21 , nC22 ];
nlP0_form = +Inf;
if ( min(C__,[],'all')>=0 );
lP0 = lP0_base - sum(gammaln(1+C__),'all'); nlP0_form = -lP0;
nlP0_form__(1+nC00,1+nC10,1+nC01,1+nC11) = nlP0_form;
tmp_label_A_ = zeros(n_u,1); tmp_label_B_ = zeros(n_u,1);
nu=0;
for nc=0:numel(C__(:))-1;
C = C__(1+nc);
tmp_label_A_(1+[nu:nu+C-1]) = mod(nc,n_label_A);
tmp_label_B_(1+[nu:nu+C-1]) = floor(nc/n_label_A);
nu=nu+C;
end;%for nc=0:numel(C__(:))-1;
assert(nu==n_u);
[nlpv_quad,nlP0_quad] = label_to_label_enrichment_quad_4(tmp_label_A_,tmp_label_B_); nlpv_quad = -nlpv_quad; nlP0_quad = -nlP0_quad;
nlP0_error = fnorm(nlP0_form-nlP0_quad)/fnorm(nlP0_form);
assert(nlP0_error<1e-3);
nlpv_quad__(1+nC00,1+nC10,1+nC01,1+nC11) = nlpv_quad;
end;%if ( min(C__,[],'all')>=0 );
na=na+1;
end;%for nC00=0:A0;
end;%for nC10=0:A1;
end;%for nC01=0:A0;
end;%for nC11=0:A1;
%%%%%%%%;
n_label_A_ = [A0;A1;A2]; num_A_ = zeros(n_u,1); num_A_(1:A0)=1+0; num_A_(A0+[1:A1])=1+1; num_A_(A0+A1+[1:A2])=1+2;
n_label_B_ = [B0;B1;B2]; num_B_ = zeros(n_u,1); num_B_(1:B0)=1+0; num_B_(B0+[1:B1])=1+1; num_B_(B0+B1+[1:B2])=1+2;
n_iteration = 256; n_T = 32; ctol = 1e-3; dt=0.1; n_timestep = ceil(n_T/dt);
rng(0);
[lN__,lP__] = ...
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
save(fname_mat);
end;%if (~exist(fname_mat,'file'));
load(fname_mat);
%%%%%%%%;
nlpv_form__ = -Inf*ones(1+A0,1+A1,1+A0,1+A1);
for nC11=0:A1; for nC01=0:A0; for nC10=0:A1; for nC00=0:A0;
nlP0_form = nlP0_form__(1+nC00,1+nC10,1+nC01,1+nC11);
tmp_index_ = efind(isfinite(nlP0_form__) & (nlP0_form__>=nlP0_form));
nlpv_form = -log(sum(exp(-nlP0_form__(1+tmp_index_))));
nlpv_form__(1+nC00,1+nC10,1+nC01,1+nC11) = nlpv_form;
end;end;end;end;%for;
%%%%%%%%;
[u_nlP0_,index_u_from_array_,index_array_from_u_] = unique(nlP0_form__(:));
index_u_from_array_ = index_u_from_array_ - 1; index_array_from_u_ = index_array_from_u_ - 1;
tmp_index_ = efind(isfinite(u_nlP0_));
u_nlP0_ = u_nlP0_(1+tmp_index_);
index_u_from_array_ = index_u_from_array_(1+tmp_index_);
index_array_from_u_ = index_array_from_u_(1+tmp_index_);
n_u_nlP0 = numel(u_nlP0_);
%%%%%%%%;
nlpv_form_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
nlP0_form = u_nlP0_(1+nu_nlP0);
tmp_index_ = efind(nlP0_form__>=nlP0_form);
nlpv_form = -log(sum(exp(-nlP0_form__(1+tmp_index_))));
nlpv_form_(1+nu_nlP0) = nlpv_form;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
nlpv_mc00_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
nlP0_form = u_nlP0_(1+nu_nlP0);
ntimestep=0;
tmp_lP_ = lP__(:,1+ntimestep);
tmp_lP_ij_ = find(tmp_lP_<=-nlP0_form);
if ( isempty(tmp_lP_ij_)); tmp_lpv = -Inf; end;
if (~isempty(tmp_lP_ij_)); tmp_lpv = log(numel(tmp_lP_ij_)/n_iteration); end;
nlpv_mc00_(1+nu_nlP0) = -tmp_lpv;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
nlpv_mcis_ = zeros(n_u_nlP0,1);
nlpv_mcis__ = zeros(n_timestep,n_u_nlP0);
for nu_nlP0=0:n_u_nlP0-1;
nlP0_form = u_nlP0_(1+nu_nlP0);
tmp_lpv_ = zeros(n_timestep,1);
for ntimestep=0:n_timestep-1;
tmp_lP_ = lP__(:,1+ntimestep);
tmp_lP_ij_ = find(tmp_lP_<=-nlP0_form);
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
n_d = (n_label_A-1)*(n_label_B-1);
C_crit__ = n_u * (n_label_A_/n_u)*(transpose(n_label_B_)/n_u);
F_nlP0_crit = F_nlP0(C_crit__);
nlpv_quad_ = zeros(n_u_nlP0,1);
for nu_nlP0=0:n_u_nlP0-1;
R = sqrt(2*max(0,u_nlP0_(1+nu_nlP0) - F_nlP0_crit));
d = (n_label_A-1) * (n_label_B-1);
lpv = lpv_erfc_iso_0(R,d);
nlpv_quad_(1+nu_nlP0) = -lpv;
end;%for nu_nlP0=0:n_u_nlP0-1;
%%%%%%%%;
figure(1);
subplot(1,3,3); plot(u_nlP0_,nlpv_form_,'k.',u_nlP0_,nlpv_quad_,'r.',u_nlP0_,nlpv_mc00_,'g.',u_nlP0_,nlpv_mcis_,'m.','MarkerSize',markersize_use); title('4d');

print('-djpeg','~/dir_bcc/dir_jamison/dir_jpg/label_to_label_enrichment_quad_4.jpg');
print('-depsc','~/dir_bcc/dir_jamison/dir_jpg/label_to_label_enrichment_quad_4.eps');

disp('returning'); return;
end;% if (nargin<1);

%%%%%%%%;
verbose=0;
if (verbose); disp(sprintf(' %% [entering label_to_label_enrichment_quad_4]')); end;
lpv = 0; lP_0 = 0; flag_method = 4;
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
if ( (n_label_A<=1) | (n_label_B<=1) );
lpv = 0; lP_0 = 0;
if (verbose); disp(sprintf(' %% lP_0 %0.2f lpv %0.2f',lP_0,lpv)); end;
end;%if ( (n_label_A<=1) | (n_label_B<=1) );
if ( (n_label_A> 1) & (n_label_B> 1) );
%%%%%%%%;
% Given a list of specific sets from A, ;
% as well as the cardinality of the sets in B, ;
% the number of ways of finding the precise intersections listed in cap_ is: ;
% lN = sum(gammaln(1+n_label_A_)) - sum(gammaln(1+cap_(:))). ;
% Similarly, the probability of observing cap_ is: ;
% nlP = + sum(gammaln(1+cap_(:))) + gammaln(1+n_X) - sum(gammaln(1+n_label_B_)) - sum(gammaln(1+n_label_A_)). ;
%%%%%%%%;
F_nlP0 = @(C__) ...
+ sum(gammaln(1+C__(:))) ...
+ gammaln(1+sum(C__,'all')) ...
- sum(gammaln(1+sum(C__,1))) ...
- sum(gammaln(1+sum(C__,2))) ...
;
C_crit__ = n_X * (n_label_A_/n_X)*(transpose(n_label_B_)/n_X);
F_nlP0_crit = F_nlP0(C_crit__);
F_nlP0_data = F_nlP0(cap_);
%%%%%%%%;
% quadratic approximation. ;
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
% Note that [Surface of sphere in n_d-dimensions] = n_d*pi^(n_d/2) / gamma(1 + n_d/2);
%%%%;
% The integral I(R;n_d) = [ \int_{R}^{\infty} exp(-0.5*r^2) r^(n_d-1) dr ] satisfies: ;
% I(R;1) = sqrt(pi/2) * erfc(R/sqrt(2));
% I(R;2) = exp(-0.5*R^2) ;
% I(R;n_d>=3) = exp(-0.5*R^2)*R^(n_d-2) + (n_d-2)*I(R;n_d-2). ;
%%%%;
% So putting this all together, we have: ;
% p(lP0) = [ n_d*pi^(n_d/2) / gamma(1 + n_d/2) ] * [ 1/sqrt(2*pi)^n_d ] * I(R;n_d) ;
%%%%%%%%;
R = sqrt(2*max(0,F_nlP0_data - F_nlP0_crit));
d = (n_label_A-1) * (n_label_B-1);
lP_0 = -F_nlP0_data;
lpv = lpv_erfc_iso_1(R,d);
if (verbose); disp(sprintf(' %% R %0.2f d %0.2f lP_0 %0.2f lpv %0.2f',R,d,lP_0,lpv)); end;
end;%if ( (n_label_A> 1) & (n_label_B> 1) );
