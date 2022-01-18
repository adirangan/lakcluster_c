function ...
[ ...
 parameter ...
,A_sqr_nlp___ ...
,A_sqr_bnlp___ ...
,A_sqr_hnlp___ ...
,A_sqr_jnlp___ ...
,A_sqr_snlp___ ...
,A_sqr_sbnlp___ ...
,A_sqr_shnlp___ ...
,A_sqr_sjnlp___ ...
] = ...
dolphin_holm_from_prm_0( ...
 parameter ...
,n_var ...
,mask__ ...
,n_shuffle ...
,A_sqr_prm___ ...
);

verbose=0;
if (verbose); disp(sprintf(' %% [entering dolphin_holm_from_prm_0]')); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_var=[]; end; na=na+1;
if (nargin<1+na); mask__=[]; end; na=na+1;
if (nargin<1+na); n_shuffle=[]; end; na=na+1;
if (nargin<1+na); A_sqr_prm___=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if isempty(mask__); mask__ = ones(n_var,n_var)-eye(n_var,n_var); end;

if ( numel(mask__==1) & (mask__=='A') ); mask__ = ones(n_var,n_var)-eye(n_var,n_var); n_bonferroni = n_var*(n_var-1); end;
if ( numel(mask__==1) & (mask__=='B') ); mask__ = ones(n_var,n_var)-eye(n_var,n_var); n_bonferroni = n_var*(n_var - 1)/2; end;

assert(size(mask__,1)==n_var);
assert(size(mask__,2)==n_var);
assert(size(A_sqr_prm___,1)==n_var);
assert(size(A_sqr_prm___,2)==n_var);
assert(size(A_sqr_prm___,3)==1+n_shuffle);

index_retain_ = efind(mask__);
n_retain = numel(index_retain_);
A_sqr_prm__ = reshape(A_sqr_prm___,[n_var^2,n_shuffle+1]);
A_prm__ = A_sqr_prm__(1+index_retain_,:);
A_avg_ = mean(A_prm__(:,2:end),2); A_std_ = std(A_prm__(:,2:end),1,2);
A_zsc__ = (A_prm__ - A_avg_)./A_std_;
A_nlp_opt__ = -z_to_p_twosided_0(A_zsc__);
%%%%%%%%;

%%%%%%%%;
% faster than loop below. ;
%%%%%%%%;
tmp_t = tic();
tmp__ = transpose(A_prm__);
[~,tmp_index_lo__] = sort(+tmp__,1,'ascend'); tmp_index_lo__ = tmp_index_lo__-1;
[~,tmp_index_lo__] = sort(tmp_index_lo__,1); tmp_index_lo__ = tmp_index_lo__-1;
tmp_index_lo__ = tmp_index_lo__ - bsxfun(@gt,tmp_index_lo__,reshape(tmp_index_lo__(1+0,:),[1,n_retain]));
A_nlo__ = tmp_index_lo__; %<-- number lower (not counting self). ;
[~,tmp_index_hi__] = sort(-tmp__,1,'ascend'); tmp_index_hi__ = tmp_index_hi__-1;
[~,tmp_index_hi__] = sort(tmp_index_hi__,1); tmp_index_hi__ = tmp_index_hi__-1;
tmp_index_hi__ = tmp_index_hi__ - bsxfun(@gt,tmp_index_hi__,reshape(tmp_index_hi__(1+0,:),[1,n_retain]));
A_nhi__ = tmp_index_hi__; %<-- number higher (not counting self). ;
tmp_t = toc(tmp_t);
clear tmp__;
A_nlo__ = transpose(A_nlo__);
A_nhi__ = transpose(A_nhi__);
if (verbose>1); disp(sprintf(' %% A_nlo__ and A_nhi__: %0.3fs',tmp_t)); end;
%{
%%%%%%%%;
% slower than vectorized version above. ;
%%%%%%%%;
tmp_t = tic();
A_nhi__ = zeros(size(A_prm__));
A_nlo__ = zeros(size(A_prm__));
for nretain=0:n_retain-1;
tmp_A_ = +A_prm__(1+nretain,:);
[~,tmp_index_] = sort(tmp_A_, 'ascend'); tmp_index_ = tmp_index_-1;
[~,tmp_index_] = sort(tmp_index_); tmp_index_ = tmp_index_-1;
tmp_tmp_index_ = efind(tmp_index_> tmp_index_(1+0));
tmp_index_(1+tmp_tmp_index_) = tmp_index_(1+tmp_tmp_index_) - 1;
A_nlo__(1+nretain,:) = tmp_index_; %<-- number lower (not counting self). ;
tmp_A_ = -A_prm__(1+nretain,:);
[~,tmp_index_] = sort(tmp_A_, 'ascend'); tmp_index_ = tmp_index_-1;
[~,tmp_index_] = sort(tmp_index_); tmp_index_ = tmp_index_-1;
tmp_tmp_index_ = efind(tmp_index_> tmp_index_(1+0));
tmp_index_(1+tmp_tmp_index_) = tmp_index_(1+tmp_tmp_index_) - 1;
A_nhi__(1+nretain,:) = tmp_index_; %<-- number higher (not counting self). ;
end;%for nretain=0:n_retain-1;
tmp_t = toc(tmp_t);
if (verbose>1); disp(sprintf(' %% A_nlo__ and A_nhi__: %0.3fs',tmp_t)); end;
%%%%%%%%;
%}
A_nlp_emp__ = 2*(0.0 + min(A_nlo__,A_nhi__))/n_shuffle; A_nlp_emp__ = -log(A_nlp_emp__);

%{
%%%%;
figure(3); plot(A_nlp_opt__(:),A_nlp_emp__(:),'.',[0,10],[0,10],'k-','LineWidth',3);
%%%%;
figure(4);
subplot(2,1,1); hist(exp(-A_nlp_opt__(:)),linspace(0,1,1024)); title('opt');
subplot(2,1,2); hist(exp(-A_nlp_emp__(:)),linspace(0,1,1024)); title('emp');
%%%%;
n_h = 128;
h__ = zeros(n_retain,n_h);
for nretain=0:n_retain-1;
h_x_ = linspace(A_avg_(1+nretain) - 2.5*A_std_(1+nretain),A_avg_(1+nretain) + 2.5*A_std_(1+nretain),n_h);
h__(1+nretain,:) = hist(A_prm__(1+nretain,:),h_x_);
end;%for nretain=0:n_retain-1;
figure(2);imagesc(h__);figbeach();
%%%%;
%}

A_nlp__ = min(A_nlp_opt__,A_nlp_emp__);
A_bnlp__ = max(0,A_nlp__-log(n_bonferroni));
A_sign__ = sign(A_prm__ - repmat(A_avg_,[1,1+n_shuffle]));
A_snlp__ = A_sign__.*A_nlp__;
A_sbnlp__ = A_sign__.*A_bnlp__;
%%%%;
[A_nlps__,index_uns_from_srt__] = sort(A_nlp__,1,'descend');
index_uns_from_srt__ = index_uns_from_srt__ - 1;
%{
%%%%%%%%;
% slower than vectorized version above. ;
%%%%%%%%;
A_nlps__ = zeros(n_retain,1+n_shuffle);
index_uns_from_srt__ = zeros(n_retain,1+n_shuffle);
for nshuffle=0:1+n_shuffle-1;
[A_nlps__(:,1+nshuffle),index_uns_from_srt__(:,1+nshuffle)] = sort(A_nlp__(:,1+nshuffle),'descend');
index_uns_from_srt__ = index_uns_from_srt__ - 1;
end;%for nshuffle=0:1+n_shuffle-1;
%}

%%%%;
tmp_t = tic();
A_hnlps_emp__ = zeros(n_retain,1+n_shuffle);
A_hnlps_opt__ = zeros(n_retain,1+n_shuffle);
for nretain=0:n_retain-1;
tmp_A_ = A_nlps__(1+nretain,:);
[ ...
 ~ ...
,tmp_g_opt_ ...
,tmp_nlp_opt_ ...
,tmp_nlp_emp_ ...
,tmp_p_opt_ ...
,tmp_p_emp_ ...
] = ...
gumbel_fit_0( ...
 [] ...
,tmp_A_(2:end) ...
,tmp_A_(1) ...
);
A_hnlps_opt__(1+nretain,:) = tmp_nlp_opt_;
A_hnlps_emp__(1+nretain,:) = tmp_nlp_emp_;
end;%for nretain=0:n_retain-1;
tmp_t = toc(tmp_t);
if (verbose>1); disp(sprintf(' %% gumbel_fit_0: %0.3fs',tmp_t)); end;
%%%%;
A_nlpuh_emp__ = zeros(n_retain,1+n_shuffle);
A_nlpuh_opt__ = zeros(n_retain,1+n_shuffle);
for nshuffle=0:1+n_shuffle-1;
A_nlpuh_emp__(1+index_uns_from_srt__(:,1+nshuffle),1+nshuffle) = A_hnlps_emp__(:,1+nshuffle);
A_nlpuh_opt__(1+index_uns_from_srt__(:,1+nshuffle),1+nshuffle) = A_hnlps_opt__(:,1+nshuffle);
end;%for nshuffle=0:1+n_shuffle-1;
%%%%;
%{
%%%%%%%%;
% slower than vectorized version below. ;
%%%%%%%%;
A_hmnlps_emp__ = A_hnlps_emp__;
A_hmnlps_opt__ = A_hnlps_opt__;
for nshuffle=0:1+n_shuffle-1;
for nretain=1:n_retain-1;
A_hmnlps_emp__(1+nretain,1+nshuffle) = min(A_hmnlps_emp__(1+nretain,1+nshuffle),A_hmnlps_emp__(1+nretain-1,1+nshuffle));
A_hmnlps_opt__(1+nretain,1+nshuffle) = min(A_hmnlps_opt__(1+nretain,1+nshuffle),A_hmnlps_opt__(1+nretain-1,1+nshuffle));
end;%for nretain=1:n_retain-1;
end;%for nshuffle=0:1+n_shuffle-1;
toc,;
%}
%%%%%%%%;
% faster than loop above. ;
%%%%%%%%;
A_hmnlps_emp__ = cummin(A_hnlps_emp__,1);
A_hmnlps_opt__ = cummin(A_hnlps_opt__,1);
%%%%;
A_hmnlp_emp__ = zeros(n_retain,1+n_shuffle);
A_hmnlp_opt__ = zeros(n_retain,1+n_shuffle);
for nshuffle=0:1+n_shuffle-1;
A_hmnlp_emp__(1+index_uns_from_srt__(:,1+nshuffle),1+nshuffle) = A_hmnlps_emp__(:,1+nshuffle);
A_hmnlp_opt__(1+index_uns_from_srt__(:,1+nshuffle),1+nshuffle) = A_hmnlps_opt__(:,1+nshuffle);
end;%for nshuffle=0:1+n_shuffle-1;
%%%%%%%%;
A_hnlp__ = min(A_hmnlp_emp__,A_hmnlp_opt__);
A_shnlp__ = A_sign__.*A_hnlp__;
A_jnlp__ = min(A_hnlp__,A_nlp__);
A_sjnlp__ = A_sign__.*A_jnlp__;
%%%%%%%%;
A_sqr_nlp__ = zeros(n_var^2,1+n_shuffle);
A_sqr_nlp__(1+index_retain_,:) = A_nlp__;
A_sqr_nlp___ = reshape(A_sqr_nlp__,[n_var,n_var,1+n_shuffle]);
%%%%;
A_sqr_snlp__ = zeros(n_var^2,1+n_shuffle);
A_sqr_snlp__(1+index_retain_,:) = A_snlp__;
A_sqr_snlp___ = reshape(A_sqr_snlp__,[n_var,n_var,1+n_shuffle]);
%%%%;
A_sqr_bnlp__ = zeros(n_var^2,1+n_shuffle);
A_sqr_bnlp__(1+index_retain_,:) = A_bnlp__;
A_sqr_bnlp___ = reshape(A_sqr_bnlp__,[n_var,n_var,1+n_shuffle]);
%%%%;
A_sqr_sbnlp__ = zeros(n_var^2,1+n_shuffle);
A_sqr_sbnlp__(1+index_retain_,:) = A_sbnlp__;
A_sqr_sbnlp___ = reshape(A_sqr_sbnlp__,[n_var,n_var,1+n_shuffle]);
%%%%;
A_sqr_hnlp__ = zeros(n_var^2,1+n_shuffle);
A_sqr_hnlp__(1+index_retain_,:) = A_hnlp__;
A_sqr_hnlp___ = reshape(A_sqr_hnlp__,[n_var,n_var,1+n_shuffle]);
%%%%;
A_sqr_shnlp__ = zeros(n_var^2,1+n_shuffle);
A_sqr_shnlp__(1+index_retain_,:) = A_shnlp__;
A_sqr_shnlp___ = reshape(A_sqr_shnlp__,[n_var,n_var,1+n_shuffle]);
%%%%;
A_sqr_jnlp__ = zeros(n_var^2,1+n_shuffle);
A_sqr_jnlp__(1+index_retain_,:) = A_jnlp__;
A_sqr_jnlp___ = reshape(A_sqr_jnlp__,[n_var,n_var,1+n_shuffle]);
%%%%;
A_sqr_sjnlp__ = zeros(n_var^2,1+n_shuffle);
A_sqr_sjnlp__(1+index_retain_,:) = A_sjnlp__;
A_sqr_sjnlp___ = reshape(A_sqr_sjnlp__,[n_var,n_var,1+n_shuffle]);
%%%%;

%%%%%%%%;
% diagnostics. ;
%%%%%%%%;
flag_check=0;
if flag_check;
%%%%%%%%;
index_mask_ = efind(ones(n_var,n_var)-eye(n_var,n_var));
x_val_ = linspace(0,1,1024); n_x_val = numel(x_val_);
tmp_0nlp_best__ = zeros(n_x_val,1+n_shuffle);
tmp_bnlp_best__ = zeros(n_x_val,1+n_shuffle);
tmp_jnlp_best__ = zeros(n_x_val,1+n_shuffle);
for nshuffle=0:1+n_shuffle-1;
tmp_A_0nlp__ =  A_sqr_nlp___(:,:,1+nshuffle); tmp_A_0nlp_ = tmp_A_0nlp__(1+index_mask_);
tmp_A_bnlp__ = A_sqr_bnlp___(:,:,1+nshuffle); tmp_A_bnlp_ = tmp_A_bnlp__(1+index_mask_);
tmp_A_jnlp__ = A_sqr_jnlp___(:,:,1+nshuffle); tmp_A_jnlp_ = tmp_A_jnlp__(1+index_mask_);
for nx_val=0:n_x_val-1;
tmp_x_val = x_val_(1+nx_val);
tmp_nlx_val = -log(tmp_x_val);
tmp_0nlp_best__(1+nx_val,1+nshuffle) = sum(tmp_A_0nlp_>=tmp_nlx_val)>0;
tmp_bnlp_best__(1+nx_val,1+nshuffle) = sum(tmp_A_bnlp_>=tmp_nlx_val)>0;
tmp_jnlp_best__(1+nx_val,1+nshuffle) = sum(tmp_A_jnlp_>=tmp_nlx_val)>0;
end;%for nx_val=0:n_x_val-1;
end;%for nshuffle=0:1+n_shuffle-1;
figure(1);clf;figsml;
linewidth_use = 3;
markersize_use = 16;
hold on;
plot(x_val_,x_val_,'-','Color',0.65*[1,1,1],'LineWidth',2);
%stairs(x_val_,mean(tmp_0nlp_best__(:,2:end),2),'b.-','LineWidth',linewidth_use,'MarkerSize',markersize_use); %<-- should be straight line. ;
stairs(x_val_,mean(tmp_bnlp_best__(:,2:end),2),'c.-','LineWidth',linewidth_use,'MarkerSize',markersize_use); %<-- should be straight line. ;
stairs(x_val_,mean(tmp_jnlp_best__(:,2:end),2),'r.-','LineWidth',linewidth_use,'MarkerSize',markersize_use); %<-- should be straight line. ;
hold off;
legend({'$P(x)=x$','$p_{b}$','$p_{h}$'},'Location','NorthWest','Interpreter','latex');
xlim([0,1]); ylim([0,1]);
set(gca,'XTick',0:0.1:1.0);
set(gca,'YTick',0:0.1:1.0);
ylabel('cumulative-distribution $P(x)$','Interpreter','latex');
xlabel('fraction $x$','Interpreter','latex');
grid on;
set(gca,'FontSize',12);
%%%%%%%%;
end;%if flag_check;

if (verbose); disp(sprintf(' %% [finished dolphin_holm_from_prm_0]')); end;
