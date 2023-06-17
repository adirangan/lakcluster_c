function ...
[ ...
 parameter ...
,A_dat_nlp__ ...
,A_dat_bnlp__ ...
,A_dat_hnlp__ ...
,A_dat_jnlp__ ...
,A_dat_snlp__ ...
,A_dat_sbnlp__ ...
,A_dat_shnlp__ ...
,A_dat_sjnlp__ ...
] = ...
holm_from_prm_0( ...
 parameter ...
,n_var ...
,mask_ ...
,n_shuffle ...
,A_dat_prm__ ...
);

verbose=0;
if (verbose); disp(sprintf(' %% [entering holm_from_prm_0]')); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_var=[]; end; na=na+1;
if (nargin<1+na); mask_=[]; end; na=na+1;
if (nargin<1+na); n_shuffle=[]; end; na=na+1;
if (nargin<1+na); A_dat_prm__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_onesided_vs_twosided'); parameter.flag_onesided_vs_twosided=0; end;
flag_onesided_vs_twosided = parameter.flag_onesided_vs_twosided;
if isempty(mask_); mask_ = ones(n_var,1); end;
n_bonferroni = sum(mask_);

assert(size(mask_,1)==n_var);
assert(size(A_dat_prm__,1)==n_var);
assert(size(A_dat_prm__,2)==1+n_shuffle);

%%%%%%%%;
index_retain_ = efind(mask_);
n_retain = numel(index_retain_);
A_prm__ = A_dat_prm__(1+index_retain_,:);
A_avg_ = mean(A_prm__(:,2:end),2); A_std_ = std(A_prm__(:,2:end),1,2);
A_zsc__ = (A_prm__ - A_avg_)./A_std_;
if flag_onesided_vs_twosided==0;
A_nlp_opt__ = -z_to_p_twosided_0(A_zsc__);
end;%if flag_onesided_vs_twosided==0;
if flag_onesided_vs_twosided==1;
A_nlp_opt__ = -z_to_p_0(A_zsc__); %<-- onesided test. ;
end;%if flag_onesided_vs_twosided==1;
%%%%%%%%;

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
if flag_onesided_vs_twosided==0;
A_nlp_emp__ = 2*(0.0 + min(A_nlo__,A_nhi__))/n_shuffle; A_nlp_emp__ = -log(A_nlp_emp__);
end;%if flag_onesided_vs_twosided==0;
if flag_onesided_vs_twosided==1;
A_nlp_emp__ = 1*(0.0 + A_nhi__)/n_shuffle; A_nlp_emp__ = -log(A_nlp_emp__);
end;%if flag_onesided_vs_twosided==1;
%%%%%%%%;
A_nlp__ = min(A_nlp_opt__,A_nlp_emp__);
A_bnlp__ = max(0,A_nlp__-log(n_bonferroni));
A_sign__ = sign(A_prm__ - repmat(A_avg_,[1,1+n_shuffle]));
A_snlp__ = A_sign__.*A_nlp__;
A_sbnlp__ = A_sign__.*A_bnlp__;
%%%%;
[A_nlps__,ij_uns_from_srt__] = sort(A_nlp__,1,'descend');
index_uns_from_srt__ = ij_uns_from_srt__ - 1;

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
A_dat_nlp__ = zeros(n_var,1+n_shuffle);
A_dat_nlp__(1+index_retain_,:) = A_nlp__;
%%%%;
A_dat_snlp__ = zeros(n_var,1+n_shuffle);
A_dat_snlp__(1+index_retain_,:) = A_snlp__;
%%%%;
A_dat_bnlp__ = zeros(n_var,1+n_shuffle);
A_dat_bnlp__(1+index_retain_,:) = A_bnlp__;
%%%%;
A_dat_sbnlp__ = zeros(n_var,1+n_shuffle);
A_dat_sbnlp__(1+index_retain_,:) = A_sbnlp__;
%%%%;
A_dat_hnlp__ = zeros(n_var,1+n_shuffle);
A_dat_hnlp__(1+index_retain_,:) = A_hnlp__;
%%%%;
A_dat_shnlp__ = zeros(n_var,1+n_shuffle);
A_dat_shnlp__(1+index_retain_,:) = A_shnlp__;
%%%%;
A_dat_jnlp__ = zeros(n_var,1+n_shuffle);
A_dat_jnlp__(1+index_retain_,:) = A_jnlp__;
%%%%;
A_dat_sjnlp__ = zeros(n_var,1+n_shuffle);
A_dat_sjnlp__(1+index_retain_,:) = A_sjnlp__;
%%%%;

if (verbose); disp(sprintf(' %% [finished holm_from_prm_0]')); end;
