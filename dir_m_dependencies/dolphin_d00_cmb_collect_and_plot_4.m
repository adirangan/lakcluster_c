function ...
[ ...
 parameter ...
 ,n_found_zzzz ...
 ,a_zzzz_o_ ...
 ,A_zzzz_o__ ...
 ,BB_zzzz_o__ ...
 ,CC_zzzz_o__ ...
 ,a_sbnlp_zzzz_ ...
 ,A_sbnlp_zzzz__ ...
 ,BB_sbnlp_zzzz__ ...
 ,CC_sbnlp_zzzz__ ...
 ,a_snlp_zzzz_ ...
 ,A_snlp_zzzz__ ...
 ,BB_snlp_zzzz__ ...
 ,CC_snlp_zzzz__ ...
 ,a_nlp_zzzz_ ...
 ,A_nlp_zzzz__ ...
 ,BB_nlp_zzzz__ ...
 ,CC_nlp_zzzz__ ...
 ,a_avg_zzzz_ ...
 ,A_avg_zzzz__ ...
 ,BB_avg_zzzz__ ...
 ,CC_avg_zzzz__ ...
 ,a_std_zzzz_ ...
 ,A_std_zzzz__ ...
 ,BB_std_zzzz__ ...
 ,CC_std_zzzz__ ...
] = ...
dolphin_d00_cmb_collect_and_plot_4( ...
 parameter ...
 ,string_dat_name_ ...
 ,dir_mat ...
 ,dir_shuffle_mat ...
 ,dir_jpg ...
 ,fname_infix_aaaa ...
 ,fname_infix_bbbb ...
 ,fname_infix_zzzz ...
);

verbose=1;
if (verbose); disp(sprintf(' %% [entering dolphin_d00_cmb_collect_and_plot_4]')); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); string_dat_name_=[]; end; na=na+1;
if (nargin<1+na); dir_mat=[]; end; na=na+1;
if (nargin<1+na); dir_shuffle_mat=[]; end; na=na+1;
if (nargin<1+na); dir_jpg=[]; end; na=na+1;
if (nargin<1+na); fname_infix_aaaa=[]; end; na=na+1;
if (nargin<1+na); fname_infix_bbbb=[]; end; na=na+1;
if (nargin<1+na); fname_infix_zzzz=[]; end; na=na+1;

% try: ;
% dolphin_dXX_load_7;
% parameter = []; string_dat_name_ = string_dat_ori_name_; dir_mat = '/home/rangan/dir_bcc/dir_dolphin/dir_mat_ind_7'; dir_shuffle_mat = '/home/rangan/dir_bcc/dir_dolphin/dir_shuffle256_mat_ind_7'; dir_jpg = '/home/rangan/dir_bcc/dir_dolphin/dir_jpg_ind_7'; fname_infix_aaaa = 'aid00_cmb'; fname_infix_bbbb = ''; fname_infix_zzzz = ''; parameter = struct('type','parameter'); parameter.index_var_retain_=transpose(0:42);
% parameter = []; string_dat_name_ = string_dat_ori_name_; dir_mat = '/home/rangan/dir_bcc/dir_dolphin/dir_mat_ind_7'; dir_shuffle_mat = '/home/rangan/dir_bcc/dir_dolphin/dir_shuffle256_mat_ind_7'; dir_jpg = '/home/rangan/dir_bcc/dir_dolphin/dir_jpg_ind_7'; fname_infix_aaaa = 'aid00_M_cmb'; fname_infix_bbbb = 'aid00_F_cmb'; fname_infix_zzzz = 'aid00_M_vs_F_cmb'; fname_infix_aaaa = 'aid00_age_7099_cmb'; fname_infix_bbbb = 'aid00_age_2070_cmb'; fname_infix_zzzz = 'aid00_age_7099_vs_2070_cmb'; parameter = struct('type','parameter'); parameter.index_var_retain_=transpose(0:42);

if (isempty(parameter)); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'flag_replot')); parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'index_var_retain_')); parameter.index_var_retain_ = []; end; %<-- parameter_bookmark. ;
flag_replot = parameter.flag_replot;
index_var_retain_ = parameter.index_var_retain_;

if (  isempty(fname_infix_bbbb) |  isempty(fname_infix_zzzz) );
flag_single = 1;
flag_double = 0;
fname_infix_zzzz = fname_infix_aaaa;
end;%if (  isempty(fname_infix_bbbb) |  isempty(fname_infix_zzzz) );

if ( ~isempty(fname_infix_bbbb) & ~isempty(fname_infix_zzzz) );
flag_single = 0;
flag_double = 1;
end;%if ( ~isempty(fname_infix_bbbb) | ~isempty(fname_infix_zzzz) );

fname_mat_aaaa_o = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix_aaaa);
fname_mat_aaaa_s = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix_aaaa);
flag_exist_aaaa = exist(fname_mat_aaaa_o,'file') & exist(fname_mat_aaaa_s,'file');
if ~flag_exist_aaaa;
disp(sprintf(' %% %s %s not found, skipping',fname_mat_aaaa_o,fname_mat_aaaa_s));
end;%if ~flag_exist_aaaa;

n_var = numel(string_dat_name_);
if (isempty(index_var_retain_)); index_var_retain_ = transpose(0:n_var-1); end;
n_var_retain = numel(index_var_retain_);
string_dat_name_retain_ = string_dat_name_(1+index_var_retain_);

flag_exist = flag_exist_aaaa;

if flag_double;
fname_mat_bbbb_o = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix_bbbb);
fname_mat_bbbb_s = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix_bbbb);
flag_exist_bbbb = exist(fname_mat_bbbb_o,'file') & exist(fname_mat_bbbb_s,'file');
if ~flag_exist_bbbb;
disp(sprintf(' %% %s %s not found, skipping',fname_mat_bbbb_o,fname_mat_bbbb_s));
end;%if ~flag_exist_bbbb;
flag_exist = flag_exist_aaaa & flag_exist_bbbb;
end;%if flag_double;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_exist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
if flag_single;
%%%%%%%%;
fname_mat_zzzz_o = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix_zzzz);
if ( exist(fname_mat_zzzz_o,'file'));
tmp_zzzz_o_ = load(fname_mat_zzzz_o);
a_zzzz_o_ = tmp_zzzz_o_.a_cmb_; n_var = numel(a_zzzz_o_);
A_zzzz_o__ = tmp_zzzz_o_.A_cmb__;
BB_inv_zzzz_o__ = tmp_zzzz_o_.BB_inv_cmb__;
CC_inv_zzzz_o__ = tmp_zzzz_o_.CC_inv_cmb__;
[ ...
 parameter ...
,a_duo_zzzz_o__ ...
,A_duo_zzzz_o___ ...
,BB_inv_duo_zzzz_o___ ...
,CC_inv_duo_zzzz_o___ ...
,index_nvar0_from_nvv_ ...
,index_nvar1_from_nvv_ ...
] = ...
dolphin_duo_from_each_0( ...
 parameter ...
,tmp_zzzz_o_.A_cmb__ ...
,tmp_zzzz_o_.BB_inv_cmb__ ...
,tmp_zzzz_o_.CC_inv_cmb__ ...
,tmp_zzzz_o_.a_each__ ...
,tmp_zzzz_o_.A_each__ ...
,tmp_zzzz_o_.BB_inv_each__ ...
,tmp_zzzz_o_.CC_inv_each__ ...
);
end;%if ( exist(fname_mat_zzzz_s,'file'));
%%%%%%%%;
a_avg_zzzz_ = zeros(n_var,1); a_std_zzzz_ = zeros(n_var,1);
A_avg_zzzz__ = zeros(n_var,n_var); A_std_zzzz__ = zeros(n_var,n_var);
BB_inv_avg_zzzz__ = zeros(n_var,n_var); BB_inv_std_zzzz__ = zeros(n_var,n_var);
CC_inv_avg_zzzz__ = zeros(n_var,n_var); CC_inv_std_zzzz__ = zeros(n_var,n_var);
n_found_zzzz = 0;
fname_mat_zzzz_s = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix_zzzz);
if ( exist(fname_mat_zzzz_s,'file'));
tmp_zzzz_s_ = load(fname_mat_zzzz_s);
n_found_zzzz = size(tmp_zzzz_s_.L_cmb_prm_,1) - 1;
n_shuffle = n_found_zzzz; n_shuffle_use = n_shuffle;
%%%%;
[ ...
 parameter ...
,a_duo_zzzz_s___ ...
,A_duo_zzzz_s____ ...
,BB_inv_duo_zzzz_s____ ...
,CC_inv_duo_zzzz_s____ ...
,index_nvar0_from_nvv_ ...
,index_nvar1_from_nvv_ ...
] = ...
dolphin_duo_from_each_prm_0( ...
 parameter ...
,tmp_zzzz_s_.A_cmb_prm___ ...
,tmp_zzzz_s_.BB_inv_cmb_prm___ ...
,tmp_zzzz_s_.CC_inv_cmb_prm___ ...
,tmp_zzzz_s_.a_each_prm___ ...
,tmp_zzzz_s_.A_each_prm___ ...
,tmp_zzzz_s_.BB_inv_each_prm___ ...
,tmp_zzzz_s_.CC_inv_each_prm___ ...
);
%%%%;
a_avg_zzzz_ = mean(tmp_zzzz_s_.a_cmb_prm__(:,2:end),2); a_std_zzzz_ = std(tmp_zzzz_s_.a_cmb_prm__(:,2:end),1,2);
A_avg_zzzz__ = mean(tmp_zzzz_s_.A_cmb_prm___(:,:,2:end),3); A_std_zzzz__ = std(tmp_zzzz_s_.A_cmb_prm___(:,:,2:end),1,3);
BB_inv_avg_zzzz__ = mean(tmp_zzzz_s_.BB_inv_cmb_prm___(:,:,2:end),3); BB_inv_std_zzzz__ = std(tmp_zzzz_s_.BB_inv_cmb_prm___(:,:,2:end),1,3);
CC_inv_avg_zzzz__ = mean(tmp_zzzz_s_.CC_inv_cmb_prm___(:,:,2:end),3); CC_inv_std_zzzz__ = std(tmp_zzzz_s_.CC_inv_cmb_prm___(:,:,2:end),1,3);
end;%if ( exist(fname_mat_zzzz_s,'file'));
a_nlp_zzzz_ = -z_to_p_twosided_0((a_zzzz_o_ - a_avg_zzzz_)./a_std_zzzz_);
A_nlp_zzzz__ = -z_to_p_twosided_0((A_zzzz_o__ - A_avg_zzzz__)./A_std_zzzz__);
BB_inv_nlp_zzzz__ = -z_to_p_twosided_0((BB_inv_zzzz_o__ - BB_inv_avg_zzzz__)./BB_inv_std_zzzz__);
CC_inv_nlp_zzzz__ = -z_to_p_twosided_0((CC_inv_zzzz_o__ - CC_inv_avg_zzzz__)./CC_inv_std_zzzz__);
a_snlp_zzzz_ = sign(a_zzzz_o_ - a_avg_zzzz_).*a_nlp_zzzz_;
A_snlp_zzzz__ = sign(A_zzzz_o__ - A_avg_zzzz__).*A_nlp_zzzz__;
BB_inv_snlp_zzzz__ = sign(BB_inv_zzzz_o__ - BB_inv_avg_zzzz__).*BB_inv_nlp_zzzz__;
CC_inv_snlp_zzzz__ = sign(CC_inv_zzzz_o__ - CC_inv_avg_zzzz__).*CC_inv_nlp_zzzz__;
a_sbnlp_zzzz_ = sign(a_zzzz_o_ - a_avg_zzzz_).*max(0,a_nlp_zzzz_-log(n_var_retain));
A_sbnlp_zzzz__ = sign(A_zzzz_o__ - A_avg_zzzz__).*max(0,A_nlp_zzzz__-log(n_var_retain*(n_var_retain-1)));
BB_inv_sbnlp_zzzz__ = sign(BB_inv_zzzz_o__ - BB_inv_avg_zzzz__).*max(0,BB_inv_nlp_zzzz__-log(n_var_retain*(n_var_retain-1)/2));
CC_inv_sbnlp_zzzz__ = sign(CC_inv_zzzz_o__ - CC_inv_avg_zzzz__).*max(0,CC_inv_nlp_zzzz__-log(n_var_retain*(n_var_retain-1)/2));
%%%%%%%%;
a_avg_zzzz_ = a_avg_zzzz_(1+index_var_retain_);
a_std_zzzz_ = a_std_zzzz_(1+index_var_retain_);
a_nlp_zzzz_ = a_nlp_zzzz_(1+index_var_retain_);
a_snlp_zzzz_ = a_snlp_zzzz_(1+index_var_retain_);
a_sbnlp_zzzz_ = a_sbnlp_zzzz_(1+index_var_retain_);
A_avg_zzzz__ = A_avg_zzzz__(1+index_var_retain_,1+index_var_retain_);
A_std_zzzz__ = A_std_zzzz__(1+index_var_retain_,1+index_var_retain_);
A_nlp_zzzz__ = A_nlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
A_snlp_zzzz__ = A_snlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
A_sbnlp_zzzz__ = A_sbnlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_avg_zzzz__ = BB_inv_avg_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_std_zzzz__ = BB_inv_std_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_nlp_zzzz__ = BB_inv_nlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_snlp_zzzz__ = BB_inv_snlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_sbnlp_zzzz__ = BB_inv_sbnlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_avg_zzzz__ = CC_inv_avg_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_std_zzzz__ = CC_inv_std_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_nlp_zzzz__ = CC_inv_nlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_snlp_zzzz__ = CC_inv_snlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_sbnlp_zzzz__ = CC_inv_sbnlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
%%%%%%%%;
end;%if flag_single;
%%%%%%%%

%%%%%%%%;
if flag_double;
%%%%%%%%;
fname_mat_aaaa_o = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix_aaaa);
if ( exist(fname_mat_aaaa_o,'file'));
tmp_aaaa_o_ = load(fname_mat_aaaa_o);
a_aaaa_o_ = tmp_aaaa_o_.a_cmb_; n_var = numel(a_aaaa_o_);
A_aaaa_o__ = tmp_aaaa_o_.A_cmb__;
BB_inv_aaaa_o__ = tmp_aaaa_o_.BB_inv_cmb__;
CC_inv_aaaa_o__ = tmp_aaaa_o_.CC_inv_cmb__;
end;%if ( exist(fname_mat_aaaa_s,'file'));
fname_mat_bbbb_o = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix_bbbb);
if ( exist(fname_mat_bbbb_o,'file'));
tmp_bbbb_o_ = load(fname_mat_bbbb_o);
a_bbbb_o_ = tmp_bbbb_o_.a_cmb_;
A_bbbb_o__ = tmp_bbbb_o_.A_cmb__;
BB_inv_bbbb_o__ = tmp_bbbb_o_.BB_inv_cmb__;
CC_inv_bbbb_o__ = tmp_bbbb_o_.CC_inv_cmb__;
end;%if ( exist(fname_mat_bbbb_s,'file'));
a_zzzz_o_ = a_aaaa_o_ - a_bbbb_o_;
A_zzzz_o__ = A_aaaa_o__ - A_bbbb_o__;
BB_inv_zzzz_o__ = BB_inv_aaaa_o__ - BB_inv_bbbb_o__;
CC_inv_zzzz_o__ = CC_inv_aaaa_o__ - CC_inv_bbbb_o__;
%%%%%%%%;
a_avg_zzzz_ = zeros(n_var,1); a_std_zzzz_ = zeros(n_var,1);
A_avg_zzzz__ = zeros(n_var,n_var); A_std_zzzz__ = zeros(n_var,n_var);
BB_inv_avg_zzzz__ = zeros(n_var,n_var); BB_inv_std_zzzz__ = zeros(n_var,n_var);
CC_inv_avg_zzzz__ = zeros(n_var,n_var); CC_inv_std_zzzz__ = zeros(n_var,n_var);
%%%%%%%%;
fname_mat_aaaa_s = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix_aaaa);
fname_mat_bbbb_s = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix_bbbb);
if ( exist(fname_mat_aaaa_s,'file') &  exist(fname_mat_bbbb_s,'file'));
tmp_aaaa_s_ = load(fname_mat_aaaa_s);
tmp_bbbb_s_ = load(fname_mat_bbbb_s);
n_found_zzzz = min(numel(tmp_aaaa_s_.L_cmb_prm_),numel(tmp_bbbb_s_.L_cmb_prm_)) - 1;
n_shuffle = n_found_zzzz; n_shuffle_use = n_shuffle^2;
tmp_ = dolphin_doublediff_0(tmp_aaaa_s_.a_cmb_prm__,tmp_bbbb_s_.a_cmb_prm__);
a_avg_zzzz_ = mean(tmp_(:,2:end),2);
a_std_zzzz_ = std(tmp_(:,2:end),1,2);
tmp_ = dolphin_doublediff_0(tmp_aaaa_s_.A_cmb_prm___,tmp_bbbb_s_.A_cmb_prm___);
A_avg_zzzz__ = mean(tmp_(:,:,2:end),3);
A_std_zzzz__ = std(tmp_(:,:,2:end),1,3);
tmp_ = dolphin_doublediff_0(tmp_aaaa_s_.BB_inv_cmb_prm___,tmp_bbbb_s_.BB_inv_cmb_prm___);
BB_inv_avg_zzzz__ = mean(tmp_(:,:,2:end),3);
BB_inv_std_zzzz__ = std(tmp_(:,:,2:end),1,3);
tmp_ = dolphin_doublediff_0(tmp_aaaa_s_.BB_inv_cmb_prm___,tmp_bbbb_s_.BB_inv_cmb_prm___);
CC_inv_avg_zzzz__ = mean(tmp_(:,:,2:end),3);
CC_inv_std_zzzz__ = std(tmp_(:,:,2:end),1,3);
clear tmp_;
end;%if ( exist(fname_mat_aaaa_s,'file') &  exist(fname_mat_bbbb_s,'file'));
%%%%;
tmp_zzzz_s_ = struct('type','tmp');
tmp_zzzz_s_.a_cmb_prm__ = dolphin_doublediff_0(tmp_aaaa_s_.a_cmb_prm__,tmp_bbbb_s_.a_cmb_prm__);
tmp_zzzz_s_.A_cmb_prm___ = dolphin_doublediff_0(tmp_aaaa_s_.A_cmb_prm___,tmp_bbbb_s_.A_cmb_prm___);
tmp_zzzz_s_.BB_inv_cmb_prm___ = dolphin_doublediff_0(tmp_aaaa_s_.BB_inv_cmb_prm___,tmp_bbbb_s_.BB_inv_cmb_prm___);
tmp_zzzz_s_.CC_inv_cmb_prm___ = dolphin_doublediff_0(tmp_aaaa_s_.CC_inv_cmb_prm___,tmp_bbbb_s_.CC_inv_cmb_prm___);
tmp_zzzz_s_.L_cmb_prm_ = dolphin_doublediff_0(tmp_aaaa_s_.L_cmb_prm_,tmp_bbbb_s_.L_cmb_prm_);
%%%%;
a_nlp_zzzz_ = -z_to_p_twosided_0((a_zzzz_o_ - a_avg_zzzz_)./a_std_zzzz_);
A_nlp_zzzz__ = -z_to_p_twosided_0((A_zzzz_o__ - A_avg_zzzz__)./A_std_zzzz__);
BB_inv_nlp_zzzz__ = -z_to_p_twosided_0((BB_inv_zzzz_o__ - BB_inv_avg_zzzz__)./BB_inv_std_zzzz__);
CC_inv_nlp_zzzz__ = -z_to_p_twosided_0((CC_inv_zzzz_o__ - CC_inv_avg_zzzz__)./CC_inv_std_zzzz__);
a_snlp_zzzz_ = sign(a_zzzz_o_ - a_avg_zzzz_).*a_nlp_zzzz_;
A_snlp_zzzz__ = sign(A_zzzz_o__ - A_avg_zzzz__).*A_nlp_zzzz__;
BB_inv_snlp_zzzz__ = sign(BB_inv_zzzz_o__ - BB_inv_avg_zzzz__).*BB_inv_nlp_zzzz__;
CC_inv_snlp_zzzz__ = sign(CC_inv_zzzz_o__ - CC_inv_avg_zzzz__).*CC_inv_nlp_zzzz__;
a_sbnlp_zzzz_ = sign(a_zzzz_o_ - a_avg_zzzz_).*max(0,a_nlp_zzzz_-log(n_var_retain));
A_sbnlp_zzzz__ = sign(A_zzzz_o__ - A_avg_zzzz__).*max(0,A_nlp_zzzz__-log(n_var_retain*(n_var_retain-1)));
BB_inv_sbnlp_zzzz__ = sign(BB_inv_zzzz_o__ - BB_inv_avg_zzzz__).*max(0,BB_inv_nlp_zzzz__-log(n_var_retain*(n_var_retain-1)/2));
CC_inv_sbnlp_zzzz__ = sign(CC_inv_zzzz_o__ - CC_inv_avg_zzzz__).*max(0,CC_inv_nlp_zzzz__-log(n_var_retain*(n_var_retain-1)/2));
%%%%%%%%;
a_avg_zzzz_ = a_avg_zzzz_(1+index_var_retain_);
a_std_zzzz_ = a_std_zzzz_(1+index_var_retain_);
a_nlp_zzzz_ = a_nlp_zzzz_(1+index_var_retain_);
a_snlp_zzzz_ = a_snlp_zzzz_(1+index_var_retain_);
a_sbnlp_zzzz_ = a_sbnlp_zzzz_(1+index_var_retain_);
A_avg_zzzz__ = A_avg_zzzz__(1+index_var_retain_,1+index_var_retain_);
A_std_zzzz__ = A_std_zzzz__(1+index_var_retain_,1+index_var_retain_);
A_nlp_zzzz__ = A_nlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
A_snlp_zzzz__ = A_snlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
A_sbnlp_zzzz__ = A_sbnlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_avg_zzzz__ = BB_inv_avg_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_std_zzzz__ = BB_inv_std_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_nlp_zzzz__ = BB_inv_nlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_snlp_zzzz__ = BB_inv_snlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
BB_inv_sbnlp_zzzz__ = BB_inv_sbnlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_avg_zzzz__ = CC_inv_avg_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_std_zzzz__ = CC_inv_std_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_nlp_zzzz__ = CC_inv_nlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_snlp_zzzz__ = CC_inv_snlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
CC_inv_sbnlp_zzzz__ = CC_inv_sbnlp_zzzz__(1+index_var_retain_,1+index_var_retain_);
%%%%%%%%;
end;%if flag_double;
%%%%%%%%;

%%%%%%%%;
disp(sprintf(' %% n_found_zzzz %d',n_found_zzzz));
%%%%%%%%;

%%%%%%%%;
% Test out holm-bonferroni test. ;
%%%%%%%%;
if (~exist('A_sqr_nlp___','var'));
[ ...
 ~ ...
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
 [] ...
,n_var_retain ...
,'A' ...
,n_shuffle_use ...
,tmp_zzzz_s_.A_cmb_prm___(1+index_var_retain_,1+index_var_retain_,:) ...
);
end;%if (~exist('A_sqr_nlp___','var'));
A_nlp_zzzz__ = A_sqr_nlp___(:,:,1+0);
A_snlp_zzzz__ = A_sqr_snlp___(:,:,1+0);
A_sbnlp_zzzz__ = A_sqr_sbnlp___(:,:,1+0);
A_sjnlp_zzzz__ = A_sqr_sjnlp___(:,:,1+0);
%%%%%%%%;
if (~exist('BB_inv_sqr_nlp___','var'));
[ ...
 ~ ...
,BB_inv_sqr_nlp___ ...
,BB_inv_sqr_bnlp___ ...
,BB_inv_sqr_hnlp___ ...
,BB_inv_sqr_jnlp___ ...
,BB_inv_sqr_snlp___ ...
,BB_inv_sqr_sbnlp___ ...
,BB_inv_sqr_shnlp___ ...
,BB_inv_sqr_sjnlp___ ...
] = ...
dolphin_holm_from_prm_0( ...
 [] ...
,n_var_retain ...
,'B' ...
,n_shuffle_use ...
,tmp_zzzz_s_.BB_inv_cmb_prm___(1+index_var_retain_,1+index_var_retain_,:) ...
);
end;%if (~exist('BB_inv_sqr_nlp___','var'));
BB_inv_nlp_zzzz__ = BB_inv_sqr_nlp___(:,:,1+0);
BB_inv_snlp_zzzz__ = BB_inv_sqr_snlp___(:,:,1+0);
BB_inv_sbnlp_zzzz__ = BB_inv_sqr_sbnlp___(:,:,1+0);
BB_inv_sjnlp_zzzz__ = BB_inv_sqr_sjnlp___(:,:,1+0);
%%%%%%%%;
%{
if (~exist('CC_inv_sqr_nlp___','var'));
[ ...
 ~ ...
,CC_inv_sqr_nlp___ ...
,CC_inv_sqr_bnlp___ ...
,CC_inv_sqr_hnlp___ ...
,CC_inv_sqr_jnlp___ ...
,CC_inv_sqr_snlp___ ...
,CC_inv_sqr_sbnlp___ ...
,CC_inv_sqr_shnlp___ ...
,CC_inv_sqr_sjnlp___ ...
] = ...
dolphin_holm_from_prm_0( ...
 [] ...
,n_var_retain ...
,'B' ...
,n_shuffle_use ...
,tmp_zzzz_s_.CC_inv_cmb_prm___(1+index_var_retain_,1+index_var_retain_,:) ...
);
end;%if (~exist('CC_inv_sqr_nlp___','var'));
CC_inv_nlp_zzzz__ = CC_inv_sqr_nlp___(:,:,1+0);
CC_inv_snlp_zzzz__ = CC_inv_sqr_snlp___(:,:,1+0);
CC_inv_sbnlp_zzzz__ = CC_inv_sqr_sbnlp___(:,:,1+0);
CC_inv_sjnlp_zzzz__ = CC_inv_sqr_sjnlp___(:,:,1+0);
%%%%%%%%;
%}

for ntype=0:2-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ntype==0; tmp_str = 'A'; X_sqr_bnlp___ = A_sqr_bnlp___; X_sqr_jnlp___ = A_sqr_jnlp___; tmp_str_title = '$A$'; end;
if ntype==1; tmp_str = 'B'; X_sqr_bnlp___ = BB_inv_sqr_bnlp___; X_sqr_jnlp___ = BB_inv_sqr_jnlp___; tmp_str_title = '$BB^{T}$'; end;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_shuffle_%s_sjnlp_Pcumsum_FIGP',dir_jpg,fname_infix_zzzz,tmp_str);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
%%%%%%%%; 
index_mask_ = efind(ones(n_var_retain,n_var_retain)-eye(n_var_retain,n_var_retain)); %<-- do not need to use, since each nshuffle will have exactly n_var_retain elements on the diagonal, all of which are 0. ;
%%%%%%%%;
figure(1);clf;figsml;
linewidth_use = 3;
markersize_use = 16;
hold on;
%%%%;
%line([0,1],[0,1],'-','Color',0.65*[1,1,1],'LineWidth',2);
%%%%;
tmp__ = reshape(X_sqr_bnlp___,[n_var_retain^2,1+n_shuffle_use]);
tmp__ = sort(tmp__,1,'descend');
tmp__ = tmp__(1:n_var_retain*(n_var_retain-1),2:end); %<-- ignore the diagonals, as well as the original data. ;
tmp_ = max(tmp__,[],1);
stairs(sort(exp(-[Inf,tmp_,0])),[0,1:n_shuffle_use,n_shuffle_use]/n_shuffle_use,'c.-','LineWidth',linewidth_use,'MarkerSize',markersize_use);
clear tmp__ tmp_;
%%%%;
tmp__ = reshape(X_sqr_jnlp___,[n_var_retain^2,1+n_shuffle_use]);
tmp__ = sort(tmp__,1,'descend');
tmp__ = tmp__(1:n_var_retain*(n_var_retain-1),2:end); %<-- ignore the diagonals, as well as the original data. ;
tmp_ = max(tmp__,[],1);
stairs(sort(exp(-[Inf,tmp_,0])),[0,1:n_shuffle_use,n_shuffle_use]/n_shuffle_use,'r.-','LineWidth',linewidth_use,'MarkerSize',markersize_use);
clear tmp__ tmp_;
%%%%;
hold off;
legend({'$p_{b}$','$p_{h}$'},'Location','NorthWest','Interpreter','latex');
xlim([0,1]); ylim([0,1]);
set(gca,'XTick',0:0.1:1.0);
set(gca,'YTick',0:0.1:1.0);
ylabel('cumulative-distribution $P(x)$','Interpreter','latex');
xlabel('fraction $x$','Interpreter','latex');
grid on;
sgtitle(sprintf('%s',fname_fig),'Interpreter','none');
set(gca,'FontSize',12);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
%%%%%%%%; 
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for ntype=0:2-1;

for ntype=0:2-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ntype==0; tmp_str = 'A'; X_sbnlp_zzzz__ = A_sbnlp_zzzz__; X_sjnlp_zzzz__ = A_sjnlp_zzzz__; tmp_str_title = '$A$'; end;
if ntype==1; tmp_str = 'B'; X_sbnlp_zzzz__ = BB_inv_sbnlp_zzzz__; X_sjnlp_zzzz__ = BB_inv_sjnlp_zzzz__; tmp_str_title = '$BB^{T}$'; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_shuffle_%s_sxnlp_FIG%s',dir_jpg,fname_infix_zzzz,tmp_str,tmp_str);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
%%%%%%%%;
figure(1);clf;figbig;fontsize_use = 12; np=0;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1); c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1); c_nlpvt2_lim_ = [-27,+27];
%%%%;
axis image; set(gca,'ydir','reverse');
l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold on;
l_x_ = +0.5 + (n_var_retain)*[zeros(1,n_var_retain);ones(1,n_var_retain)];
l_y_ = +0.5 + repmat([0:n_var_retain-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
%%%%;
n_side = 4; r_rad = 1/sqrt(2);
tmp_w_ = linspace(0,2*pi,1+n_side); tmp_w_ = tmp_w_ + 0.5*mean(diff(tmp_w_)); n_w = numel(tmp_w_);
tmp_rcw_ = r_rad*cos(tmp_w_);
tmp_rsw_ = r_rad*sin(tmp_w_);
tmp_x__ = zeros(n_w,n_var_retain*(n_var_retain-1));
tmp_y__ = zeros(n_w,n_var_retain*(n_var_retain-1));
tmp_c___ = zeros(1,n_var_retain*(n_var_retain-1),3);
na=0;
for nvar0=0:n_var_retain-1;
for nvar1=0:n_var_retain-1;
if (nvar0~=nvar1);
tmp_x_ = 1.0 + nvar1 + tmp_rsw_;
tmp_y_ = 1.0 + nvar0 + tmp_rcw_;
tmp_val = X_sbnlp_zzzz__(1+nvar0,1+nvar1);
if (abs(tmp_val)>=-log(0.05));
tmp_val = (tmp_val - min(c_nlpvt2_lim_))/diff(c_nlpvt2_lim_);
nc_nlpvt2 = max(0,min(n_c_nlpvt2-1,floor(n_c_nlpvt2*tmp_val)));
tmp_c_ = c_nlpvt2__(1+nc_nlpvt2,:);
tmp_x__(:,1+na) = tmp_x_;
tmp_y__(:,1+na) = tmp_y_;
tmp_c___(1,1+na,:) = tmp_c_;
na=na+1;
end;%if (abs(tmp_val)>=-log(0.05));
end;%if (nvar0~=nvar1);
end;%for nvar1=0:n_var_retain-1;
end;%for nvar0=0:n_var_retain-1;
n_a = na;
tmp_x__ = tmp_x__(:,1:n_a);
tmp_y__ = tmp_y__(:,1:n_a);
tmp_c___ = tmp_c___(:,1:n_a,:);
hold on;
p = patch(tmp_x__,tmp_y__,tmp_c___);
set(p,'EdgeColor','k');
hold off;
%%%%;
n_side = 8; r_rad = 0.5/sqrt(2);
tmp_w_ = linspace(0,2*pi,1+n_side); n_w = numel(tmp_w_);
tmp_rcw_ = r_rad*cos(tmp_w_);
tmp_rsw_ = r_rad*sin(tmp_w_);
tmp_x__ = zeros(n_w,n_var_retain*(n_var_retain-1));
tmp_y__ = zeros(n_w,n_var_retain*(n_var_retain-1));
tmp_c___ = zeros(1,n_var_retain*(n_var_retain-1),3);
na=0;
for nvar0=0:n_var_retain-1;
for nvar1=0:n_var_retain-1;
if (nvar0~=nvar1);
tmp_x_ = 1.0 + nvar1 + tmp_rsw_;
tmp_y_ = 1.0 + nvar0 + tmp_rcw_;
tmp_val = X_sjnlp_zzzz__(1+nvar0,1+nvar1);
if (abs(tmp_val)>=-log(0.05));
tmp_val = (tmp_val - min(c_nlpvt2_lim_))/diff(c_nlpvt2_lim_);
nc_nlpvt2 = max(0,min(n_c_nlpvt2-1,floor(n_c_nlpvt2*tmp_val)));
tmp_c_ = c_nlpvt2__(1+nc_nlpvt2,:);
tmp_x__(:,1+na) = tmp_x_;
tmp_y__(:,1+na) = tmp_y_;
tmp_c___(1,1+na,:) = tmp_c_;
na=na+1;
end;%if (abs(tmp_val)>=-log(0.05));
end;%if (nvar0~=nvar1);
end;%for nvar1=0:n_var_retain-1;
end;%for nvar0=0:n_var_retain-1;
n_a = na;
tmp_x__ = tmp_x__(:,1:n_a);
tmp_y__ = tmp_y__(:,1:n_a);
tmp_c___ = tmp_c___(:,1:n_a,:);
hold on;
p = patch(tmp_x__,tmp_y__,tmp_c___);
set(p,'EdgeColor','k');
hold off;
%%%%;
r_rad = 0.5/sqrt(2);
tmp_x_plus__ = r_rad * [ [ +1 ; -1 ] , [  0 ;  0 ] ];
tmp_y_plus__ = r_rad * [ [  0 ;  0 ] , [ +1 ; -1 ] ];
tmp_z_plus__ = 1e-3 * ones(2,2);
tmp_x_subt__ = r_rad * [ [  0 ;  0 ] ];
tmp_y_subt__ = r_rad * [ [ +1 ; -1 ] ];
tmp_z_subt__ = 1e-3 * ones(2,1);
hold on;
for nvar0=0:n_var_retain-1;
for nvar1=0:n_var_retain-1;
if (nvar0~=nvar1);
tmp_val = X_sjnlp_zzzz__(1+nvar0,1+nvar1);
if (abs(tmp_val)>=-log(0.05));
if (tmp_val> 0); l=line( 1.0 + nvar1 + tmp_y_plus__ , 1 + nvar0 + tmp_x_plus__ , tmp_z_plus__ ); set(l,'Color','k','LineWidth',1); end;
if (tmp_val< 0); l=line( 1.0 + nvar1 + tmp_y_subt__ , 1 + nvar0 + tmp_x_subt__ , tmp_z_subt__ ); set(l,'Color','k','LineWidth',1); end;
end;%if (abs(tmp_val)>=-log(0.05));
end;%if (nvar0~=nvar1);
end;%for nvar1=0:n_var_retain-1;
end;%for nvar0=0:n_var_retain-1;
hold off;
%%%%;
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title(sprintf('signed %s p-value',tmp_str_title),'Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
%%%%;
sgtitle(sprintf('%s',fname_fig),'Interpreter','none');
colormap(c_nlpvt2__);
cb = colorbar();
cb_lim_ = get(cb,'limits');
tmp_tick_ = linspace(min(cb_lim_),max(cb_lim_),numel(-27:3:+27));
set(cb,'TickLength',[0.0],'Ticks',tmp_tick_,'TickLabels',[-27:3:+27]);
disp(sprintf(' %% writing %s',fname_fig));
%%%%%%%%;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for ntype=0:2-1;

flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_%s_shuffle_A_sjnlp_FIGA',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 12; np=0;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1); c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1); c_nlpvt2_lim_ = [-27,+27];
%%%%%%%%;
imagesc(clear_diag_0(A_sqr_sjnlp___(:,:,1+0)),c_nlpvt2_lim_); axis image;
l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ p-value (holm)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var_retain)*[zeros(1,n_var_retain);ones(1,n_var_retain)];
l_y_ = +0.5 + repmat([0:n_var_retain-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
sgtitle(sprintf('%s',fname_fig),'Interpreter','none');
colormap(c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_%s_shuffle_B_sjnlp_FIGB',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 12; np=0;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1); c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1); c_nlpvt2_lim_ = [-27,+27];
%%%%%%%%;
imagesc(clear_diag_0(BB_inv_sqr_sjnlp___(:,:,1+0)),c_nlpvt2_lim_); axis image;
l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ p-value (holm)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var_retain)*[zeros(1,n_var_retain);ones(1,n_var_retain)];
l_y_ = +0.5 + repmat([0:n_var_retain-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
sgtitle(sprintf('%s',fname_fig),'Interpreter','none');
colormap(c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_%s_shuffle_AB2_FIGAB2',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 12;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(1,2,np);np=np+1; imagesc(clear_diag_0(A_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ p-value (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var_retain)*[zeros(1,n_var_retain);ones(1,n_var_retain)];
l_y_ = +0.5 + repmat([0:n_var_retain-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
subplot_(1,2) = subplot(1,2,np);np=np+1; imagesc(clear_diag_0(BB_inv_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ p-value (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var_retain)*[zeros(1,n_var_retain);ones(1,n_var_retain)];
l_y_ = +0.5 + repmat([0:n_var_retain-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_zzzz),'Interpreter','none');
colormap(subplot_(1,1),c_nlpvt2__);
colormap(subplot_(1,2),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_%s_shuffle_A_FIGA',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;set(gcf,'Position',1+[0,0,1024,1024]);fontsize_use = 12;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(1,1,np);np=np+1; imagesc(clear_diag_0(A_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ z-score (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var_retain)*[zeros(1,n_var_retain);ones(1,n_var_retain)];
l_y_ = +0.5 + repmat([0:n_var_retain-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_zzzz),'Interpreter','none');
colormap(subplot_(1,1),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_%s_shuffle_B_FIGB',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;set(gcf,'Position',1+[0,0,1024,1024]);fontsize_use = 12;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(1,1,np);np=np+1; imagesc(clear_diag_0(BB_inv_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ z-score (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var_retain)*[zeros(1,n_var_retain);ones(1,n_var_retain)];
l_y_ = +0.5 + repmat([0:n_var_retain-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_zzzz),'Interpreter','none');
colormap(subplot_(1,1),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_%s_shuffle_C_FIGC',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;set(gcf,'Position',1+[0,0,1024,1024]);fontsize_use = 12;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(1,1,np);np=np+1; imagesc(clear_diag_0(CC_inv_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $CC^{T}$ z-score (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var_retain)*[zeros(1,n_var_retain);ones(1,n_var_retain)];
l_y_ = +0.5 + repmat([0:n_var_retain-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_zzzz),'Interpreter','none');
colormap(subplot_(1,1),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_%s_shuffle_aA_FIGaA',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 8;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(2,3,np);np=np+1; 
tmp_x_ = zeros(4,n_var_retain); tmp_y_ = zeros(4,n_var_retain); tmp_c_ = zeros(1,n_var_retain,3);
for nvar=0:n_var_retain-1; 
tmp_y = a_avg_zzzz_(1+nvar); tmp_x_(:,1+nvar) = -0.5+1+nvar+[0;0;1;1]; tmp_y_(:,1+nvar) = [0;tmp_y;tmp_y;0]; 
tmp_nc = max(0,min(n_c_beach-1,floor(n_c_beach*(tmp_y-min(c_beach_lim_))/diff(c_beach_lim_))));
tmp_c_(1,1+nvar,:) = c_beach__(1+tmp_nc,:);
end;%for nvar=0:n_var_retain-1; 
patch(tmp_x_,tmp_y_,tmp_c_,'EdgeColor','k');
xlim([0,n_var_retain+1]);ylim(0.125*c_beach_lim_);
%set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source/target'); ylabel('magnitude'); title('$a$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
subplot_(1,2) = subplot(2,3,np);np=np+1; 
tmp_x_ = zeros(4,n_var_retain); tmp_y_ = zeros(4,n_var_retain); tmp_c_ = zeros(1,n_var_retain,3);
for nvar=0:n_var_retain-1; 
tmp_y = a_std_zzzz_(1+nvar); tmp_x_(:,1+nvar) = -0.5+1+nvar+[0;0;1;1]; tmp_y_(:,1+nvar) = [0;tmp_y;tmp_y;0]; 
tmp_nc = max(0,min(n_c_beach-1,floor(n_c_beach*(tmp_y-min(c_beach_lim_))/diff(c_beach_lim_))));
tmp_c_(1,1+nvar,:) = c_beach__(1+tmp_nc,:);
end;%for nvar=0:n_var_retain-1; 
patch(tmp_x_,tmp_y_,tmp_c_,'EdgeColor','k');
xlim([0,n_var_retain+1]);ylim(0.125*c_beach_lim_);
%set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source/target'); ylabel('magnitude'); title('$a$ std','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
subplot_(1,3) = subplot(2,3,np);np=np+1;
tmp_x_ = zeros(4,n_var_retain); tmp_y_ = zeros(4,n_var_retain); tmp_c_ = zeros(1,n_var_retain,3);
for nvar=0:n_var_retain-1; 
tmp_y = a_snlp_zzzz_(1+nvar); tmp_x_(:,1+nvar) = -0.5+1+nvar+[0;0;1;1]; tmp_y_(:,1+nvar) = [0;tmp_y;tmp_y;0]; 
tmp_nc = max(0,min(n_c_nlpvt2-1,floor(n_c_nlpvt2*(tmp_y-min(c_nlpvt2_lim_))/diff(c_nlpvt2_lim_))));
tmp_c_(1,1+nvar,:) = c_nlpvt2__(1+tmp_nc,:);
end;%for nvar=0:n_var_retain-1; 
patch(tmp_x_,tmp_y_,tmp_c_,'EdgeColor','k');
xlim([0,n_var_retain+1]);ylim(c_nlpvt2_lim_);
%set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source/target'); ylabel('magnitude'); title('signed $a$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
%%%%%%%%;
subplot_(2,1) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(A_avg_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$A$ avg','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,2) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(A_std_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$A$ std','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,3) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(A_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_zzzz),'Interpreter','none');
colormap(subplot_(2,1),colormap_beach);
colormap(subplot_(2,2),colormap_beach);
colormap(subplot_(2,3),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_%s_shuffle_BC_FIGBC',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 8;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(BB_inv_avg_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$BB^{T}$ avg','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(1,2) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(BB_inv_std_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$BB^{T}$ std','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(1,3) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(BB_inv_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(2,1) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(CC_inv_avg_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$CC^{T}$ avg','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,2) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(CC_inv_std_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$CC^{T}$ std','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,3) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(CC_inv_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $CC^{T}$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_zzzz),'Interpreter','none');
colormap(subplot_(1,1),colormap_beach);
colormap(subplot_(1,2),colormap_beach);
colormap(subplot_(1,3),c_nlpvt2__);
colormap(subplot_(2,1),colormap_beach);
colormap(subplot_(2,2),colormap_beach);
colormap(subplot_(2,3),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_%s_shuffle_AB_FIGAB',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 8;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(2,2,np);np=np+1; imagesc(clear_diag_0(A_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(1,2) = subplot(2,2,np);np=np+1; imagesc(clear_diag_0(BB_inv_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(2,1) = subplot(2,2,np);np=np+1; imagesc(clear_diag_0(A_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ z-score (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(2,2) = subplot(2,2,np);np=np+1; imagesc(clear_diag_0(BB_inv_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var_retain); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var_retain,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var_retain,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ z-score (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_zzzz),'Interpreter','none');
colormap(subplot_(1,1),c_nlpvt2__);
colormap(subplot_(1,2),c_nlpvt2__);
colormap(subplot_(2,1),c_nlpvt2__);
colormap(subplot_(2,2),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% check distribution of shuffle coefficients. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_single;
%%%%%%%%;
fname_mat_zzzz_s = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix_zzzz);
if ( exist(fname_mat_zzzz_s,'file'));
n_found_zzzz = size(tmp_zzzz_s_.L_cmb_prm_,1) - 1;
a_zzzz__ = tmp_zzzz_s_.a_cmb_prm__(1+index_var_retain_,2:end);
A_zzzz__ = reshape(tmp_zzzz_s_.A_cmb_prm___(1+index_var_retain_,1+index_var_retain_,2:end),[n_var_retain^2,n_found_zzzz]);
BB_inv_zzzz__ = reshape(tmp_zzzz_s_.BB_inv_cmb_prm___(1+index_var_retain_,1+index_var_retain_,2:end),[n_var_retain^2,n_found_zzzz]);
CC_inv_zzzz__ = reshape(tmp_zzzz_s_.CC_inv_cmb_prm___(1+index_var_retain_,1+index_var_retain_,2:end),[n_var_retain^2,n_found_zzzz]);
disp(sprintf(' %% n_found_zzzz %d',n_found_zzzz'));
a_zzzz_nrm__ = (a_zzzz__ - repmat(mean(a_zzzz__,2),[1,n_found_zzzz]))./repmat(std(a_zzzz__,[],2),[1,n_found_zzzz]);
A_zzzz_nrm__ = (A_zzzz__ - repmat(mean(A_zzzz__,2),[1,n_found_zzzz]))./repmat(std(A_zzzz__,[],2),[1,n_found_zzzz]);
BB_inv_zzzz_nrm__ = (BB_inv_zzzz__ - repmat(mean(BB_inv_zzzz__,2),[1,n_found_zzzz]))./repmat(std(BB_inv_zzzz__,[],2),[1,n_found_zzzz]);
CC_inv_zzzz_nrm__ = (CC_inv_zzzz__ - repmat(mean(CC_inv_zzzz__,2),[1,n_found_zzzz]))./repmat(std(CC_inv_zzzz__,[],2),[1,n_found_zzzz]);
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_shuffle_coefficient_distribution',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 8;
subplot(2,2,1); hist(a_zzzz_nrm__(:),128); title('$a$','Interpreter','latex');
subplot(2,2,2); hist(A_zzzz_nrm__(:),128); title('$A$','Interpreter','latex');
subplot(2,2,3); hist(BB_inv_zzzz_nrm__(:),128); title('$BB^{T}$','Interpreter','latex');
subplot(2,2,4); hist(CC_inv_zzzz_nrm__(:),128); title('$CC^{T}$','Interpreter','latex');
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if ( exist(fname_mat_zzzz_s,'file'));
%%%%%%%%;
end;%if flag_single;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_exist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (verbose); disp(sprintf(' %% [finished dolphin_d00_cmb_collect_and_plot_4]')); end;
