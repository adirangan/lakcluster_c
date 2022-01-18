clear;
dXX = 0; n_shuffle = 256; age_percentile_range_ = []; sex_constrain = ''; isU_constrain = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (isempty(dXX)); dXX = 0; end;
if (isempty(n_shuffle)); n_shuffle = 32; end;
if (isempty(sex_constrain)); sex_constrain = ''; end;
if (isempty(isU_constrain)); isU_constrain = ''; end;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_dolphin',string_root);
dir_jpg = sprintf('%s/dir_jpg_ind_7',dir_trunk);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
dir_mat = sprintf('%s/dir_mat_ind_7',dir_trunk);
if (~exist(dir_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_mat)); mkdir(dir_mat); end;
dir_shuffle_mat = sprintf('%s/dir_shuffle%d_mat_ind_7',dir_trunk,n_shuffle);
if (~exist(dir_shuffle_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_shuffle_mat)); mkdir(dir_shuffle_mat); end;

tolerance_master = 1e-6;
flag_replot=0;
date_diff_threshold = 0.5;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;

%%%%;
str_sex_constrain = '';
if ~isempty(sex_constrain);
str_sex_constrain = sprintf('_%s',sex_constrain);
end;%if ~isempty(sex_constrain);
%%%%;
str_isU_constrain = '';
if ~isempty(isU_constrain);
str_isU_constrain = sprintf('_U%s',isU_constrain);
end;%if ~isempty(isU_constrain);
%%%%;
str_age_constrain = '';
if (~isempty(age_percentile_range_));
str_age_constrain = sprintf('_age_%.2d%.2d',min(age_percentile_range_),min(99,max(age_percentile_range_)));
end;%if (~isempty(age_percentile_range_));
%%%%;
fname_infix = ...
sprintf( ...
 'aid%.2d%s%s%s' ...
,dXX ...
,str_age_constrain ...
,str_sex_constrain ...
,str_isU_constrain ...
);
disp(sprintf(' %% fname_infix %s',fname_infix));

dolphin_dXX_load_7;
string_dat_name_ = string_dat_ori_name_;
n_var = n_var_ori;
isM_ = isM_ori_;
isF_ = isF_ori_;
isU_bbar_lo_ = isU_bbar_lo_ori_;
isU_bbar_hi_ = isU_bbar_hi_ori_;
aid_ = aid_ori_;
age_ = age_ori_;
dat__ = res_ind_nrm__;
tmp_index_use_ = efind( 1 ...
&isfinite(age_) ...
& (sum(isfinite(dat__),2)> 0) ...
);
isM_ = isM_(1+tmp_index_use_);
isF_ = isF_(1+tmp_index_use_);
isU_bbar_lo_ = isU_bbar_lo_(1+tmp_index_use_);
isU_bbar_hi_ = isU_bbar_hi_(1+tmp_index_use_);
aid_ = aid_(1+tmp_index_use_);
age_ = age_(1+tmp_index_use_);
dat__ = dat__(1+tmp_index_use_,:);

%%%%%%%%;
% select sex. ;
%%%%%%%%;
if ~isempty(sex_constrain);
if strcmp(sex_constrain,'M'); tmp_index_ = efind(isM_); end;%if strcmp(sex_constrain,'M');
if strcmp(sex_constrain,'F'); tmp_index_ = efind(isF_); end;%if strcmp(sex_constrain,'F');
isM_ = isM_(1+tmp_index_);
isF_ = isF_(1+tmp_index_);
isU_bbar_lo_ = isU_bbar_lo_(1+tmp_index_);
isU_bbar_hi_ = isU_bbar_hi_(1+tmp_index_);
aid_ = aid_(1+tmp_index_);
age_ = age_(1+tmp_index_);
dat__ = dat__(1+tmp_index_,:);
end;%if ~isempty(sex_constrain);

%%%%%%%%;
% select isU. ;
%%%%%%%%;
if ~isempty(isU_constrain);
if strcmp(isU_constrain,'lo'); tmp_index_ = efind(isU_bbar_lo_); end;%if strcmp(isU_constrain,'lo');
if strcmp(isU_constrain,'hi'); tmp_index_ = efind(isU_bbar_hi_); end;%if strcmp(isU_constrain,'hi');
isM_ = isM_(1+tmp_index_);
isF_ = isF_(1+tmp_index_);
isU_bbar_lo_ = isU_bbar_lo_(1+tmp_index_);
isU_bbar_hi_ = isU_bbar_hi_(1+tmp_index_);
aid_ = aid_(1+tmp_index_);
age_ = age_(1+tmp_index_);
dat__ = dat__(1+tmp_index_,:);
end;%if ~isempty(isU_constrain);

%%%%%%%%;
% select age range. ;
%%%%%%%%;
if (~isempty(age_percentile_range_));
tmp_index_ = efind(age_>=prctile(age_,min(age_percentile_range_)) & age_<=prctile(age_,max(age_percentile_range_)));
disp(sprintf(' %% age_percentile_range_ [%0.2f %0.2f] --> [%0.2f %0.2f] <-- %d found',age_percentile_range_,prctile(age_,age_percentile_range_),numel(tmp_index_)));
isM_ = isM_(1+tmp_index_);
isF_ = isF_(1+tmp_index_);
isU_bbar_lo_ = isU_bbar_lo_(1+tmp_index_);
isU_bbar_hi_ = isU_bbar_hi_(1+tmp_index_);
aid_ = aid_(1+tmp_index_);
age_ = age_(1+tmp_index_);
dat__ = dat__(1+tmp_index_,:);
end;%if (~isempty(age_percentile_range_));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now repeat for several permuted trials. ;
%%%%%%%%;
nvar0_use = efind(strcmp(string_dat_ori_name_,'Mg'));
nvar1_use = efind(strcmp(string_dat_ori_name_,'GFR'));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
%%%%%%%%;
a_cmb_prm__ = zeros(n_var,1+n_shuffle);
A_cmb_prm___ = zeros(n_var,n_var,1+n_shuffle);
BB_inv_cmb_prm___ = zeros(n_var,n_var,1+n_shuffle);
CC_inv_cmb_prm___ = zeros(n_var,n_var,1+n_shuffle);
L_cmb_prm_ = zeros(1+n_shuffle,1);
a_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
A_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
BB_inv_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
CC_inv_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
L_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
for nshuffle=0:1;%for nshuffle=0:n_shuffle;
if (nshuffle==0); dat_prm__ = dat__; end;
if (nshuffle> 0); dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle); end;
if (mod(nshuffle,32)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
%%%%%%%%;
a_cmb_ = zeros(n_var,1);
A_cmb__ = zeros(n_var,n_var);
BB_inv_cmb__ = zeros(n_var,n_var);
CC_inv_cmb__ = zeros(n_var,n_var);
L_cmb = 0;
a_each__ = zeros(n_var,n_var);
A_each__ = zeros(n_var,n_var);
BB_inv_each__ = zeros(n_var,n_var);
CC_inv_each__ = zeros(n_var,n_var);
L_each__ = zeros(n_var,n_var);
for nvar0=nvar0_use;%for nvar0=0:n_var-1;
for nvar1=nvar1_use;%for nvar1=nvar0+1:n_var-1;
dat_pair__ = dat_prm__(:,1+[nvar0,nvar1]);
[ ...
 parameter ...
,a_pair_ ...
,A_pair__ ...
,BB_inv_pair__ ...
,CC_inv_pair__ ...
,L_pair ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_pair__ ...
);
a_cmb_(1+[nvar0;nvar1]) = a_cmb_(1+[nvar0;nvar1]) + a_pair_(:);
A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + A_pair__;
BB_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = BB_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + BB_inv_pair__;
CC_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = CC_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + CC_inv_pair__;
L_cmb = L_cmb + L_pair;
a_each__(1+nvar0,1+nvar1) = a_pair_(1+0); a_each__(1+nvar1,1+nvar0) = a_pair_(1+1);
A_each__(1+nvar0,1+nvar1) = A_pair__(1+0,1+0); A_each__(1+nvar1,1+nvar0) = A_pair__(1+1,1+1);
BB_inv_each__(1+nvar0,1+nvar1) = BB_inv_pair__(1+0,1+0); BB_inv_each__(1+nvar1,1+nvar0) = BB_inv_pair__(1+1,1+1);
CC_inv_each__(1+nvar0,1+nvar1) = CC_inv_pair__(1+0,1+0); CC_inv_each__(1+nvar1,1+nvar0) = CC_inv_pair__(1+1,1+1);
L_each__(1+nvar0,1+nvar1) = L_pair; L_each__(1+nvar1,1+nvar0) = L_pair;
end;%for nvar1=nvar0:n_var-1;
end;%for nvar0=0:n_var-1;
a_cmb_prm__(:,1+nshuffle) = a_cmb_;
A_cmb_prm___(:,:,1+nshuffle) = A_cmb__;
BB_inv_cmb_prm___(:,:,1+nshuffle) = BB_inv_cmb__;
CC_inv_cmb_prm___(:,:,1+nshuffle) = CC_inv_cmb__;
L_cmb_prm_(1+nshuffle) = L_cmb;
a_each_prm___(:,:,1+nshuffle) = a_each__;
A_each_prm___(:,:,1+nshuffle) = A_each__;
BB_inv_each_prm___(:,:,1+nshuffle) = BB_inv_each__;
CC_inv_each_prm___(:,:,1+nshuffle) = CC_inv_each__;
L_each_prm___(:,:,1+nshuffle) = L_each__;
end;%for nshuffle=0:n_shuffle;

%{
dat_out__ = dolphin_impute_aA_0(aid_,age_,dat_pair__);
tmp_dat_out__ = fliplr(dolphin_impute_aA_0(aid_,age_,fliplr(dat_pair__)));
disp(sprintf(' %% ~isfinite: %d',numel(efind(~isfinite(dat_out__(:))))));
disp(sprintf(' %% ~isfinite: %d',numel(efind(~isfinite(tmp_dat_out__(:))))));
plot(dat_out__(:,2),tmp_dat_out__(:,2),'.');
%}
