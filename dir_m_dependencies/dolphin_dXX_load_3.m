fname_data = sprintf('%s/Dolphin_Data_201918755_s2.xlsx',dir_trunk);
T_ = readtable(fname_data);
aid_ori_ = T_{:,1};
age_ori_ = T_{:,4};
nvar_base = 7;
dat_ori__ = T_{:,(1+nvar_base):end};
n_smp_ori = size(dat_ori__,1);
n_var_ori = size(dat_ori__,2);
string_dat_ori_name_ = cell(n_var_ori,1);
for nvar=0:n_var_ori-1;
string_dat_ori_name_{1+nvar} = T_.Properties.VariableNames{1+nvar_base+nvar};
end;%for nvar=0:n_var_ori-1;

%%%%%%%%;
% log-normalizing SED60 and GGT. ;
%%%%%%%%;
nvar = efind(strcmp(string_dat_ori_name_,'SED60'));
dat_ori__(:,1+nvar) = log(1+dat_ori__(:,1+nvar));
string_dat_ori_name_{1+nvar} = sprintf('log(1+%s)',string_dat_ori_name_{1+nvar});
nvar = efind(strcmp(string_dat_ori_name_,'GGT'));
dat_ori__(:,1+nvar) = log(1+dat_ori__(:,1+nvar));
string_dat_ori_name_{1+nvar} = sprintf('log(1+%s)',string_dat_ori_name_{1+nvar});
%%%%%%%%;
% removing RBCDist and Albumin. ;
%%%%%%%%;
dat_ret__ = dat_ori__;
mss_var_threshold = 500;
index_var_rmv_ = efind(sum(~isfinite(dat_ret__),1)>mss_var_threshold);
index_var_rmv_ = union(index_var_rmv_,efind(strcmp(string_dat_ori_name_,'RBCDist')));
index_var_rmv_ = union(index_var_rmv_,efind(strcmp(string_dat_ori_name_,'Albumin')));
index_var_ret_ = setdiff(0:n_var_ori-1,index_var_rmv_);
dat_ret__ = dat_ret__(:,1+index_var_ret_);
n_var_ret = size(dat_ret__,2);
string_dat_ret_name_ = cell(n_var_ret,1);
for nvar=0:n_var_ret-1;
string_dat_ret_name_{1+nvar} = string_dat_ori_name_{1+index_var_ret_(1+nvar)};
end;%for nvar=0:n_var_ret-1;
%%%%%%%%;
% missingness threshold. ;
%%%%%%%%;
mss_smp_threshold = 7;
index_smp_ret_ = setdiff(0:n_smp_ori-1,efind(sum(~isfinite(dat_ret__),2)>mss_smp_threshold));
dat_ret__ = dat_ret__(1+index_smp_ret_,:);
n_smp_ret = size(dat_ret__,1);

aid_ret_ = aid_ori_(1+index_smp_ret_);
age_ret_ = age_ori_(1+index_smp_ret_);

aid_ = aid_ret_;
age_ = age_ret_;
dat_nrm__ = dat_ret__;
[n_smp,n_var] = size(dat_nrm__);
string_dat_name_ = string_dat_ret_name_;
prctile_threshold = 0.01;
for nvar=0:n_var-1;
tmp_index_ = efind(isfinite(dat_nrm__(:,1+nvar)));
tmp_lo = prctile(dat_nrm__(1+tmp_index_,1+nvar),100*prctile_threshold);
tmp_hi = prctile(dat_nrm__(1+tmp_index_,1+nvar),100*(1-prctile_threshold));
dat_nrm__(1+tmp_index_,1+nvar) = max(tmp_lo,min(tmp_hi,dat_nrm__(1+tmp_index_,1+nvar)));
end;%for nvar=0:n_var-1;

n_step = 1; dt_all_ = dolphin_dt_all_0(aid_,age_,n_step);
%%%%%%%%;
% correcting for drift. ;
%%%%%%%%;
[a_drift_,b_drift_] = dolphin_estimate_ab_0(aid_,age_,dat_nrm__);
%%%%;
u_aid_ = unique(aid_);
n_aid = numel(u_aid_);
n_aid_ = zeros(n_aid,1);
index_aid__ = cell(n_aid,1);
for naid=0:n_aid-1;
tmp_index_aid_ = efind(aid_==u_aid_(1+naid));
index_aid__{1+naid} = tmp_index_aid_;
n_aid_(1+naid) = numel(index_aid__{1+naid});
tmp_dat_nrm__ = dat_nrm__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
if (flag_drift0==1); tmp_dat_nrm__ = tmp_dat_nrm__ - (ones(numel(tmp_age_),1)*reshape(b_drift_,[1,n_var]) + tmp_age_*reshape(a_drift_,[1,n_var])); end;
if (flag_drift0==0); tmp_dat_nrm__ = tmp_dat_nrm__; end;
if (flag_drift0==0 & flag_center==1); tmp_dat_nrm__ = normalize(tmp_dat_nrm__,'zscore'); end;
if (flag_drift0==0 & flag_center==0); tmp_dat_nrm__ = tmp_dat_nrm__; end;
dat_nrm__(1+tmp_index_aid_,:) = tmp_dat_nrm__;
end;%for naid=0:n_aid-1;
if (flag_drift0==1 & flag_center==1); dat_nrm__ = normalize(dat_nrm__,'zscore'); end;
if (flag_drift0==1 & flag_center==0); dat_nrm__ = dat_nrm__; end;
%%%%;
n_step = 1; relative_variation = dolphin_relative_variation_0(aid_,age_,dat_nrm__,n_step);

%%%%%%%%;
% select dolphin. ;
%%%%%%%%;
if dXX>0; 
naid = efind(u_aid_==dXX);
if (numel(naid)==0); disp(sprintf(' %% dXX %d not found',dXX)); return; end;
disp(sprintf(' %% dXX %d, naid %d',dXX,naid));
tmp_index_aid_ = efind(aid_==u_aid_(1+naid));
aid_ = aid_(1+tmp_index_aid_);
age_ = age_(1+tmp_index_aid_);
dat_nrm__ = dat_nrm__(1+tmp_index_aid_,:);
end;%if dXX>0; 

dat__ = dat_nrm__;

