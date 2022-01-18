% dolphin_dXX_load_7;

%%%%%%%%;
% Now represent a 2-variable case-study using res_ind. ;
% Here we plot A__. ;
%%%%%%%%;
nvar_pair_ = ...
{ ...
  { 'SED60' , 'CPK' } ...
 ,{ 'SED60' , 'Albumin' } ...
 ,{ 'SED60' , 'HGB' } ...
 ,{ 'CPK' , 'Creatinine' } ...
 ,{ 'GFR' , 'Creatinine' } ...
 ,{ 'Iron' , 'MCH' } ...
 ,{ 'Bilirubin' , 'MCHC' } ...
 ,{ 'Bilirubin' , 'HCT' } ...
 ,{ 'CPK' , 'Creatinine' } ...
 ,{ 'HGB' , 'MCH' } ...
 ,{ 'RBCDist' , 'MCH' } ...
 ,{ 'MPV' , 'UricAcid' } ...
 ,{ 'RBC' , 'Platelets' } ...
};
n_pair = numel(nvar_pair_);
%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npair=0:n_pair-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

str_pair_0 = nvar_pair_{1+npair}{1+0};
str_pair_1 = nvar_pair_{1+npair}{1+1};
nvar_ori_0 = efind(strcmp(string_dat_ori_name_,str_pair_0));
nvar_ori_1 = efind(strcmp(string_dat_ori_name_,str_pair_1));
res_ind_nrm_sub__ = res_ind_nrm__(:,1+[nvar_ori_0,nvar_ori_1]);
n_var_sub = 2;
%%%%%%%%;
n_shuffle = 0;
a_ind_sub_prm__ = zeros(n_var_sub,1+n_shuffle);
A_ind_sub_prm___ = zeros(n_var_sub,n_var_sub,1+n_shuffle);
BB_inv_ind_sub_prm___ = zeros(n_var_sub,n_var_sub,1+n_shuffle);
CC_inv_ind_sub_prm___ = zeros(n_var_sub,n_var_sub,1+n_shuffle);
for nshuffle=0:n_shuffle;
if (mod(nshuffle,16)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
res_ind_nrm_sub_prm__ = res_ind_nrm_sub__;
if (nshuffle>0);
[res_ind_nrm_sub_prm__] = dolphin_permute_0(aid_ori_,age_ori_,res_ind_nrm_sub__,nshuffle);
end;%if (nshuffle>0);
parameter = struct('type','parameter');
[ ...
 parameter ...
,a_ind_sub_prm_ ...
,A_ind_sub_prm__ ...
,BB_inv_ind_sub_prm__ ...
,CC_inv_ind_sub_prm__ ...
,L_ind_sub_prm ...
,niteration_ind_sub_prm ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ori_ ...
,age_ori_ ...
,res_ind_nrm_sub_prm__ ...
);
a_ind_sub_prm__(:,1+nshuffle) = a_ind_sub_prm_;
A_ind_sub_prm___(:,:,1+nshuffle) = A_ind_sub_prm__;
BB_inv_ind_sub_prm___(:,:,1+nshuffle) = BB_inv_ind_sub_prm__;
CC_inv_ind_sub_prm___(:,:,1+nshuffle) = CC_inv_ind_sub_prm__;
end;%for nshuffle=0:n_shuffle;
%%%%%%%%;
a_ind_sub_avg_ = mean(a_ind_sub_prm__(:,2:end),2);
a_ind_sub_std_ = std(a_ind_sub_prm__(:,2:end),1,2);
A_ind_sub_avg__ = mean(A_ind_sub_prm___(:,:,2:end),3);
A_ind_sub_std__ = std(A_ind_sub_prm___(:,:,2:end),1,3);
BB_inv_ind_sub_avg__ = mean(BB_inv_ind_sub_prm___(:,:,2:end),3);
BB_inv_ind_sub_std__ = std(BB_inv_ind_sub_prm___(:,:,2:end),1,3);
CC_inv_ind_sub_avg__ = mean(CC_inv_ind_sub_prm___(:,:,2:end),3);
CC_inv_ind_sub_std__ = std(CC_inv_ind_sub_prm___(:,:,2:end),1,3);
a_ind_sub_z__ = (a_ind_sub_prm__ - repmat(a_ind_sub_avg_,[1,1+n_shuffle]))./repmat(a_ind_sub_std_,[1,1+n_shuffle]);
a_ind_sub_nlp__ = -z_to_p_twosided_0(a_ind_sub_z__);
a_ind_sub_snlp__ = sign(a_ind_sub_prm__ - repmat(a_ind_sub_avg_,[1,1+n_shuffle])).*a_ind_sub_nlp__;
A_ind_sub_z___ = (A_ind_sub_prm___ - repmat(A_ind_sub_avg__,[1,1,1+n_shuffle]))./repmat(A_ind_sub_std__,[1,1,1+n_shuffle]);
A_ind_sub_nlp___ = -z_to_p_twosided_0(A_ind_sub_z___);
A_ind_sub_snlp___ = sign(A_ind_sub_prm___ - repmat(A_ind_sub_avg__,[1,1,1+n_shuffle])).*A_ind_sub_nlp___;
BB_inv_ind_sub_z___ = (BB_inv_ind_sub_prm___ - repmat(BB_inv_ind_sub_avg__,[1,1,1+n_shuffle]))./repmat(BB_inv_ind_sub_std__,[1,1,1+n_shuffle]);
BB_inv_ind_sub_nlp___ = -z_to_p_twosided_0(BB_inv_ind_sub_z___);
BB_inv_ind_sub_snlp___ = sign(BB_inv_ind_sub_prm___ - repmat(BB_inv_ind_sub_avg__,[1,1,1+n_shuffle])).*BB_inv_ind_sub_nlp___;
CC_inv_ind_sub_z___ = (CC_inv_ind_sub_prm___ - repmat(CC_inv_ind_sub_avg__,[1,1,1+n_shuffle]))./repmat(CC_inv_ind_sub_std__,[1,1,1+n_shuffle]);
CC_inv_ind_sub_nlp___ = -z_to_p_twosided_0(CC_inv_ind_sub_z___);
CC_inv_ind_sub_snlp___ = sign(CC_inv_ind_sub_prm___ - repmat(CC_inv_ind_sub_avg__,[1,1,1+n_shuffle])).*CC_inv_ind_sub_nlp___;
%%%%%%%%;
%%%%;
%%%%%%%%;
a_ind_sub_ = a_ind_sub_prm__(:,1);
A_ind_sub__ = A_ind_sub_prm___(:,:,1);
BB_inv_ind_sub__ = BB_inv_ind_sub_prm___(:,:,1);
CC_inv_ind_sub__ = CC_inv_ind_sub_prm___(:,:,1);
%%%%;
n_step = 1;
n_a = 0;
for nu_aid_ori=0:n_u_aid_ori-1;
tmp_index_aid_ori_ = index_aid_ori__{1+nu_aid_ori};
tmp_res__ = res_ind_nrm_sub__(1+tmp_index_aid_ori_,:);
tmp_age_ = age_ori_(1+tmp_index_aid_ori_);
tmp_n_a = max(0, numel( efind( isfinite(tmp_age_) & isfinite(sum(tmp_res__,2)) ) ) - 1);
n_a = n_a + tmp_n_a*n_step;
end;%for nu_aid_ori=0:n_u_aid_ori-1;

%%%%;
tmp_dt_all_ = zeros(n_a,1);
tmp_dY_all__ = zeros(n_var_sub,n_a);
sum_1 = 0;
sum_dt = 0;
sum_dtdt = 0;
sum_DD_inv_j__ = zeros(n_var_sub,n_var_sub);
sum_DD_inv_j_dt__ = zeros(n_var_sub,n_var_sub);
l2_DD_inv__ = zeros(n_var_sub,n_var_sub);
l2_BC__ = zeros(n_var_sub,n_var_sub);
na=0;
for nu_aid_ori=0:n_u_aid_ori-1;
tmp_index_aid_ori_ = index_aid_ori__{1+nu_aid_ori};
tmp_res__ = res_ind_nrm_sub__(1+tmp_index_aid_ori_,:);
tmp_age_ = age_ori_(1+tmp_index_aid_ori_);
index_use_ = efind( isfinite(tmp_age_) & isfinite(sum(tmp_res__,2)) );
tmp_res_use__ = tmp_res__(1+index_use_,:);
tmp_age_use_ = tmp_age_(1+index_use_);
[tmp_age_sort_,tmp_index_age_srt_] = sort(tmp_age_use_,'ascend'); tmp_index_age_srt_ = tmp_index_age_srt_-1;
n_age = numel(tmp_index_age_srt_);
for nstep=0:n_step-1;
if (n_age>n_step+1);
tmp_dt_sort_ = tmp_age_sort_(1+(1+nstep):end)-tmp_age_sort_(1:end-(1+nstep));
for nage=(1+nstep):n_age-1;
tmp_dt = tmp_age_use_(1+tmp_index_age_srt_(1+nage+0)) - tmp_age_use_(1+tmp_index_age_srt_(1+nage-(1+nstep)));
assert(fnorm(tmp_dt_sort_(1+nage-(1+nstep))-tmp_dt)<1e-6);
tmp_dt_all_(1+na) = tmp_dt;
Y_pre_ = reshape(tmp_res_use__(1+tmp_index_age_srt_(1+nage-(1+nstep)),:),[n_var_sub,1]);
Y_pos_ = reshape(tmp_res_use__(1+tmp_index_age_srt_(1+nage+0),:),[n_var_sub,1]);
dY_ = (Y_pos_ - Y_pre_) - tmp_dt*(a_ind_sub_ + A_ind_sub__*Y_pre_);
tmp_dY_all__(:,1+na) = dY_;
DD_inv_j__ = dY_*transpose(dY_);
sum_1 = sum_1 + 1;
sum_dt = sum_dt + tmp_dt;
sum_dtdt = sum_dtdt + tmp_dt.^2;
sum_DD_inv_j__ = sum_DD_inv_j__ + DD_inv_j__;
sum_DD_inv_j_dt__ = sum_DD_inv_j_dt__ + DD_inv_j__*tmp_dt;
l2_DD_inv__ = l2_DD_inv__ + DD_inv_j__.^2;
na=na+1;
end;%for nage=1:n_age-1;
end;%if (n_age>n_step+1);
end;%for nstep=0:n_step-1;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
%%%%;
tmp_BB_inv_ind_sub__ = zeros(n_var_sub,n_var_sub);
tmp_CC_inv_ind_sub__ = zeros(n_var_sub,n_var_sub);
tmp_LHS__ = [ sum_dtdt , sum_dt ; sum_dt , sum_1 ];
for nvar0=0:n_var_sub-1;
for nvar1=nvar0:n_var_sub-1;
tmp_RHS_ = [ sum_DD_inv_j_dt__(1+nvar0,1+nvar1) ; sum_DD_inv_j__(1+nvar0,1+nvar1) ];
tmp_BC_ = tmp_LHS__ \ tmp_RHS_;
tmp_BB_inv_ind_sub__(1+nvar0,1+nvar1) = tmp_BC_(1+0);
tmp_CC_inv_ind_sub__(1+nvar0,1+nvar1) = tmp_BC_(1+1);
tmp_BB_inv_ind_sub__(1+nvar1,1+nvar0) = tmp_BB_inv_ind_sub__(1+nvar0,1+nvar1);
tmp_CC_inv_ind_sub__(1+nvar1,1+nvar0) = tmp_CC_inv_ind_sub__(1+nvar0,1+nvar1);
end;%for nvar1=nvar0:n_var_sub-1;
end;%for nvar0=0:n_var_sub-1;
%%%%;

%%%%%%%%;
fname_fig = sprintf('%s/dolphin_res_ind_nrm_sub_dt_nvar%.2d_%s_nvar%.2d_%s_FIGH',dir_jpg_base,nvar_ori_0,string_dat_ori_name_{1+nvar_ori_0},nvar_ori_1,string_dat_ori_name_{1+nvar_ori_1});
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;
p_row = 2; p_col = 5; ns=0;
linewidth_use = 8;
markersize_use = 4;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
prctile_ = linspace(0,100,1+p_row*p_col);
tmp_dt_all_pXX_ = prctile(tmp_dt_all_,prctile_);
ylim_ = prctile(tmp_dY_all__,[ 1,99],'all');
ylim_ = mean(ylim_) + 1.25*0.5*diff(ylim_)*[-1,+1];
%%%%%%%%;
for ns=0:p_row*p_col-1;
subplot(p_row,p_col,1+ns);cla;
tmp_dt_all_min = tmp_dt_all_pXX_(1+ns+0);
tmp_dt_all_max = tmp_dt_all_pXX_(1+ns+1);
tmp_DD_inv_min__ = tmp_dt_all_min*BB_inv_ind_sub__ + CC_inv_ind_sub__;
tmp_DD_inv_max__ = tmp_dt_all_max*BB_inv_ind_sub__ + CC_inv_ind_sub__;
tmp_index_use_ = efind( (tmp_dt_all_>=tmp_dt_all_min) & (tmp_dt_all_< tmp_dt_all_max) );
tmp_dt_all_use_ = tmp_dt_all_(1+tmp_index_use_);
tmp_dY_all_use__ = tmp_dY_all__(:,1+tmp_index_use_);
patch([ylim_,flip(ylim_)],[ylim_(1),ylim_,ylim_(2)],-1*ones(1,4),0.85*[1,1,1]);
dolphin_ellipse_from_DD_inv_0(tmp_DD_inv_min__,[1,1,0.85],linewidth_use,2.0);
dolphin_ellipse_from_DD_inv_0(tmp_DD_inv_min__,[1,1,0.25],linewidth_use,1.0);
hold on;
plot(tmp_dY_all_use__(1,:),tmp_dY_all_use__(2,:),'ko','MarkerFaceColor',1.0*[0.8,0.8,1],'MarkerSize',markersize_use);
xlim(ylim_); ylim(ylim_);
xlabel(sprintf('%s increment',string_dat_ori_name_{1+nvar_ori_0}),'Interpreter','none');
ylabel(sprintf('%s increment',string_dat_ori_name_{1+nvar_ori_1}),'Interpreter','none');
title(sprintf('dt in [%.2d,%.2d] percentile',prctile_(1+ns+0),prctile_(1+ns+1)),'Interpreter','none');
grid on;
axis square;
hold off;
end;%for ns=0:p_row*p_col-1;
%%%%%%%%;
sgtitle(sprintf('%s vs %s',string_dat_ori_name_{1+nvar_ori_0},string_dat_ori_name_{1+nvar_ori_1}));
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for npair=0:n_pair-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;















