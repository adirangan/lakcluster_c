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
index_mss_ = efind(~isfinite(res_ind_nrm_sub__));
res_ind_nrm_sub_imp__ = dolphin_impute_aA_0(aid_ori_,age_ori_,res_ind_nrm_sub__,a_ind_sub_,A_ind_sub__,index_mss_);
%%%%%%%%;
[ ...
 parameter ...
,res_ind_nrm_sub_fil__ ...
,P_pos_nnv___ ...
] = ...
dolphin_estimate_kalman_filter_0( ...
 parameter ...
,aid_ori_ ...
,age_ori_ ...
,res_ind_nrm_sub_imp__ ...
,a_ind_sub_ ...
,A_ind_sub__ ...
,1*BB_inv_ind_sub__ ...
,1*CC_inv_ind_sub__ ...
);
[ ...
 parameter ...
,res_ind_nrm_sub_fil__ ...
,P_pos_nnv___ ...
] = ...
dolphin_estimate_kalman_filter_0( ...
 parameter ...
,aid_ori_ ...
,age_ori_ ...
,res_ind_nrm_sub_imp__ ...
,a_ind_sub_ ...
,A_ind_sub__ ...
,1*BB_inv_ind_sub__ ...
,1*CC_inv_ind_sub__ ...
,P_pos_nnv___ ...
);
%%%%%%%%;
%%%%;
%%%%%%%%;
parameter_gen = struct('type','parameter');
parameter_gen.tolerance_master = 1e-6;
parameter_gen.dt_avg = 1/64;
parameter_gen.T_max = age_upb;
parameter_gen.rseed = 0;
[ ...
 parameter_gen ...
,age_ind_gen_ ...
,X_ind_gen__ ...
,Y_ind_gen__ ...
] = ...
SDE_generate_data_0( ...
 parameter_gen ...
,a_ind_sub_ ...
,A_ind_sub__ ...
,[] ...
,BB_inv_ind_sub__ ...
,[] ...
,CC_inv_ind_sub__ ...
,zeros(2,1) ...
);
%%%%%%%%;
tmp_res_ind_gen__ = transpose(X_ind_gen__);
glim_ = prctile(X_ind_gen__,[ 1,99],'all');
glim_ = 0*mean(glim_) + 1.25*0.5*diff(glim_)*[-1,+1];
%%%%%%%%;
n_index_use_ = zeros(n_u_aid_ori,1);
for nu_aid_ori=0:n_u_aid_ori-1;
index_aid_ori_ = intersect( index_aid_ori__{1+nu_aid_ori} , index_smp_TT_ );
index_use_ = efind( isfinite(age_ori_(1+index_aid_ori_)) & isfinite(sum(res_ind_nrm_sub_fil__(1+index_aid_ori_,:),2)) );
n_index_use = numel(index_use_);
n_index_use_(1+nu_aid_ori) = n_index_use;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
[~,tmp_index_srt_] = sort(n_index_use_,'descend'); tmp_index_srt_ = tmp_index_srt_-1;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_res_ind_nrm_sub_fil_nvar%.2d_%s_nvar%.2d_%s_FIGF',dir_jpg_base,nvar_ori_0,string_dat_ori_name_{1+nvar_ori_0},nvar_ori_1,string_dat_ori_name_{1+nvar_ori_1});
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;
%p_row = 3; p_col = 6; p_pers = 1; ns=0;
%p_big_ = [ 0 , 1 , 2 , 1*p_col+0 , 1*p_col+1 , 1*p_col+2 , 2*p_col+0 , 2*p_col+1 , 2*p_col+2 ];
p_row = 2; p_col = 4; p_pers = 1; ns=0;
p_big_ = [ 0 , 1 , 1*p_col+0 , 1*p_col+1 ];
p_sml_ = setdiff([0:p_row*p_col-1],p_big_);
linewidth_use = 3;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
%%%%%%%%;
parameter_gen = struct('type','parameter');
parameter_gen.x0_lim_ = glim_;
parameter_gen.x1_lim_ = glim_;
subplot(p_row,p_col,1 + p_big_);cla;
vector_field_0(parameter_gen,A_ind_sub__); axisnotick; axis image;
hold on;
s=surfline_0(tmp_res_ind_gen__(:,1+0),tmp_res_ind_gen__(:,1+1),age_ind_gen_);
set(s,'LineWidth',linewidth_use);
hold off;
set(gca,'Clim',[age_lob,age_upb]);
xlabel(string_dat_ori_name_{1+nvar_ori_0});
ylabel(string_dat_ori_name_{1+nvar_ori_1});
title('SDE sample');
cb = colorbar;
set(cb,'TickLength',[0]);
set(cb,'Ticks',[age_lob,age_upb]);
%%%%%%%%;
for nl=0:n_u_aid_ori-1;
nu_aid_ori = tmp_index_srt_(1+nl);
if ns<numel(p_sml_);
subplot(p_row,p_col,1+p_sml_(1+ns));
%%%%%%%%;
index_aid_ori_ = intersect( index_aid_ori__{1+nu_aid_ori} , index_smp_TT_ );
index_use_ = efind( isfinite(age_ori_(1+index_aid_ori_)) & isfinite(sum(res_ind_nrm_sub_fil__(1+index_aid_ori_,:),2)) );
[tmp_age_srt_,tmp_ij] = sort(age_ori_(1+index_aid_ori_(1+index_use_)),'ascend');
tmp_res_ind_srt__ = res_ind_nrm_sub_fil__(1+index_aid_ori_(1+index_use_(tmp_ij)),:);
tmp_age_srt_ = tmp_age_srt_(2:end); tmp_res_ind_srt__ = tmp_res_ind_srt__(2:end,:);
%%%%;
xlim_ = prctile(tmp_res_ind_srt__(:,1+0),[5,95],'all');
xlim_ = mean(xlim_) + 1.25*0.5*diff(xlim_)*[-1,+1];
ylim_ = prctile(tmp_res_ind_srt__(:,1+1),[5,95],'all');
ylim_ = mean(ylim_) + 1.25*0.5*diff(ylim_)*[-1,+1];
hold on;
s=surfline_0(tmp_res_ind_srt__(:,1+0),tmp_res_ind_srt__(:,1+1),tmp_age_srt_);
set(s,'LineWidth',linewidth_use);
hold off;
set(gca,'Clim',[age_lob,age_upb]);
%%%%%%%%;
xlim(xlim_);xlabel(string_dat_ori_name_{1+nvar_ori_0});
ylim(ylim_);ylabel(string_dat_ori_name_{1+nvar_ori_1});
axis square; axis image;
title(sprintf('ID %d',u_aid_ori_(1+nu_aid_ori)));
grid on;
end;%if ns<numel(p_sml_)-1;
if (mod(nl,p_pers)==0); ns=ns+1; end;
end;%for nl=0:n_u_aid_ori-1;
fig80s;
sgtitle(sprintf('%s vs %s',string_dat_ori_name_{1+nvar_ori_0},string_dat_ori_name_{1+nvar_ori_1}));
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for npair=0:n_pair-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;















