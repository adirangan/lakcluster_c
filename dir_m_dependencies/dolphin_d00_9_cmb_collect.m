
clear;
dolphin_dXX_load_7;
dXX = 0;
n_shuffle = 256;
index_var_retain_ = setdiff(transpose(0:n_var_ori-1),efind(strcmp(string_dat_ori_name_,'GFR')));

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
flag_replot=1;
date_diff_threshold = 0.5;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;

for nU=0:2;
if nU==0; isU_constrain = ''; end;
if nU==1; isU_constrain = 'lo'; end;
if nU==2; isU_constrain = 'hi'; end;
for nS=0:2;
if nS==0; sex_constrain = ''; end;
if nS==1; sex_constrain = 'M'; end;
if nS==2; sex_constrain = 'F'; end;
for nA=0:3;
if nA==0; age_percentile_range_ = [      ]; end;
if nA==1 & isempty(sex_constrain); age_percentile_range_ = [ 0, 20]; end;
if nA==2 & isempty(sex_constrain); age_percentile_range_ = [20, 70]; end;
if nA==3 & isempty(sex_constrain); age_percentile_range_ = [70,100]; end;
if nA==1 & strcmp(sex_constrain,'M'); age_percentile_range_ = [ 0, 21]; end;
if nA==2 & strcmp(sex_constrain,'M'); age_percentile_range_ = [21, 73]; end;
if nA==3 & strcmp(sex_constrain,'M'); age_percentile_range_ = [73,100]; end;
if nA==1 & strcmp(sex_constrain,'F'); age_percentile_range_ = [ 0, 14]; end;
if nA==2 & strcmp(sex_constrain,'F'); age_percentile_range_ = [14, 68]; end;
if nA==3 & strcmp(sex_constrain,'F'); age_percentile_range_ = [68,100]; end;
%%%%%%%%;
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
%{
fname_infix = ...
sprintf( ...
 'aid%.2d%s%s%s' ...
,dXX ...
,str_age_constrain ...
,str_sex_constrain ...
,str_isU_constrain ...
);
disp(sprintf(' %% fname_infix %s',fname_infix));
%%%%;
parameter = struct('type','parameter');
parameter.flag_replot = flag_replot;
parameter.index_var_retain_ = index_var_retain_;
dolphin_d00_collect_and_plot_4( ...
 parameter ...
 ,string_dat_ori_name_ ...
 ,dir_mat ...
 ,dir_shuffle_mat ...
 ,dir_jpg ...
 ,fname_infix ...
);
%}
%%%%;
fname_infix = ...
sprintf( ...
 'aid%.2d%s%s%s_cmb' ...
,dXX ...
,str_age_constrain ...
,str_sex_constrain ...
,str_isU_constrain ...
);
disp(sprintf(' %% fname_infix %s',fname_infix));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_replot = flag_replot;
parameter.index_var_retain_ = index_var_retain_;
dolphin_d00_cmb_collect_and_plot_4( ...
 parameter ...
 ,string_dat_ori_name_ ...
 ,dir_mat ...
 ,dir_shuffle_mat ...
 ,dir_jpg ...
 ,fname_infix ...
);
%%%%%%%%;
end;%for nA=0:3;
end;%if nS=0:2;
end;%for nU=0:2;


na=0;
fname_infix_abz__ = cell(3,0);
%%%%;
fname_infix_abz__{1+0,1+na} = 'aid00_age_7099_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_2070_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_7099_vs_2070_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_age_7099_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_2070_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_7099_vs_2070_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_M_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_F_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_M_vs_F_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_age_7399_M_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_2173_M_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_7399_vs_2173_M_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_age_6899_F_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_1468_F_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_6899_vs_1468_F_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_Ulo_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_Uhi_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_Ulo_vs_Uhi_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_age_2070_Ulo_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_2070_Uhi_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_2070_Ulo_vs_Uhi_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_age_7099_Ulo_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_7099_Uhi_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_7099_Ulo_vs_Uhi_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_age_7399_M_Ulo_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_7399_M_Uhi_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_7399_M_Ulo_vs_Uhi_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_age_2173_M_Ulo_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_2173_M_Uhi_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_2173_M_Ulo_vs_Uhi_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_age_6899_F_Ulo_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_6899_F_Uhi_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_6899_F_Ulo_vs_Uhi_cmb';
na=na+1;
fname_infix_abz__{1+0,1+na} = 'aid00_age_1468_F_Ulo_cmb';
fname_infix_abz__{1+1,1+na} = 'aid00_age_1468_F_Uhi_cmb';
fname_infix_abz__{1+2,1+na} = 'aid00_age_1468_F_Ulo_vs_Uhi_cmb';
na=na+1;
%%%%;
n_a = na;
for na=0:n_a-1;
disp(sprintf(' %% na %d/%d --> %s %s --> %s',na,n_a,fname_infix_abz__{1+0,1+na},fname_infix_abz__{1+1,1+na},fname_infix_abz__{1+2,1+na}));
parameter = struct('type','parameter');
parameter.flag_replot = flag_replot;
parameter.index_var_retain_ = index_var_retain_;
dolphin_d00_cmb_collect_and_plot_4( ...
 parameter ...
 ,string_dat_ori_name_ ...
 ,dir_mat ...
 ,dir_shuffle_mat ...
 ,dir_jpg ...
 ,fname_infix_abz__{1+0,1+na} ...
 ,fname_infix_abz__{1+1,1+na} ...
 ,fname_infix_abz__{1+2,1+na} ...
);
end;%for na=0:n_a-1;

disp('returning'); return;

