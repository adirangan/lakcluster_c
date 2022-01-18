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
dolphin_d00_collect_and_plot_3( ...
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
if (verbose); disp(sprintf(' %% [entering dolphin_d00_collect_and_plot_3]')); end;

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
% parameter = []; string_dat_name_ = string_dat_ori_name_; dir_mat = '/home/rangan/dir_bcc/dir_dolphin/dir_mat_ind_7'; dir_shuffle_mat = '/home/rangan/dir_bcc/dir_dolphin/dir_shuffle256_mat_ind_7'; dir_jpg = '/home/rangan/dir_bcc/dir_dolphin/dir_jpg_ind_7'; fname_infix_aaaa = 'aid00_M'; fname_infix_bbbb = 'aid00_F'; fname_infix_zzzz = 'aid00_M_vs_F'; fname_infix_aaaa = 'aid00'; fname_infix_bbbb = ''; fname_infix_zzzz = '';

if (isempty(parameter)); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'flag_replot')); parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
flag_replot = parameter.flag_replot;


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
a_zzzz_o_ = tmp_zzzz_o_.a_; n_var = numel(a_zzzz_o_);
A_zzzz_o__ = tmp_zzzz_o_.A__;
BB_inv_zzzz_o__ = tmp_zzzz_o_.BB_inv__;
CC_inv_zzzz_o__ = tmp_zzzz_o_.CC_inv__;
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
n_found_zzzz = size(tmp_zzzz_s_.L_prm_,1) - 1;
n_shuffle = n_found_zzzz;
%%%%;
a_avg_zzzz_ = mean(tmp_zzzz_s_.a_prm__(:,2:end),2); a_std_zzzz_ = std(tmp_zzzz_s_.a_prm__(:,2:end),1,2);
A_avg_zzzz__ = mean(tmp_zzzz_s_.A_prm___(:,:,2:end),3); A_std_zzzz__ = std(tmp_zzzz_s_.A_prm___(:,:,2:end),1,3);
BB_inv_avg_zzzz__ = mean(tmp_zzzz_s_.BB_inv_prm___(:,:,2:end),3); BB_inv_std_zzzz__ = std(tmp_zzzz_s_.BB_inv_prm___(:,:,2:end),1,3);
CC_inv_avg_zzzz__ = mean(tmp_zzzz_s_.CC_inv_prm___(:,:,2:end),3); CC_inv_std_zzzz__ = std(tmp_zzzz_s_.CC_inv_prm___(:,:,2:end),1,3);
end;%if ( exist(fname_mat_zzzz_s,'file'));
a_nlp_zzzz_ = -z_to_p_twosided_0((a_zzzz_o_ - a_avg_zzzz_)./a_std_zzzz_);
A_nlp_zzzz__ = -z_to_p_twosided_0((A_zzzz_o__ - A_avg_zzzz__)./A_std_zzzz__);
BB_inv_nlp_zzzz__ = -z_to_p_twosided_0((BB_inv_zzzz_o__ - BB_inv_avg_zzzz__)./BB_inv_std_zzzz__);
CC_inv_nlp_zzzz__ = -z_to_p_twosided_0((CC_inv_zzzz_o__ - CC_inv_avg_zzzz__)./CC_inv_std_zzzz__);
a_snlp_zzzz_ = sign(a_zzzz_o_ - a_avg_zzzz_).*a_nlp_zzzz_;
A_snlp_zzzz__ = sign(A_zzzz_o__ - A_avg_zzzz__).*A_nlp_zzzz__;
BB_inv_snlp_zzzz__ = sign(BB_inv_zzzz_o__ - BB_inv_avg_zzzz__).*BB_inv_nlp_zzzz__;
CC_inv_snlp_zzzz__ = sign(CC_inv_zzzz_o__ - CC_inv_avg_zzzz__).*CC_inv_nlp_zzzz__;
a_sbnlp_zzzz_ = sign(a_zzzz_o_ - a_avg_zzzz_).*max(0,a_nlp_zzzz_-log(n_var));
A_sbnlp_zzzz__ = sign(A_zzzz_o__ - A_avg_zzzz__).*max(0,A_nlp_zzzz__-log(n_var*(n_var-1)));
BB_inv_sbnlp_zzzz__ = sign(BB_inv_zzzz_o__ - BB_inv_avg_zzzz__).*max(0,BB_inv_nlp_zzzz__-log(n_var*(n_var-1)/2));
CC_inv_sbnlp_zzzz__ = sign(CC_inv_zzzz_o__ - CC_inv_avg_zzzz__).*max(0,CC_inv_nlp_zzzz__-log(n_var*(n_var-1)/2));
%%%%%%%%;
end;%if flag_single;
%%%%%%%%

%%%%%%%%;
if flag_double;
%%%%%%%%;
fname_mat_aaaa_o = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix_aaaa);
if ( exist(fname_mat_aaaa_o,'file'));
tmp_aaaa_o_ = load(fname_mat_aaaa_o);
a_aaaa_o_ = tmp_aaaa_o_.a_; n_var = numel(a_aaaa_o_);
A_aaaa_o__ = tmp_aaaa_o_.A__;
BB_inv_aaaa_o__ = tmp_aaaa_o_.BB_inv__;
CC_inv_aaaa_o__ = tmp_aaaa_o_.CC_inv__;
end;%if ( exist(fname_mat_aaaa_s,'file'));
fname_mat_bbbb_o = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix_bbbb);
if ( exist(fname_mat_bbbb_o,'file'));
tmp_bbbb_o_ = load(fname_mat_bbbb_o);
a_bbbb_o_ = tmp_bbbb_o_.a_;
A_bbbb_o__ = tmp_bbbb_o_.A__;
BB_inv_bbbb_o__ = tmp_bbbb_o_.BB_inv__;
CC_inv_bbbb_o__ = tmp_bbbb_o_.CC_inv__;
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
n_found_zzzz = min(numel(tmp_aaaa_s_.L_prm_),numel(tmp_bbbb_s_.L_prm_)) - 1;
n_shuffle = n_found_zzzz;
a_avg_zzzz_ = mean(tmp_aaaa_s_.a_prm__(:,2:end) - tmp_bbbb_s_.a_prm__(:,2:end),2);
a_std_zzzz_ = std(tmp_aaaa_s_.a_prm__(:,2:end) - tmp_bbbb_s_.a_prm__(:,2:end),1,2);
A_avg_zzzz__ = mean(tmp_aaaa_s_.A_prm___(:,:,2:end) - tmp_bbbb_s_.A_prm___(:,:,2:end),3);
A_std_zzzz__ = std(tmp_aaaa_s_.A_prm___(:,:,2:end) - tmp_bbbb_s_.A_prm___(:,:,2:end),1,3);
BB_inv_avg_zzzz__ = mean(tmp_aaaa_s_.BB_inv_prm___(:,:,2:end) - tmp_bbbb_s_.BB_inv_prm___(:,:,2:end),3);
BB_inv_std_zzzz__ = std(tmp_aaaa_s_.BB_inv_prm___(:,:,2:end) - tmp_bbbb_s_.BB_inv_prm___(:,:,2:end),1,3);
CC_inv_avg_zzzz__ = mean(tmp_aaaa_s_.CC_inv_prm___(:,:,2:end) - tmp_bbbb_s_.CC_inv_prm___(:,:,2:end),3);
CC_inv_std_zzzz__ = std(tmp_aaaa_s_.CC_inv_prm___(:,:,2:end) - tmp_bbbb_s_.CC_inv_prm___(:,:,2:end),1,3);
end;%if ( exist(fname_mat_aaaa_s,'file') &  exist(fname_mat_bbbb_s,'file'));
%%%%;
tmp_zzzz_s_ = struct('type','tmp');
tmp_zzzz_s_.a_prm__ = tmp_aaaa_s_.a_prm__ - tmp_bbbb_s_.a_prm__;
tmp_zzzz_s_.A_prm___ = tmp_aaaa_s_.A_prm___ - tmp_bbbb_s_.A_prm___;
tmp_zzzz_s_.BB_inv_prm___ = tmp_aaaa_s_.BB_inv_prm___ - tmp_bbbb_s_.BB_inv_prm___;
tmp_zzzz_s_.CC_inv_prm___ = tmp_aaaa_s_.CC_inv_prm___ - tmp_bbbb_s_.CC_inv_prm___;
tmp_zzzz_s_.L_prm_ = tmp_aaaa_s_.L_prm_ - tmp_bbbb_s_.L_prm_;
%%%%;
a_nlp_zzzz_ = -z_to_p_twosided_0((a_zzzz_o_ - a_avg_zzzz_)./a_std_zzzz_);
A_nlp_zzzz__ = -z_to_p_twosided_0((A_zzzz_o__ - A_avg_zzzz__)./A_std_zzzz__);
BB_inv_nlp_zzzz__ = -z_to_p_twosided_0((BB_inv_zzzz_o__ - BB_inv_avg_zzzz__)./BB_inv_std_zzzz__);
CC_inv_nlp_zzzz__ = -z_to_p_twosided_0((CC_inv_zzzz_o__ - CC_inv_avg_zzzz__)./CC_inv_std_zzzz__);
a_snlp_zzzz_ = sign(a_zzzz_o_ - a_avg_zzzz_).*a_nlp_zzzz_;
A_snlp_zzzz__ = sign(A_zzzz_o__ - A_avg_zzzz__).*A_nlp_zzzz__;
BB_inv_snlp_zzzz__ = sign(BB_inv_zzzz_o__ - BB_inv_avg_zzzz__).*BB_inv_nlp_zzzz__;
CC_inv_snlp_zzzz__ = sign(CC_inv_zzzz_o__ - CC_inv_avg_zzzz__).*CC_inv_nlp_zzzz__;
a_sbnlp_zzzz_ = sign(a_zzzz_o_ - a_avg_zzzz_).*max(0,a_nlp_zzzz_-log(n_var));
A_sbnlp_zzzz__ = sign(A_zzzz_o__ - A_avg_zzzz__).*max(0,A_nlp_zzzz__-log(n_var*(n_var-1)));
BB_inv_sbnlp_zzzz__ = sign(BB_inv_zzzz_o__ - BB_inv_avg_zzzz__).*max(0,BB_inv_nlp_zzzz__-log(n_var*(n_var-1)/2));
CC_inv_sbnlp_zzzz__ = sign(CC_inv_zzzz_o__ - CC_inv_avg_zzzz__).*max(0,CC_inv_nlp_zzzz__-log(n_var*(n_var-1)/2));
%%%%%%%%;
end;%if flag_double;
%%%%%%%%;

%%%%%%%%;
disp(sprintf(' %% n_found_zzzz %d',n_found_zzzz));
%%%%%%%%;

%%%%%%%%;
% Test out holm-bonferroni test. ;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_shuffle_A_sjnlp_FIGA',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
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
,n_var ...
,'A' ...
,n_shuffle ...
,tmp_zzzz_s_.A_prm___ ...
);
figure(1);clf;figbig;fontsize_use = 12; np=0;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1); c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1); c_nlpvt2_lim_ = [-27,+27];
%%%%%%%%;
imagesc(clear_diag_0(A_sqr_sjnlp___(:,:,1+0)),c_nlpvt2_lim_); axis image;
l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ p-value (holm)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
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

fname_fig = sprintf('%s/dolphin_%s_shuffle_B_sjnlp_FIGB',dir_jpg,fname_infix_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
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
,n_var ...
,'A' ...
,n_shuffle ...
,tmp_zzzz_s_.BB_inv_prm___ ...
);
figure(1);clf;figbig;fontsize_use = 12; np=0;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1); c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1); c_nlpvt2_lim_ = [-27,+27];
%%%%%%%%;
imagesc(clear_diag_0(BB_inv_sqr_sjnlp___(:,:,1+0)),c_nlpvt2_lim_); axis image;
l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ p-value (holm)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
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
subplot_(1,1) = subplot(1,2,np);np=np+1; imagesc(clear_diag_0(A_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ p-value (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
subplot_(1,2) = subplot(1,2,np);np=np+1; imagesc(clear_diag_0(BB_inv_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ p-value (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
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
subplot_(1,1) = subplot(1,1,np);np=np+1; imagesc(clear_diag_0(A_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ z-score (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
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
subplot_(1,1) = subplot(1,1,np);np=np+1; imagesc(clear_diag_0(BB_inv_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ z-score (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
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
subplot_(1,1) = subplot(1,1,np);np=np+1; imagesc(clear_diag_0(CC_inv_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $CC^{T}$ z-score (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
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
tmp_x_ = zeros(4,n_var); tmp_y_ = zeros(4,n_var); tmp_c_ = zeros(1,n_var,3);
for nvar=0:n_var-1; 
tmp_y = a_avg_zzzz_(1+nvar); tmp_x_(:,1+nvar) = -0.5+1+nvar+[0;0;1;1]; tmp_y_(:,1+nvar) = [0;tmp_y;tmp_y;0]; 
tmp_nc = max(0,min(n_c_beach-1,floor(n_c_beach*(tmp_y-min(c_beach_lim_))/diff(c_beach_lim_))));
tmp_c_(1,1+nvar,:) = c_beach__(1+tmp_nc,:);
end;%for nvar=0:n_var-1; 
patch(tmp_x_,tmp_y_,tmp_c_,'EdgeColor','k');
xlim([0,n_var+1]);ylim(0.125*c_beach_lim_);
%set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source/target'); ylabel('magnitude'); title('$a$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
subplot_(1,2) = subplot(2,3,np);np=np+1; 
tmp_x_ = zeros(4,n_var); tmp_y_ = zeros(4,n_var); tmp_c_ = zeros(1,n_var,3);
for nvar=0:n_var-1; 
tmp_y = a_std_zzzz_(1+nvar); tmp_x_(:,1+nvar) = -0.5+1+nvar+[0;0;1;1]; tmp_y_(:,1+nvar) = [0;tmp_y;tmp_y;0]; 
tmp_nc = max(0,min(n_c_beach-1,floor(n_c_beach*(tmp_y-min(c_beach_lim_))/diff(c_beach_lim_))));
tmp_c_(1,1+nvar,:) = c_beach__(1+tmp_nc,:);
end;%for nvar=0:n_var-1; 
patch(tmp_x_,tmp_y_,tmp_c_,'EdgeColor','k');
xlim([0,n_var+1]);ylim(0.125*c_beach_lim_);
%set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source/target'); ylabel('magnitude'); title('$a$ std','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
subplot_(1,3) = subplot(2,3,np);np=np+1;
tmp_x_ = zeros(4,n_var); tmp_y_ = zeros(4,n_var); tmp_c_ = zeros(1,n_var,3);
for nvar=0:n_var-1; 
tmp_y = a_snlp_zzzz_(1+nvar); tmp_x_(:,1+nvar) = -0.5+1+nvar+[0;0;1;1]; tmp_y_(:,1+nvar) = [0;tmp_y;tmp_y;0]; 
tmp_nc = max(0,min(n_c_nlpvt2-1,floor(n_c_nlpvt2*(tmp_y-min(c_nlpvt2_lim_))/diff(c_nlpvt2_lim_))));
tmp_c_(1,1+nvar,:) = c_nlpvt2__(1+tmp_nc,:);
end;%for nvar=0:n_var-1; 
patch(tmp_x_,tmp_y_,tmp_c_,'EdgeColor','k');
xlim([0,n_var+1]);ylim(c_nlpvt2_lim_);
%set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source/target'); ylabel('magnitude'); title('signed $a$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
%%%%%%%%;
subplot_(2,1) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(A_avg_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$A$ avg','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,2) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(A_std_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$A$ std','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,3) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(A_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
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
subplot_(1,1) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(BB_inv_avg_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$BB^{T}$ avg','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(1,2) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(BB_inv_std_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$BB^{T}$ std','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(1,3) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(BB_inv_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(2,1) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(CC_inv_avg_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$CC^{T}$ avg','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,2) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(CC_inv_std_zzzz__),c_beach_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('$CC^{T}$ std','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,3) = subplot(2,3,np);np=np+1; imagesc(clear_diag_0(CC_inv_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
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
subplot_(1,1) = subplot(2,2,np);np=np+1; imagesc(clear_diag_0(A_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(1,2) = subplot(2,2,np);np=np+1; imagesc(clear_diag_0(BB_inv_snlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $BB^{T}$ z-score','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(2,1) = subplot(2,2,np);np=np+1; imagesc(clear_diag_0(A_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed $A$ z-score (bonferroni)','Interpreter','latex');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(2,2) = subplot(2,2,np);np=np+1; imagesc(clear_diag_0(BB_inv_sbnlp_zzzz__),c_nlpvt2_lim_); axis image; l = slash_diag_0(n_var); set(l,'Color',0.85*[1,1,1],'LineWidth',0.5);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% check distribution of shuffle coefficients. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_mat_zzzz_s = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix_zzzz);
if ( exist(fname_mat_zzzz_s,'file'));
n_found_zzzz = size(tmp_zzzz_s_.L_prm_,1) - 1;
a_zzzz__ = tmp_zzzz_s_.a_prm__(:,2:end);
A_zzzz__ = reshape(tmp_zzzz_s_.A_prm___(:,:,2:end),[n_var^2,n_found_zzzz]);
BB_inv_zzzz__ = reshape(tmp_zzzz_s_.BB_inv_prm___(:,:,2:end),[n_var^2,n_found_zzzz]);
CC_inv_zzzz__ = zeros(n_var,n_found_zzzz);
for nfound=0:n_found_zzzz-1;
CC_inv_zzzz__(:,1+nfound)= diag(tmp_zzzz_s_.CC_inv_prm___(:,:,1+1+nfound));
end;%for nfound=0:n_found_zzzz-1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_exist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (verbose); disp(sprintf(' %% [finished dolphin_d00_collect_and_plot_3]')); end;
