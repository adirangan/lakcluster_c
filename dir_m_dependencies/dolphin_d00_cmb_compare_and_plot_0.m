function ...
[ ...
 parameter ...
 ,n_found_yyyy_vs_zzzz ...
 ,a_yyyy_vs_zzzz_o_ ...
 ,A_yyyy_vs_zzzz_o__ ...
 ,BB_yyyy_vs_zzzz_o__ ...
 ,CC_yyyy_vs_zzzz_o__ ...
 ,a_sbnlp_yyyy_vs_zzzz_ ...
 ,A_sbnlp_yyyy_vs_zzzz__ ...
 ,BB_sbnlp_yyyy_vs_zzzz__ ...
 ,CC_sbnlp_yyyy_vs_zzzz__ ...
 ,a_snlp_yyyy_vs_zzzz_ ...
 ,A_snlp_yyyy_vs_zzzz__ ...
 ,BB_snlp_yyyy_vs_zzzz__ ...
 ,CC_snlp_yyyy_vs_zzzz__ ...
 ,a_nlp_yyyy_vs_zzzz_ ...
 ,A_nlp_yyyy_vs_zzzz__ ...
 ,BB_nlp_yyyy_vs_zzzz__ ...
 ,CC_nlp_yyyy_vs_zzzz__ ...
 ,a_avg_yyyy_vs_zzzz_ ...
 ,A_avg_yyyy_vs_zzzz__ ...
 ,BB_avg_yyyy_vs_zzzz__ ...
 ,CC_avg_yyyy_vs_zzzz__ ...
 ,a_std_yyyy_vs_zzzz_ ...
 ,A_std_yyyy_vs_zzzz__ ...
 ,BB_std_yyyy_vs_zzzz__ ...
 ,CC_std_yyyy_vs_zzzz__ ...
] = ...
dolphin_d00_cmb_compare_and_plot_0( ...
 parameter ...
 ,n_var ...
 ,string_dat_name_ ...
 ,dir_mat ...
 ,dir_shuffle_mat ...
 ,dir_jpg ...
 ,fname_infix_yyyy ...
 ,fname_infix_zzzz ...
 ,fname_infix_yyyy_vs_zzzz ...
);

verbose=1;
if (verbose); disp(sprintf(' %% [entering dolphin_d00_cmb_compare_and_plot_0]')); end;

if (isempty(parameter)); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'flag_replot')); parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
flag_replot = parameter.flag_replot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% yyyy_vs_zzzz ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_mat_yyyy_o = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix_yyyy);
if ( exist(fname_mat_yyyy_o,'file'));
tmp_ = load(fname_mat_yyyy_o);
a_yyyy_o_ = tmp_.a_cmb_;
A_yyyy_o__ = tmp_.A_cmb__;
BB_yyyy_o__ = tmp_.BB_cmb__;
CC_yyyy_o__ = tmp_.CC_cmb__;
end;%if ( exist(fname_mat_yyyy_s,'file'));
fname_mat_zzzz_o = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix_zzzz);
if ( exist(fname_mat_zzzz_o,'file'));
tmp_ = load(fname_mat_zzzz_o);
a_zzzz_o_ = tmp_.a_cmb_;
A_zzzz_o__ = tmp_.A_cmb__;
BB_zzzz_o__ = tmp_.BB_cmb__;
CC_zzzz_o__ = tmp_.CC_cmb__;
end;%if ( exist(fname_mat_zzzz_s,'file'));
a_yyyy_vs_zzzz_o_ = a_yyyy_o_ - a_zzzz_o_;
A_yyyy_vs_zzzz_o__ = A_yyyy_o__ - A_zzzz_o__;
BB_yyyy_vs_zzzz_o__ = BB_yyyy_o__ - BB_zzzz_o__;
CC_yyyy_vs_zzzz_o__ = CC_yyyy_o__ - CC_zzzz_o__;
%%%%%%%%;
a_avg_yyyy_vs_zzzz_ = zeros(n_var,1); a_std_yyyy_vs_zzzz_ = zeros(n_var,1);
A_avg_yyyy_vs_zzzz__ = zeros(n_var,n_var); A_std_yyyy_vs_zzzz__ = zeros(n_var,n_var);
BB_avg_yyyy_vs_zzzz__ = zeros(n_var,n_var); BB_std_yyyy_vs_zzzz__ = zeros(n_var,n_var);
CC_avg_yyyy_vs_zzzz__ = zeros(n_var,n_var); CC_std_yyyy_vs_zzzz__ = zeros(n_var,n_var);
%%%%;
fname_mat_yyyy_s = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix_yyyy);
fname_mat_zzzz_s = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix_zzzz);
if ( exist(fname_mat_yyyy_s,'file') &  exist(fname_mat_zzzz_s,'file'));
tmp_yyyy_ = load(fname_mat_yyyy_s);
tmp_zzzz_ = load(fname_mat_zzzz_s);
n_found_yyyy_vs_zzzz = min(numel(tmp_yyyy_.L_cmb_prm_),numel(tmp_zzzz_.L_cmb_prm_)) - 1;
a_avg_yyyy_vs_zzzz_ = mean(tmp_yyyy_.a_cmb_prm__(:,2:end) - tmp_zzzz_.a_cmb_prm__(:,2:end),2);
a_std_yyyy_vs_zzzz_ = std(tmp_yyyy_.a_cmb_prm__(:,2:end) - tmp_zzzz_.a_cmb_prm__(:,2:end),1,2);
A_avg_yyyy_vs_zzzz__ = mean(tmp_yyyy_.A_cmb_prm___(:,:,2:end) - tmp_zzzz_.A_cmb_prm___(:,:,2:end),3);
A_std_yyyy_vs_zzzz__ = std(tmp_yyyy_.A_cmb_prm___(:,:,2:end) - tmp_zzzz_.A_cmb_prm___(:,:,2:end),1,3);
BB_avg_yyyy_vs_zzzz__ = mean(tmp_yyyy_.BB_cmb_prm___(:,:,2:end) - tmp_zzzz_.BB_cmb_prm___(:,:,2:end),3);
BB_std_yyyy_vs_zzzz__ = std(tmp_yyyy_.BB_cmb_prm___(:,:,2:end) - tmp_zzzz_.BB_cmb_prm___(:,:,2:end),1,3);
CC_avg_yyyy_vs_zzzz__ = mean(tmp_yyyy_.CC_cmb_prm___(:,:,2:end) - tmp_zzzz_.CC_cmb_prm___(:,:,2:end),3);
CC_std_yyyy_vs_zzzz__ = std(tmp_yyyy_.CC_cmb_prm___(:,:,2:end) - tmp_zzzz_.CC_cmb_prm___(:,:,2:end),1,3);
end;%if ( exist(fname_mat_yyyy_s,'file') &  exist(fname_mat_zzzz_s,'file'));
%%%%;
a_nlp_yyyy_vs_zzzz_ = -z_to_p_twosided_0((a_yyyy_vs_zzzz_o_ - a_avg_yyyy_vs_zzzz_)./a_std_yyyy_vs_zzzz_);
A_nlp_yyyy_vs_zzzz__ = -z_to_p_twosided_0((A_yyyy_vs_zzzz_o__ - A_avg_yyyy_vs_zzzz__)./A_std_yyyy_vs_zzzz__);
BB_nlp_yyyy_vs_zzzz__ = -z_to_p_twosided_0((BB_yyyy_vs_zzzz_o__ - BB_avg_yyyy_vs_zzzz__)./BB_std_yyyy_vs_zzzz__);
CC_nlp_yyyy_vs_zzzz__ = -z_to_p_twosided_0((CC_yyyy_vs_zzzz_o__ - CC_avg_yyyy_vs_zzzz__)./CC_std_yyyy_vs_zzzz__);
a_snlp_yyyy_vs_zzzz_ = sign(a_yyyy_vs_zzzz_o_ - a_avg_yyyy_vs_zzzz_).*a_nlp_yyyy_vs_zzzz_;
A_snlp_yyyy_vs_zzzz__ = sign(A_yyyy_vs_zzzz_o__ - A_avg_yyyy_vs_zzzz__).*A_nlp_yyyy_vs_zzzz__;
BB_snlp_yyyy_vs_zzzz__ = sign(BB_yyyy_vs_zzzz_o__ - BB_avg_yyyy_vs_zzzz__).*BB_nlp_yyyy_vs_zzzz__;
CC_snlp_yyyy_vs_zzzz__ = sign(CC_yyyy_vs_zzzz_o__ - CC_avg_yyyy_vs_zzzz__).*CC_nlp_yyyy_vs_zzzz__;
a_sbnlp_yyyy_vs_zzzz_ = sign(a_yyyy_vs_zzzz_o_ - a_avg_yyyy_vs_zzzz_).*max(0,a_nlp_yyyy_vs_zzzz_-log(n_var));
A_sbnlp_yyyy_vs_zzzz__ = sign(A_yyyy_vs_zzzz_o__ - A_avg_yyyy_vs_zzzz__).*max(0,A_nlp_yyyy_vs_zzzz__-log(n_var*(n_var-1)));
BB_sbnlp_yyyy_vs_zzzz__ = sign(BB_yyyy_vs_zzzz_o__ - BB_avg_yyyy_vs_zzzz__).*max(0,BB_nlp_yyyy_vs_zzzz__-log(n_var*(n_var-1)/2));
CC_sbnlp_yyyy_vs_zzzz__ = sign(CC_yyyy_vs_zzzz_o__ - CC_avg_yyyy_vs_zzzz__).*max(0,CC_nlp_yyyy_vs_zzzz__-log(n_var*(n_var-1)/2));
%%%%%%%%;
disp(sprintf(' %% n_found_yyyy_vs_zzzz %d',n_found_yyyy_vs_zzzz'));
%%%%%%%%;

fname_fig = sprintf('%s/dolphin_%s_shuffle_AB2_FIGAB2',dir_jpg,fname_infix_yyyy_vs_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 12;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(1,2,np);np=np+1; imagesc(A_sbnlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed A_nlp_yyyy_vs_zzzz__ (bonferroni)','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
subplot_(1,2) = subplot(1,2,np);np=np+1; imagesc(BB_sbnlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed BB_nlp_yyyy_vs_zzzz__ (bonferroni)','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_yyyy_vs_zzzz),'Interpreter','none');
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

fname_fig = sprintf('%s/dolphin_%s_shuffle_A_FIGA',dir_jpg,fname_infix_yyyy_vs_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;set(gcf,'Position',1+[0,0,1024,1024]);fontsize_use = 12;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(1,1,np);np=np+1; imagesc(A_sbnlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed A_nlp_yyyy_vs_zzzz__ (bonferroni)','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_yyyy_vs_zzzz),'Interpreter','none');
colormap(subplot_(1,1),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/dolphin_%s_shuffle_B_FIGB',dir_jpg,fname_infix_yyyy_vs_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;set(gcf,'Position',1+[0,0,1024,1024]);fontsize_use = 12;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(1,1,np);np=np+1; imagesc(BB_sbnlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed BB_nlp_yyyy_vs_zzzz__ (bonferroni)','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_yyyy_vs_zzzz),'Interpreter','none');
colormap(subplot_(1,1),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/dolphin_%s_shuffle_C_FIGC',dir_jpg,fname_infix_yyyy_vs_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;set(gcf,'Position',1+[0,0,1024,1024]);fontsize_use = 12;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(1,1,np);np=np+1; imagesc(CC_sbnlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed CC_nlp_yyyy_vs_zzzz__ (bonferroni)','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
hold on;
l_x_ = +0.5 + (n_var)*[zeros(1,n_var);ones(1,n_var)];
l_y_ = +0.5 + repmat([0:n_var-1],[2,1]);
line(l_x_,l_y_,'Color',0.85*[1,1,1],'LineWidth',0.5);
line(l_y_,l_x_,'Color',0.85*[1,1,1],'LineWidth',0.5);
hold off;
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_yyyy_vs_zzzz),'Interpreter','none');
colormap(subplot_(1,1),c_nlpvt2__);
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/dolphin_%s_shuffle_aA_FIGaA',dir_jpg,fname_infix_yyyy_vs_zzzz);
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
tmp_y = a_avg_yyyy_vs_zzzz_(1+nvar); tmp_x_(:,1+nvar) = -0.5+1+nvar+[0;0;1;1]; tmp_y_(:,1+nvar) = [0;tmp_y;tmp_y;0]; 
tmp_nc = max(0,min(n_c_beach-1,floor(n_c_beach*(tmp_y-min(c_beach_lim_))/diff(c_beach_lim_))));
tmp_c_(1,1+nvar,:) = c_beach__(1+tmp_nc,:);
end;%for nvar=0:n_var-1; 
patch(tmp_x_,tmp_y_,tmp_c_,'EdgeColor','k');
xlim([0,n_var+1]);ylim(0.125*c_beach_lim_);
%set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source/target'); ylabel('magnitude'); title('a_avg_yyyy_vs_zzzz_','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
subplot_(1,2) = subplot(2,3,np);np=np+1; 
tmp_x_ = zeros(4,n_var); tmp_y_ = zeros(4,n_var); tmp_c_ = zeros(1,n_var,3);
for nvar=0:n_var-1; 
tmp_y = a_std_yyyy_vs_zzzz_(1+nvar); tmp_x_(:,1+nvar) = -0.5+1+nvar+[0;0;1;1]; tmp_y_(:,1+nvar) = [0;tmp_y;tmp_y;0]; 
tmp_nc = max(0,min(n_c_beach-1,floor(n_c_beach*(tmp_y-min(c_beach_lim_))/diff(c_beach_lim_))));
tmp_c_(1,1+nvar,:) = c_beach__(1+tmp_nc,:);
end;%for nvar=0:n_var-1; 
patch(tmp_x_,tmp_y_,tmp_c_,'EdgeColor','k');
xlim([0,n_var+1]);ylim(0.125*c_beach_lim_);
%set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source/target'); ylabel('magnitude'); title('a_std_yyyy_vs_zzzz_','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
subplot_(1,3) = subplot(2,3,np);np=np+1;
tmp_x_ = zeros(4,n_var); tmp_y_ = zeros(4,n_var); tmp_c_ = zeros(1,n_var,3);
for nvar=0:n_var-1; 
tmp_y = a_snlp_yyyy_vs_zzzz_(1+nvar); tmp_x_(:,1+nvar) = -0.5+1+nvar+[0;0;1;1]; tmp_y_(:,1+nvar) = [0;tmp_y;tmp_y;0]; 
tmp_nc = max(0,min(n_c_nlpvt2-1,floor(n_c_nlpvt2*(tmp_y-min(c_nlpvt2_lim_))/diff(c_nlpvt2_lim_))));
tmp_c_(1,1+nvar,:) = c_nlpvt2__(1+tmp_nc,:);
end;%for nvar=0:n_var-1; 
patch(tmp_x_,tmp_y_,tmp_c_,'EdgeColor','k');
xlim([0,n_var+1]);ylim(c_nlpvt2_lim_);
%set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source/target'); ylabel('magnitude'); title('signed a_nlp_yyyy_vs_zzzz_','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
%%%%%%%%;
subplot_(2,1) = subplot(2,3,np);np=np+1; imagesc(A_avg_yyyy_vs_zzzz__,c_beach_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('A_avg_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,2) = subplot(2,3,np);np=np+1; imagesc(A_std_yyyy_vs_zzzz__,c_beach_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('A_std_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,3) = subplot(2,3,np);np=np+1; imagesc(A_snlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed A_nlp_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_yyyy_vs_zzzz),'Interpreter','none');
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

fname_fig = sprintf('%s/dolphin_%s_shuffle_BC_FIGBC',dir_jpg,fname_infix_yyyy_vs_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 8;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(2,3,np);np=np+1; imagesc(BB_avg_yyyy_vs_zzzz__,c_beach_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('BB_avg_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(1,2) = subplot(2,3,np);np=np+1; imagesc(BB_std_yyyy_vs_zzzz__,c_beach_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('BB_std_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(1,3) = subplot(2,3,np);np=np+1; imagesc(BB_snlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed BB_nlp_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(2,1) = subplot(2,3,np);np=np+1; imagesc(CC_avg_yyyy_vs_zzzz__,c_beach_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('CC_avg_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,2) = subplot(2,3,np);np=np+1; imagesc(CC_std_yyyy_vs_zzzz__,c_beach_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('CC_std_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
subplot_(2,3) = subplot(2,3,np);np=np+1; imagesc(CC_snlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed CC_nlp_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_yyyy_vs_zzzz),'Interpreter','none');
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

fname_fig = sprintf('%s/dolphin_%s_shuffle_AB_FIGAB',dir_jpg,fname_infix_yyyy_vs_zzzz);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;fontsize_use = 8;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
c_beach_lim_ = [-1,+1];
c_nlpvt2__ = colormap_nlpvt_twosided(); n_c_nlpvt2 = size(c_nlpvt2__,1);
c_nlpvt2_lim_ = [-27,+27];
np = 1;
%%%%%%%%;
subplot_(1,1) = subplot(2,2,np);np=np+1; imagesc(A_snlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed A_nlp_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(1,2) = subplot(2,2,np);np=np+1; imagesc(BB_snlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed BB_nlp_yyyy_vs_zzzz__','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(2,1) = subplot(2,2,np);np=np+1; imagesc(A_sbnlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed A_nlp_yyyy_vs_zzzz__ (bonferroni)','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
subplot_(2,2) = subplot(2,2,np);np=np+1; imagesc(BB_sbnlp_yyyy_vs_zzzz__,c_nlpvt2_lim_); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
xlabel('source'); ylabel('target'); title('signed BB_nlp_yyyy_vs_zzzz__ (bonferroni)','Interpreter','none');
set(gca,'TickLength',[0,0]); set(gca,'fontsize',fontsize_use);
colorbar();
%%%%%%%%;
sgtitle(sprintf('%s',fname_infix_yyyy_vs_zzzz),'Interpreter','none');
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

if (verbose); disp(sprintf(' %% [finished dolphin_d00_cmb_compare_and_plot_0]')); end;
