function ...
[ ...
 parameter ...
] = ...
label_tree_compilation_plot_1( ...
 parameter ...
,dir_base ...
,str_infix ...
,str_sgtitle ...
,local_fname_output_label ...
,local_fname_nlpbra_label ...
,local_fname_nlpnex_label ...
,X_ ...
,rank_estimate_X ...
,label_A_ ...
,fname_louvain_mat ...
,fname_umap_mat ...
);

verbose=1;
if (verbose); disp(sprintf(' %% [entering label_tree_compilation_plot_1]')); end;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
if (~isfield(parameter,'date_diff_threshold')); parameter.date_diff_threshold = 0.5; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_force_create_mat')); parameter.flag_force_create_mat = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_force_create_tmp')); parameter.flag_force_create_tmp = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_replot')); parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
date_diff_threshold = parameter.date_diff_threshold;
flag_force_create_mat = parameter.flag_force_create_mat;
flag_force_create_tmp = parameter.flag_force_create_tmp;
flag_replot = parameter.flag_replot;

dir_jpg = sprintf('%s/dir_jpg',dir_base);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

local_output_label__ = [];
local_nlpbra_label__ = [];
local_nlpnex_label__ = [];
flag_exist = exist(local_fname_output_label,'file') & exist(local_fname_nlpbra_label,'file') & exist(local_fname_nlpnex_label,'file');
if (~flag_exist);
disp(sprintf(' %% Warning, ~flag_exist in label_tree_compilation_plot_1'));
end;%if (~flag_exist);
if (flag_exist);
local_fname_pre = sprintf('%s/%s_postprocess',dir_base,str_infix);
[local_flag_skip,local_fname_mat] = open_fname_tmp(local_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~local_flag_skip;
fp = fopen(local_fname_output_label); local_output_label__ = textscan(fp,'%s','Delimiter','\n'); local_output_label__ = local_output_label__{1}; fclose(fp);
fp = fopen(local_fname_nlpbra_label); local_nlpbra_label__ = textscan(fp,'%s','Delimiter','\n'); local_nlpbra_label__ = local_nlpbra_label__{1}; fclose(fp);
fp = fopen(local_fname_nlpnex_label); local_nlpnex_label__ = textscan(fp,'%s','Delimiter','\n'); local_nlpnex_label__ = local_nlpnex_label__{1}; fclose(fp);
[lpv] = label_to_label_enrichment_quad_4(label_A_,label_str_to_enum_0(local_output_label__));
%%%%%%%%;
try;
label_tree ...
= ...
label_to_tree_2( ...
 0 ...
,[] ...
,[] ...
,[] ...
,local_output_label__ ...
,local_nlpbra_label__ ...
,local_nlpnex_label__ ...
);
disp(sprintf(' %% tree compatible with label_to_tree_2, not trying label_to_tree_1'));
catch;
disp(sprintf(' %% tree not compatible with label_to_tree_2, trying label_to_tree_1'));
label_tree ...
= ...
label_to_tree_1( ...
 0 ...
,[] ...
,[] ...
,[] ...
,local_output_label__ ...
,local_nlpbra_label__ ...
,local_nlpnex_label__ ...
);
end;%try;
%%%%%%%%%%%%%%%%;
[U_,S_,V_] = svds(X_,rank_estimate_X);
US_X_sg__ = U_*S_; clear U_ S_ V_;
[ ...
 parameter ...
,nlpnex_threshold_ ...
,psplit_threshold_ ...
,~ ...
,lpv_nA__ ...
,n_cluster_n_ ...
,label_rearrange_n_iteration ...
,lpv_rearrange_nAi___ ...
,n_cluster_rearrange_ni__ ...
] = ...
label_to_tree_node_by_node_enrichment_2( ...
 parameter ...
 ,label_tree ...
 ,{label_A_} ...
 ,US_X_sg__ ...
);
%%%%%%%%;
save(local_fname_mat ...
     ,'nlpnex_threshold_' ...
     ,'psplit_threshold_' ...
     ,'lpv_nA__' ...
     ,'n_cluster_n_' ...
     ,'label_rearrange_n_iteration' ...
     ,'lpv_rearrange_nAi___' ...
     ,'n_cluster_rearrange_ni__' ...
     );
close_fname_tmp(local_fname_pre);
end;%if ~local_flag_skip;
%%%%%%%%;
n_cluster_louvain00_default = 0;
lpv_louvain00_default = 0;
if ( exist(fname_louvain_mat,'file'));
tmp_ = load(sprintf('%s',fname_louvain_mat));
tmp_label_B_ = tmp_.label_B__{1};
n_cluster_louvain00_default = numel(unique(tmp_label_B_));
lpv_louvain00_default = label_to_label_enrichment_quad_4(label_A_,tmp_label_B_); %<-- should equal lpv. ;
clear tmp_ tmp_label_B_;
end;%if ( exist(fname_louvain_mat,'file'));
%%%%%%%%;
n_cluster_umap00_default = 0;
lpv_umap00_default = 0;
if ( exist(fname_umap_mat,'file'));
tmp_ = load(sprintf('%s',fname_umap_mat));
tmp_label_B_ = tmp_.label_B__{1};
n_cluster_umap00_default = numel(unique(tmp_label_B_));
lpv_umap00_default = label_to_label_enrichment_quad_4(label_A_,tmp_label_B_); %<-- should equal lpv. ;
clear tmp_ tmp_label_B_;
end;%if ( exist(fname_umap_mat,'file'));
%%%%%%%%;
if exist(local_fname_mat,'file');
load(local_fname_mat);
n_cluster_rearrange_n0_ = n_cluster_rearrange_ni__(:,1+0);
n_cluster_rearrange_ne_ = n_cluster_rearrange_ni__(:,end);
lpv_rearrange_nA0__ = lpv_rearrange_nAi___(:,:,1+0);
lpv_rearrange_nAe__ = lpv_rearrange_nAi___(:,:,end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_fig = sprintf('%s/local_lpv_nlpnex_threshold_%s_FIGA',dir_jpg,str_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;
%%%%%%%%;
subplot(1,1,1); hold on;
plot(nlpnex_threshold_,-lpv_nA__(:,1+0),'m.-');
plot(nlpnex_threshold_,-lpv_rearrange_nA0__(:,1+0),'mo-');
plot(nlpnex_threshold_,-lpv_rearrange_nAe__(:,1+0),'ms-');
plot(nlpnex_threshold_,-lpv_louvain00_default*ones(numel(nlpnex_threshold_),1),'rx-');
plot(nlpnex_threshold_,-lpv_umap00_default*ones(numel(nlpnex_threshold_),1),'k^-');
xlabel('-log(p)'); ylabel('-log(P)');
title('all','Interpreter','none');
%%%%%%%%;
sgtitle(sprintf('%s',str_sgtitle),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_fig = sprintf('%s/local_lpv_nlpnex_threshold_%s_FIGB',dir_jpg,str_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;
%%%%%%%%;
subplot(1,1,1); hold on;
xlabel('-log(p)');
%%%%;
yyaxis right;
hold on;
stairs(nlpnex_threshold_,n_cluster_rearrange_n0_,'m-','LineWidth',2);
%stairs(nlpnex_threshold_,n_cluster_rearrange_ne_,'ms-','LineWidth',2);
hold off;
ylabel('n_cluster','Interpreter','none');
%%%%;
yyaxis left;
hold on;
edge_ = [nlpnex_threshold_;3];
lim_ = [0,max(-lpv_rearrange_nA0__(:,1+0))];
p = stairfall_1(edge_,-lpv_rearrange_nA0__(:,1+0),lim_,0.85*[1,1,1]); set(p,'EdgeColor','none');
%p = stairfall_1(edge_,-lpv_rearrange_nAe__(:,1+0),lim_,0.85*[1,1,1]); set(p,'EdgeColor',0.65*[1,1,1],'FaceColor','none');
plot(edge_,-lpv_louvain00_default*ones(numel(edge_),1),'r--','LineWidth',2);
plot(edge_,-lpv_umap00_default*ones(numel(edge_),1),'k--','LineWidth',2);
hold off;
ylabel('-log(P)');
%%%%%%%%;
sgtitle(sprintf('%s',str_sgtitle),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_fig = sprintf('%s/local_lpv_nlpnex_threshold_%s_FIGC',dir_jpg,str_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;
%%%%%%%%;
subplot(1,1,1); hold on;
xlabel('n_cluster','Interpreter','none');
%%%%;
yyaxis right;
hold on;
stairs(n_cluster_rearrange_n0_,-lpv_rearrange_nA0__(:,1+0),'m-','LineWidth',2);
plot(n_cluster_louvain00_default,-lpv_louvain00_default,'rx','LineWidth',2,'MarkerSize',16);
plot(n_cluster_louvain00_default,-lpv_louvain00_default,'ro','LineWidth',2,'MarkerSize',16);
plot(n_cluster_umap00_default,-lpv_umap00_default,'kx','LineWidth',2,'MarkerSize',16);
plot(n_cluster_umap00_default,-lpv_umap00_default,'ko','LineWidth',2,'MarkerSize',16);
hold off;
ylabel('-log(P)','Interpreter','none');
%%%%;
yyaxis left;
hold on;
edge_ = [n_cluster_rearrange_n0_;max(n_cluster_rearrange_n0_)+0];
lim_ = [0,max(nlpnex_threshold_)];
p = stairfall_1(edge_,nlpnex_threshold_,lim_,0.85*[1,1,1]); set(p,'EdgeColor','none');
hold off;
ylabel('-log(p)');
%%%%%%%%;
sgtitle(sprintf('%s',str_sgtitle),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_fig = sprintf('%s/local_lpv_nlpnex_threshold_%s_FIGD',dir_jpg,str_infix);
if (1 | flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
if isempty(local_output_label__); fp = fopen(local_fname_output_label); local_output_label__ = textscan(fp,'%s','Delimiter','\n'); local_output_label__ = local_output_label__{1}; fclose(fp); end;%if isempty(local_output_label__);
if isempty(local_nlpbra_label__); fp = fopen(local_fname_nlpbra_label); local_nlpbra_label__ = textscan(fp,'%s','Delimiter','\n'); local_nlpbra_label__ = local_nlpbra_label__{1}; fclose(fp); end;%if isempty(local_nlpbra_label__);
if isempty(local_nlpnex_label__);fp = fopen(local_fname_nlpnex_label); local_nlpnex_label__ = textscan(fp,'%s','Delimiter','\n'); local_nlpnex_label__ = local_nlpnex_label__{1}; fclose(fp); end;%if isempty(local_nlpnex_label__);
figure(1);clf;figbig;
nlplim_ = [0,27];
subplot(1,6,[6:6]);
hold on;
edge_ = [nlpnex_threshold_;3];
%stairs([0;-lpv_rearrange_nA0__(:,1+0)],edge_,'m-','LineWidth',2);
tmp_x_ = [0;-lpv_rearrange_nA0__(:,1+0)]; tmp_x_ = repmat(transpose(tmp_x_),[2,1]); tmp_x_ = tmp_x_(:);
tmp_y_ = edge_; tmp_y_ = repmat(transpose(tmp_y_),[2,1]); tmp_y_ = tmp_y_(:); tmp_y_ = circshift(tmp_y_,+1);
p = patch(tmp_x_,tmp_y_,0.85*[1,1,1]); set(p,'EdgeColor','none');
plot(-lpv_louvain00_default*ones(numel(edge_),1),edge_,'r--','LineWidth',2);
plot(-lpv_umap00_default*ones(numel(edge_),1),edge_,'k--','LineWidth',2);
text(-lpv_louvain00_default,min(edge_)-1,num2str(n_cluster_louvain00_default),'Color','r');
text(-lpv_umap00_default,min(edge_)-2,num2str(n_cluster_umap00_default),'Color','k');
hold off;
xlabel('-log(P)'); ylim([0,max(nlplim_)+1]);
set(gca,'YTick',0:max(nlplim_)+1,'YTickLabel',{});
grid on;
subplot(1,6,[1:5]);
tmp_parameter = struct('type','parameter');
tmp_parameter.flag_ylabel = 0;
label_plot_recursive_nlpscale_3( ...
 tmp_parameter ...
,local_output_label__ ...
,local_nlpbra_label__ ...
,local_nlpnex_label__ ...
,[] ...
,nlplim_ ...
,[] ...
);
%%%%%%%%;
sgtitle(sprintf('%s',str_sgtitle),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if exist(local_fname_mat,'file');
%%%%%%%%;
end;%if (flag_exist);

if (verbose); disp(sprintf(' %% [finished label_tree_compilation_plot_1]')); end;
