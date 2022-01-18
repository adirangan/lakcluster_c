function ...
label_tree_SVD_discat_0( ...
 depth ...
,depth_upb ...
,label_tree ...
,alpha_over_beta ...
,X_n__ ...
,label_X_ ...
,dir_jpg ...
,str_infix ...
,str_sgtitle ...
,flag_replot ...
);

if (depth<=depth_upb);
item_index_X_ = [];
if (isfield(label_tree,'item_index_')); item_index_X_ = label_tree.item_index_; end;
n_item_X = numel(item_index_X_);
if (n_item_X>1);
if ( ~isempty(label_tree.childA) & ~isempty(label_tree.childB) );
item_index_A_ = label_tree.childA.item_index_ ;
item_index_B_ = label_tree.childB.item_index_ ;
[item_index_AX_,ij_nX_from_nAX_,ij_nA_from_nAX_] = intersect(item_index_X_,item_index_A_,'stable');
[item_index_BX_,ij_nX_from_nBX_,ij_nB_from_nBX_] = intersect(item_index_X_,item_index_B_,'stable');
assert( numel(item_index_AX_) == numel(item_index_A_) );
assert( numel(item_index_BX_) == numel(item_index_B_) );
assert( numel(item_index_A_) + numel(item_index_B_) == n_item_X );
label_Y_enum_ = zeros(n_item_X,1);
label_Y_enum_(ij_nX_from_nAX_) = 1;
label_Y_enum_(ij_nX_from_nBX_) = 2;
label_X_enum_ = label_num_to_enum_0(label_X_(1+item_index_X_));
%%%%;
[~,USL_n_,~,~] = SVD_discat_0(X_n__(1+item_index_X_,:),label_Y_enum_,2,alpha_over_beta);
%%%%;

%%%%%%%%;
fname_fig = sprintf('%s/%s_USL01',dir_jpg,str_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
figure(1);clf;figmed; colormap('lines');
markersize_use = 24;
subplot(1,2,1); scatter(USL_n_(:,1),USL_n_(:,2),markersize_use,label_X_enum_,'filled');
xlabel('pc1');ylabel('pc2'); axis equal;
title('ground-truth');
subplot(1,2,2); scatter(USL_n_(:,1),USL_n_(:,2),markersize_use,label_Y_enum_,'filled');
xlabel('pc1');ylabel('pc2'); axis equal;
title('branch-split');
disp(sprintf(' %% writing %s',fname_fig));
sgtitle(str_sgtitle,'Interpreter','none');
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
fname_fig = sprintf('%s/%s_USL0',dir_jpg,str_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(1);clf;set(gcf,'Position',1+[0,0,512,1024]);
%%%%;
%ulim_ = [ min(USL_n_(:,1)) , max(USL_n_(:,1)) ];
ulim_ = [ prctile(USL_n_(:,1),1) , prctile(USL_n_(:,1),99) ];
ulim_ = mean(ulim_) + 1.25*0.5*diff(ulim_)*[-1,+1];
n_h = 64;
h_bin_ = linspace(min(ulim_),max(ulim_),n_h);
n_label_X = max(label_X_enum_);
n_label_Y = max(label_Y_enum_);
%%%%;
h_nx__ = zeros(n_h,n_label_X);
for nlabel_X=0:n_label_X-1;
tmp_index_ = efind(label_X_enum_==1+nlabel_X);
h_nx__(:,1+nlabel_X) = hist(USL_n_(1+tmp_index_,1),h_bin_);
end;%for nlabel_X=0:n_label_X-1;
h_nx__ = h_nx__ / n_item_X;
h_x_max_ = max(h_nx__,[],1);
h_x_sum_ = sum(h_nx__,1);
[~,tmp_ij] = sort(h_x_sum_,'ascend');
h_nx__ = h_nx__(:,tmp_ij);
%%%%;
h_ny__ = zeros(n_h,n_label_Y);
for nlabel_Y=0:n_label_Y-1;
tmp_index_ = efind(label_Y_enum_==1+nlabel_Y);
h_ny__(:,1+nlabel_Y) = hist(USL_n_(1+tmp_index_,1),h_bin_);
end;%for nlabel_Y=0:n_label_Y-1;
h_ny__ = h_ny__ / n_item_X;
h_y_max_ = max(h_ny__,[],1);
h_y_sum_ = sum(h_ny__,1);
[~,tmp_ij] = sort(h_y_sum_,'ascend');
h_ny__ = h_ny__(:,tmp_ij);
%%%%;
tmp_Y1_ = USL_n_(1+efind(label_Y_enum_==1),1);
tmp_Y2_ = USL_n_(1+efind(label_Y_enum_==2),1);
T = linear_discriminator_tpfp_0( tmp_Y1_ , tmp_Y2_ );
dv = max(h_nx__,[],'all')/4;
%%%%;
subplot(1,1,1);
hold on;
tmp_v = 0;
%%%%;
for nlabel_X=0:n_label_X-1;
%stairs(h_bin_,tmp_v + h_nx__(:,1+nlabel_X),'LineWidth',2,'Color','k');
p = stairfall_1(linspace(min(ulim_),max(ulim_),1+n_h),tmp_v + h_nx__(:,1+nlabel_X),tmp_v + [0,max(h_nx__(:,1+nlabel_X))],0.85*[1,1,1]); set(p,'EdgeColor','k');
tmp_v = tmp_v + dv + max(h_nx__(:,1+nlabel_X));
end;%for nlabel_X=0:n_label_X-1;
%%%%;
for nlabel_Y=0:n_label_Y-1;
if (nlabel_Y==0); tmp_str = 'r'; tmp_c_ = [1,0,0]; else tmp_str = 'b'; tmp_c_ = [0,0,1]; end;
%stairs(h_bin_,tmp_v + h_ny__(:,1+nlabel_Y),'LineWidth',2,'Color',tmp_str);
p = stairfall_1(linspace(min(ulim_),max(ulim_),1+n_h),tmp_v + h_ny__(:,1+nlabel_Y),tmp_v + [0,max(h_ny__(:,1+nlabel_Y))],tmp_c_); set(p,'EdgeColor','k');
tmp_v = tmp_v + dv + max(h_ny__(:,1+nlabel_Y));
end;%for nlabel_Y=0:n_label_Y-1;
%%%%;
plot([T,T],[0,tmp_v],'k-','LineWidth',2);
hold off;
xlim(ulim_);
xlabel('pc'); ylabel('histogram');
title('branch-split');
ylim([0,tmp_v]);
set(gca,'YTick',[],'YTickLabel',[]);
set(gca,'XTick',[],'XTickLabel',[]);
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
sgtitle(str_sgtitle,'Interpreter','none');
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%;
str_infix_A = sprintf('%s_A',str_infix);
str_infix_B = sprintf('%s_B',str_infix);
%%%%;
label_tree_SVD_discat_0( ...
 1+depth ...
,depth_upb ...
,label_tree.childA ...
,alpha_over_beta ...
,X_n__ ...
,label_X_ ...
,dir_jpg ...
,str_infix_A ...
,str_sgtitle ...
,flag_replot ...
);
%%%%;
label_tree_SVD_discat_0( ...
 1+depth ...
,depth_upb ...
,label_tree.childB ...
,alpha_over_beta ...
,X_n__ ...
,label_X_ ...
,dir_jpg ...
,str_infix_B ...
,str_sgtitle ...
,flag_replot ...
);
%%%%;
end;%if ( ~isempty(label_tree.childA) & ~isempty(label_tree.childB) );
end;%if (n_item_X>1);
end;%if (depth<=depth_upb);
