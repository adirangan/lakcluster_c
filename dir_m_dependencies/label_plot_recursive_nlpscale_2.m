function label_plot_recursive_nlpscale_2(output_label_,nlpbra_label_,nlpnex_label_,item_label_,nlplim_,nlp_set);
% use on labels generated by test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_3. ;
na=0;
if nargin<1+na; output_label_ = []; end; na=na+1;
if nargin<1+na; nlpbra_label_ = []; end; na=na+1;
if nargin<1+na; nlpnex_label_ = []; end; na=na+1;
if nargin<1+na; item_label_ = []; end; na=na+1;
if nargin<1+na; nlplim_ = []; end; na=na+1;
if nargin<1+na; nlp_set = []; end; na=na+1;

if (isempty(output_label_)); output_label_ = []; end;
if (isempty(nlpbra_label_)); nlpbra_label_ = []; end;
if (isempty(nlpnex_label_)); nlpnex_label_ = nlpbra_label_; end;
if (isempty(item_label_)); item_label_ = []; end;
if (isempty(nlplim_)); nlplim_ = [0,27]; end;
if (isempty(nlp_set)); nlp_set = 3; end;

verbose=1;
n_item = numel(output_label_);
assert(numel(nlpbra_label_)==n_item);
assert(numel(nlpnex_label_)==n_item);
if (~isempty(item_label_)); assert(numel(item_label_)==n_item); end;
[u_output_label_,index_u_from_o_] = unique(output_label_); index_u_from_o_ = index_u_from_o_-1;
u_nlpbra_label_ = nlpbra_label_(1+index_u_from_o_);
u_nlpnex_label_ = nlpnex_label_(1+index_u_from_o_);
n_cluster = numel(u_output_label_);
if (verbose); disp(sprintf(' %% n_item %d --> n_cluster %d',n_item,n_cluster)); end;
cluster_index__ = cell(n_cluster,1);
for ncluster=0:n_cluster-1;
u_output_label = u_output_label_{1+ncluster};
cluster_index__{1+ncluster} = efind(strcmp(output_label_,u_output_label));
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
if (verbose);
for ncluster=0:n_cluster-1;
u_output_label = u_output_label_{1+ncluster};
u_nlpbra_label = u_nlpbra_label_{1+ncluster};
u_nlpnex_label = u_nlpnex_label_{1+ncluster};
cluster_index_ = cluster_index__{1+ncluster};
n_item_sub = numel(cluster_index_);
disp(sprintf(' %% ncluster %.2d/%.2d : length %.3d --> %s %s %s',ncluster,n_cluster,n_item_sub,u_output_label,u_nlpbra_label,u_nlpnex_label));
for nitem_sub=0:n_item_sub-1;
nitem = cluster_index_(1+nitem_sub);
if (~isempty(item_label_)); disp(sprintf(' %% %% %s',item_label_{1+nitem})); end;
end;%for nitem_sub=0:n_item_sub-1;
disp(sprintf(' %% '));
end;%for ncluster=0:n_cluster-1;
end;%if (verbose);
%%%%%%%%;
xgap = max(3,ceil(0.01*n_item)); xw = ceil(xgap/3);
n_x = n_item + (n_cluster-1)*xgap;
xpos_from_item_ = zeros(n_item,1);
ypos_top_from_item_ = zeros(n_item,1);
ypos_bot_from_item_ = zeros(n_item,1);
item_from_xpos_ = zeros(n_x,1);
xtick_from_xpos_ = zeros(n_x,1);
xpos=1; %<-- one based so that we can efind later. ;
for ncluster=0:n_cluster-1;
u_output_label = u_output_label_{1+ncluster};
ypos_top = length(u_output_label); %<-- vertical placement by level. ;
ypos_bot = length(u_output_label)+1; %<-- vertical placement by level. ;
u_nlpnex_label = u_nlpnex_label_{1+ncluster};
u_nlpnex_val_ = str2num(u_nlpnex_label);
ypos_top = u_nlpnex_val_(end); %<-- vertical placement by nlp. ;
ypos_bot = 0; %<-- vertical placement by nlp. ;
cluster_index_ = cluster_index__{1+ncluster};
n_item_sub = numel(cluster_index_);
for nitem_sub=0:n_item_sub-1;
nitem = cluster_index_(1+nitem_sub);
xpos_from_item_(1+nitem) = xpos;
ypos_top_from_item_(1+nitem) = ypos_top;
ypos_bot_from_item_(1+nitem) = ypos_bot;
item_from_xpos_(1+xpos) = nitem;
xtick_from_xpos_(1+xpos) = xpos;
xpos=xpos+1;
end;%for nitem_sub=0:n_item_sub-1;
xpos=xpos+xgap;
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
node_label_ = cell(0,0);
node_isroot_ = zeros(0,0);
node_ischildA_ = zeros(0,0);
node_ischildB_ = zeros(0,0);
node_parent_label_ = cell(0,0);
node_parent_nlpbra_ = zeros(0,0);
node_parent_nlpnex_ = zeros(0,0);
node_isleaf_ = zeros(0,0);
node_childA_nlpbra_ = zeros(0,0);
node_childA_nlpnex_ = zeros(0,0);
node_childB_nlpbra_ = zeros(0,0);
node_childB_nlpnex_ = zeros(0,0);
node_childA_label_ = cell(0,0);
node_childB_label_ = cell(0,0);
node_nlpbra_ = zeros(0,0);
node_nlpnex_ = zeros(0,0);
node_ij_ = zeros(0,0);
node_leaf_distance_ = zeros(0,0);
node_root_distance_ = zeros(0,0);
n_node=0;
for ncluster=0:n_cluster-1;
u_output_label = u_output_label_{1+ncluster};
u_nlpbra_label = u_nlpbra_label_{1+ncluster};
u_nlpnex_label = u_nlpnex_label_{1+ncluster};
u_nlpbra_ = str2num(u_nlpbra_label);
u_nlpnex_ = str2num(u_nlpnex_label);
n_y = length(u_output_label)-1;
n_0 = strfind(u_output_label,'0');
assert(n_0==n_y+1); %<-- 0 termination. ;
assert(numel(u_nlpbra_)==n_y);
assert(numel(u_nlpnex_)==n_y);
for ny=-1:n_y-1;
node_root_distance = 1+ny;
node_leaf_distance = n_y-1-ny;
node_label = sprintf('%s0',u_output_label(1+(0:ny)));
node_isroot = 0;
node_ischildA = 0;
node_ischildB = 0;
node_parent_nlpbra = 0;
node_parent_nlpnex = 0;
if (ny==-1);
node_isroot = 1;
node_ischildA = 0;
node_ischildB = 0;
node_parent_label = '';
node_parent_nlpbra = max(min(nlplim_),min(max(nlplim_),+Inf));
node_parent_nlpnex = max(min(nlplim_),min(max(nlplim_),+Inf));
end;%if (ny==-1);
if (ny>=0);
node_isroot = 0;
node_parent_label = sprintf('%s0',u_output_label(1+(0:ny-1)));
node_ischildA = (u_output_label(1+ny)=='A');
node_ischildB = (u_output_label(1+ny)=='B');
node_parent_nlpbra = max(min(nlplim_),min(max(nlplim_),u_nlpbra_(1+0+ny)));
node_parent_nlpnex = max(min(nlplim_),min(max(nlplim_),u_nlpnex_(1+0+ny)));
end;%if (ny>=0);
node_isleaf = 1;
node_childA_label = '';
node_childB_label = '';
node_childA_nlpbra = max(min(nlplim_),min(max(nlplim_),nlp_set));
node_childA_nlpnex = max(min(nlplim_),min(max(nlplim_),nlp_set));
node_childB_nlpbra = max(min(nlplim_),min(max(nlplim_),nlp_set));
node_childB_nlpnex = max(min(nlplim_),min(max(nlplim_),nlp_set));
node_nlpbra = max(min(nlplim_),min(max(nlplim_),nlp_set));
node_nlpnex = max(min(nlplim_),min(max(nlplim_),nlp_set));
if (ny<n_y-1);
node_isleaf = 0;
node_childA_label = sprintf('%sA0',u_output_label(1+(0:ny)));
node_childB_label = sprintf('%sB0',u_output_label(1+(0:ny)));
node_childA_nlpbra = 0; %<-- fill in later. ;
node_childA_nlpnex = 0; %<-- fill in later. ;
node_childB_nlpbra = 0; %<-- fill in later. ;
node_childB_nlpnex = 0; %<-- fill in later. ;
node_nlpbra = max(min(nlplim_),min(max(nlplim_),u_nlpbra_(1+1+ny)));
node_nlpnex = max(min(nlplim_),min(max(nlplim_),u_nlpnex_(1+1+ny)));
end;%if (ny<n_y-1);
node_ij = 0; if (ny==n_y-1); node_ij = 1+ncluster; end;
if (isempty(intersect(node_label_,{node_label})));
node_label_{1+n_node} = node_label;
node_isroot_(1+n_node) = node_isroot;
node_ischildA_(1+n_node) = node_ischildA;
node_ischildB_(1+n_node) = node_ischildB;
node_parent_label_{1+n_node} = node_parent_label;
node_parent_nlpbra_(1+n_node) = node_parent_nlpbra;
node_parent_nlpnex_(1+n_node) = node_parent_nlpnex;
node_isleaf_(1+n_node) = node_isleaf;
node_childA_label_{1+n_node} = node_childA_label;
node_childB_label_{1+n_node} = node_childB_label;
node_childA_nlpbra_(1+n_node) = node_childA_nlpbra;
node_childA_nlpnex_(1+n_node) = node_childA_nlpnex;
node_childB_nlpbra_(1+n_node) = node_childB_nlpbra;
node_childB_nlpnex_(1+n_node) = node_childB_nlpnex;
node_nlpbra_(1+n_node) = node_nlpbra;
node_nlpnex_(1+n_node) = node_nlpnex;
node_ij_(1+n_node) = node_ij;
node_root_distance_(1+n_node) = node_root_distance;
node_leaf_distance_(1+n_node) = node_leaf_distance;
n_node=n_node+1;
end;%if (isempty(intersect(node_label_,{node_label})));
end;%for ny=-1:n_y-1;
end;%for ncluster=0:n_cluster-1;
assert(n_node==numel(node_label_));
assert(n_node==numel(node_parent_nlpbra_));
assert(n_node==numel(node_parent_nlpnex_));
assert(n_node==numel(node_nlpbra_));
assert(n_node==numel(node_childA_nlpnex_));
assert(n_node==numel(node_childA_nlpbra_));
assert(n_node==numel(node_childB_nlpnex_));
assert(n_node==numel(node_childB_nlpbra_));
assert(n_node==numel(node_nlpnex_));
assert(n_node==numel(node_ij_));
%%%%%%%%;
node_n_a_ = zeros(n_node,1);
node_parent_index_ = zeros(n_node,1);
node_childA_index_ = zeros(n_node,1);
node_childB_index_ = zeros(n_node,1);
node_xpos_avg_ = zeros(n_node,1);
node_xpos_max_ = zeros(n_node,1);
node_xpos_min_ = zeros(n_node,1);
node_ypos_top_ = zeros(n_node,1);
node_ypos_bot_ = zeros(n_node,1);
for nnode=0:n_node-1;
node_label = node_label_{1+nnode};
node_isroot = node_isroot_(1+nnode);
node_ischildA = node_ischildA_(1+nnode);
node_ischildB = node_ischildB_(1+nnode);
if (~node_isroot);
node_parent_index_(1+nnode) = efind(strcmp(node_label_,node_parent_label_{1+nnode}));
if (node_ischildA); %<-- childA. ;
node_childA_nlpbra_(1+node_parent_index_(1+nnode)) = max(node_childA_nlpbra_(1+node_parent_index_(1+nnode)),node_nlpbra_(1+nnode));
node_childA_nlpnex_(1+node_parent_index_(1+nnode)) = max(node_childA_nlpnex_(1+node_parent_index_(1+nnode)),node_nlpnex_(1+nnode));
end;%if (node_ischildA); %<-- childA. ;
if (node_ischildB); %<-- childB. ;
node_childB_nlpbra_(1+node_parent_index_(1+nnode)) = max(node_childB_nlpbra_(1+node_parent_index_(1+nnode)),node_nlpbra_(1+nnode));
node_childB_nlpnex_(1+node_parent_index_(1+nnode)) = max(node_childB_nlpnex_(1+node_parent_index_(1+nnode)),node_nlpnex_(1+nnode));
end;%if (node_ischildB); %<-- childB. ;
end;%if (~node_isroot); %<-- inefficient. ;
node_isleaf = node_isleaf_(1+nnode);
if (~node_isleaf); node_childA_index_(1+nnode) = efind(strcmp(node_label_,node_childA_label_{1+nnode})); end; %<-- inefficient. ;
if (~node_isleaf); node_childB_index_(1+nnode) = efind(strcmp(node_label_,node_childB_label_{1+nnode})); end; %<-- inefficient. ;
n_y = numel(node_label);
if (n_y==1);
n_a = n_item;
tmp_index_ = 0:n_item-1;
else;%if (n_y==1);
n_a=0;
for nitem=0:n_item-1; 
output_label = output_label_{1+nitem};
if (min(strfind(output_label,node_label(1:end-1)))==1);
n_a = n_a+1;
end;%if (min(strfind(output_label,node_label(1:end-1)))==1);
end;%for nitem=0:n_item-1; 
tmp_index_ = zeros(n_a,1);
na=0;
for nitem=0:n_item-1; 
output_label = output_label_{1+nitem};
if (min(strfind(output_label,node_label(1:end-1)))==1);
tmp_index_(1+na) = nitem;
na=na+1;
end;%if (min(strfind(output_label,node_label(1:end-1)))==1);
end;%for nitem=0:n_item-1; 
end;%if (n_y==1);
node_n_a_(1+nnode) = n_a;
node_xpos_avg_(1+nnode) = mean(xpos_from_item_(1+tmp_index_));
node_xpos_max_(1+nnode) = max(xpos_from_item_(1+tmp_index_));
node_xpos_min_(1+nnode) = min(xpos_from_item_(1+tmp_index_));
%node_ypos_(1+nnode) = n_y; %<-- vertical placement by tree level. ;
node_ypos_top_(1+nnode) = node_parent_nlpbra_(1+nnode); %<-- vertical placement by p-value. ;
node_ypos_bot_(1+nnode) = node_nlpbra_(1+nnode); %<-- vertical placement by p-value. ;
end;%for nnode=0:n_node-1;
%%%%%%%%;
node_xpos_parent_ = zeros(n_node,1);
node_xpos_childA_ = zeros(n_node,1);
node_xpos_childB_ = zeros(n_node,1);
[~,nnode_root_last_] = sort(node_root_distance_,'descend'); nnode_root_last_ = nnode_root_last_ - 1;
for nl=0:n_node-1;
nnode = nnode_root_last_(1+nl);
node_label = node_label_{1+nnode};
node_isroot = node_isroot_(1+nnode);
node_parent_index = node_parent_index_(1+nnode);
node_isleaf = node_isleaf_(1+nnode);
node_childA_index = node_childA_index_(1+nnode);
node_childB_index = node_childB_index_(1+nnode);
if ( node_isleaf);
node_xpos_childA_(1+nnode) = node_xpos_min_(1+nnode);
node_xpos_childB_(1+nnode) = node_xpos_max_(1+nnode);
end;%if ( node_isleaf);
if (~node_isleaf);
node_xpos_childA_(1+nnode) = node_xpos_parent_(1+node_childA_index);
node_xpos_childB_(1+nnode) = node_xpos_parent_(1+node_childB_index);
end;%if (~node_isleaf);
node_xpos_parent_(1+nnode) = 0.5*(node_xpos_childA_(1+nnode) + node_xpos_childB_(1+nnode));
end;%for nl=0:n_node-1;
%%%%%%%%;
if (verbose);
for nnode=0:n_node-1;
node_label = node_label_{1+nnode};
node_isroot = node_isroot_(1+nnode);
node_parent_label = node_parent_label_{1+nnode};
node_parent_index = node_parent_index_(1+nnode);
node_isleaf = node_isleaf_(1+nnode);
node_childA_label = node_childA_label_{1+nnode};
node_childA_index = node_childA_index_(1+nnode);
node_childB_label = node_childB_label_{1+nnode};
node_childB_index = node_childB_index_(1+nnode);
node_nlpbra = node_nlpbra_(1+nnode);
node_nlpnex = node_nlpnex_(1+nnode);
ypos_top = node_ypos_top_(1+nnode);
ypos_bot = node_ypos_bot_(1+nnode);
xpos_avg = node_xpos_avg_(1+nnode);
xpos_min = node_xpos_min_(1+nnode);
xpos_max = node_xpos_max_(1+nnode);
node_xpos_childA = node_xpos_childA_(1+nnode);
node_xpos_childB = node_xpos_childB_(1+nnode);
node_xpos_parent = node_xpos_parent_(1+nnode);
n_a = node_n_a_(1+nnode);
disp(sprintf(' %% node %.3d/%.3d --> %s (r%dl%d) (parent %d %s, children %d %s %d %s) [nlpbra %+0.2f nlpnex %+0.2f] a %.3d --> [x,y] = [(%+0.2f,%+0.2f,%+0.2f),%+0.2f,%+0.2f]',nnode,n_node,node_label,node_isroot,node_isleaf,node_parent_index,node_parent_label,node_childA_index,node_childA_label,node_childB_index,node_childB_label,node_nlpbra,node_nlpnex,n_a,xpos_min,xpos_avg,xpos_max,ypos_top,ypos_bot));
end;%for nnode=0:n_node-1;
end;%if (verbose);
%%%%%%%%;

n_item_threshold = 1024;
n_node_threshold = 128;

cla;
hold on;
%ypos_upb = max(ypos_from_item_); %<-- vertical position from tree level. ;
ypos_upb = max(nlplim_); %<-- vertical position from nlp. ;
xpos_upb = max(xpos_from_item_);
xpos_lob = min(xpos_from_item_);
c_ = colormap_nlpvt(64,(3-min(nlplim_))/diff(nlplim_)); n_c = size(c_,1);
%%%%%%%%;
for nnode=0:n_node-1;
node_isleaf = node_isleaf_(1+nnode);
node_nlpbra = node_nlpbra_(1+nnode);
node_nlpnex = node_nlpnex_(1+nnode); if (node_isleaf); node_nlpnex = node_nlpnex_(1+node_parent_index_(1+nnode)); end;
node_xpos_childA = node_xpos_childA_(1+nnode);
node_xpos_childB = node_xpos_childB_(1+nnode);
node_xpos_parent = node_xpos_parent_(1+nnode);
node_ij = node_ij_(1+nnode);
if ((node_nlpnex>=10) & (node_isleaf==0)); str_tmp = sprintf('%.2d',min(99,floor(node_nlpnex))); end;
if ((node_nlpnex< 10) & (node_isleaf==0)); str_tmp = sprintf('%.1f',floor(10*node_nlpnex)/10); end;
if (node_isleaf==1); str_tmp = sprintf('[%d]',node_ij); end;
xpos_avg = node_xpos_avg_(1+nnode);
xpos_min = node_xpos_min_(1+nnode);
xpos_max = node_xpos_max_(1+nnode);
ypos_top = node_ypos_top_(1+nnode);
ypos_bot = node_ypos_bot_(1+nnode);
ypos_parent = node_parent_nlpnex_(1+nnode);
ypos_00self = node_nlpnex_(1+nnode);
ypos_childA = node_childA_nlpnex_(1+nnode);
ypos_childB = node_childB_nlpnex_(1+nnode);
nc = max(0,min(n_c-1,floor(n_c*(node_nlpnex - min(nlplim_))/diff(nlplim_))));
if ( node_isleaf); 
% triangle: ;
tmp_x_ = [ ...
 node_xpos_childA - xw ...
;node_xpos_childB + xw ...
;node_xpos_parent + xw ...
;node_xpos_parent - xw ...
;node_xpos_childA - xw...
];
tmp_y_ = [ ...
 ypos_00self ...
;ypos_00self ...
;ypos_parent ...
;ypos_parent ...
;ypos_00self ...
];
patch(-0.5+1+tmp_x_,tmp_y_,c_(1+nc,:),'LineStyle','-','LineWidth',1,'EdgeColor',c_(1+nc,:));
% box: ;
tmp_x_ = [ ...
 node_xpos_childA - xw ...
;node_xpos_childB + xw ...
;node_xpos_childB + xw ...
;node_xpos_childA - xw...
];
tmp_y_ = [ ...
 0 ...
;0 ...
;ypos_00self ...
;ypos_00self ...
];
patch(-0.5+1+tmp_x_,tmp_y_,c_(1+nc,:),'LineStyle','-','LineWidth',1,'EdgeColor',c_(1+nc,:));
end;%if ( node_isleaf); 
if (~node_isleaf); 
node_childA_index = node_childA_index_(1+nnode);
node_childA_xpos_childA = node_xpos_childA_(1+node_childA_index);
node_childA_xpos_childB = node_xpos_childB_(1+node_childA_index);
node_childB_index = node_childB_index_(1+nnode);
node_childB_xpos_childA = node_xpos_childA_(1+node_childB_index);
node_childB_xpos_childB = node_xpos_childB_(1+node_childB_index);
% staple: ;
tmp_x_ = [ ...
 node_xpos_childA - xw ...
;node_xpos_childA + xw ...
;node_xpos_childA + xw ...
;node_xpos_childB - xw ...
;node_xpos_childB - xw ...
;node_xpos_childB + xw ...
;node_xpos_childB + xw ...
;node_xpos_childA - xw...
];
tmp_yw = min( 2*xw*ypos_upb/xpos_upb , ypos_00self - max(ypos_childA,ypos_childB) );
tmp_yw = 2*xw*ypos_upb/xpos_upb;
tmp_y_ = [ ...
 ypos_childA ...
;ypos_childA ...
;ypos_00self - tmp_yw ...
;ypos_00self - tmp_yw ...
;ypos_childB ...
;ypos_childB ...
;ypos_00self ...
;ypos_00self ...
];
%patch(-0.5+1+tmp_x_,tmp_y_,c_(1+nc,:),'LineStyle','-','LineWidth',1,'EdgeColor',c_(1+nc,:));
% sawhorse;
tmp_x_ = [ ...
 node_childA_xpos_childA - xw ...
;node_childA_xpos_childB + xw ...
;node_xpos_childA + xw ...
;node_xpos_childB - xw ...
;node_childB_xpos_childA - xw ...
;node_childB_xpos_childB + xw ...
;node_xpos_childB + xw ...
;node_xpos_childA - xw...
];
tmp_yw = min( 2*xw*ypos_upb/xpos_upb , ypos_00self - max(ypos_childA,ypos_childB) );
tmp_yw = 2*xw*ypos_upb/xpos_upb;
tmp_y_ = [ ...
 ypos_childA ...
;ypos_childA ...
;ypos_00self - tmp_yw ...
;ypos_00self - tmp_yw ...
;ypos_childB ...
;ypos_childB ...
;ypos_00self ...
;ypos_00self ...
];
patch(-0.5+1+tmp_x_,tmp_y_,c_(1+nc,:),'LineStyle','-','LineWidth',1,'EdgeColor',c_(1+nc,:));
end;%if (~node_isleaf); 
end;%for nnode=0:n_node-1;
%%%%%%%%;
% colorbar;
%%%%%%%%;
colorbar_width = (1.0/n_c) * (1+1+xpos_upb + 2*xw - xpos_lob) ;
tmp_x_ = 1+1+xpos_upb+xw + 2*colorbar_width + colorbar_width*repmat([0;1;1;0],[1,n_c]);
tmp_y_ = linspace(0,1+ypos_upb,n_c+1); tmp_y_ = [ tmp_y_(1:n_c+0) ; tmp_y_(1:n_c+0) ; tmp_y_(2:n_c+1) ; tmp_y_(2:n_c+1) ] ;
patch(tmp_x_,tmp_y_,reshape(c_,[1,n_c,3]),'LineStyle','none');
n_l = 9;
for nl=0:n_l;
text(1+1+xpos_upb+xw + 3*colorbar_width + colorbar_width*0.25,(1+ypos_upb)*nl/n_l,sprintf('%0.2f',min(nlplim_) + diff(nlplim_)*nl/n_l));
tmp_ = text(1+1+xpos_upb+xw + 1.25*colorbar_width,(1+ypos_upb)/2,'-log(p)'); set(tmp_,'Rotation',90);
end;%for nl=0:n_l;
%%%%%%%%;
set(gca,'YTick',[]);
tmp_ = text(xpos_lob-xw-2*colorbar_width,(1+ypos_upb)/2,'nlp'); set(tmp_,'Rotation',90);
for nl=0:1+ypos_upb;
tmp_ = text(xpos_lob-xw-1*colorbar_width,(1+ypos_upb)*nl/(1+ypos_upb),sprintf('%d',nl));
end;%for nl=0:1+ypos_upb;
%%%%%%%%;
hold off;
set(gca,'TickLength',[0,0]); box off;
if ( isempty(item_label_));
tmp_ = gca; tmp_.XColor = 'w'; tmp_.YColor = 'w';
set(gca,'XTick',[]); set(gca,'TickLength',[0,0]);
end;%if ( isempty(item_label_));
if (~isempty(item_label_));
%tmp_ = gca; tmp_.XColor = 'w'; tmp_.YColor = 'w';
set(gca,'XTick',1+efind(xtick_from_xpos_),'XTickLabel',item_label_(1+item_from_xpos_(1+efind(xtick_from_xpos_)))); xtickangle(90);
end;%if (~isempty(item_label_)); 
set(gca,'FontSize',6);
xlim([xpos_lob-xw - 2.5*colorbar_width,1+1+xpos_upb+xw + 3*colorbar_width]);
ylim([0,1+ypos_upb]);
set(gcf,'Color',[1,1,1]);
figbig;



