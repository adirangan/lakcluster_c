function ...
[ ...
 parameter ...
,f_z_ ...
,f_zp__ ...
] = ...
affine_point_match_p_0( ...
 parameter ...
,X_dx__ ...
,X_m_ ...
,Y_dy__ ...
,Y_m_ ...
);
%%%%%%%%;
% Here we assume that the labels unique(X_m_) and unique(Y_m_) are the same. ;
% We then calculate the proportion of nearest neighbors of X_dx__ (in Y_dy__) ;
% that have the same label. ;
% We also calculate the proportion of nearest neighbors of Y_dy__ (in X_dx__) ;
% that have the same label. ;
%%%%%%%%;
str_thisfunction = 'affine_point_match_p_0';

na=0;
if nargin<1+na; parameter=[]; end; na=na+1;
if nargin<1+na; X_dx__=[]; end; na=na+1;
if nargin<1+na; X_m_=[]; end; na=na+1;
if nargin<1+na; Y_dy__=[]; end; na=na+1;
if nargin<1+na; Y_m_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
flag_disp = parameter.flag_disp;
if ~isfield(parameter,'n_shuffle'); parameter.n_shuffle = 1024; end;
n_shuffle = parameter.n_shuffle;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_d = size(X_dx__,1); assert(size(Y_dy__,1)==n_d);
n_x = size(X_dx__,2); n_y = size(Y_dy__,2);
n_z = max(n_x,n_y);

if flag_disp;
figure(1);clf;figbig;
markersize_big = 12;
markersize_sml = 8;
fontsize_use = 12;
subplot(1,1,1);
hold on;
tmp_index_ = efind(Y_m_==1);
plot(Y_dy__(1+0,1+tmp_index_),Y_dy__(1+1,1+tmp_index_),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.0,1.0,1.0]);
tmp_index_ = efind(Y_m_==2);
plot(Y_dy__(1+0,1+tmp_index_),Y_dy__(1+1,1+tmp_index_),'ko','MarkerSize',markersize_big,'MarkerFaceColor',[0.5,0.0,0.5]);
tmp_index_ = efind(X_m_==1);
plot(X_dx__(1+0,1+tmp_index_),X_dx__(1+1,1+tmp_index_),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.0,0.8,0.8]);
tmp_index_ = efind(X_m_==2);
plot(X_dx__(1+0,1+tmp_index_),X_dx__(1+1,1+tmp_index_),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',[0.8,0.0,0.8]);
legend({'Y ctrl','Y case','X ctrl','X case'});
hold off;
end;%if flag_disp;

[ij_x_from_y_yk__,d_x_from_y_yk__] = knnsearch(transpose(X_dx__),transpose(Y_dy__),'K',n_x);
[ij_y_from_x_xk__,d_y_from_x_xk__] = knnsearch(transpose(Y_dy__),transpose(X_dx__),'K',n_y);
f_z_ = affine_point_match_label_0(X_m_,Y_m_,ij_x_from_y_yk__,ij_y_from_x_xk__);
f_zp__ = zeros(n_z,n_shuffle);
for nshuffle=0:n_shuffle-1;
if (flag_verbose>-1); if mod(nshuffle,64)==0; disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end; end;
rng(nshuffle);
tmp_p_x_ = randperm(n_x); tmp_p_y_ = randperm(n_y);
f_zp__(:,1+nshuffle) = affine_point_match_label_0(X_m_(tmp_p_x_),Y_m_(tmp_p_y_),ij_x_from_y_yk__,ij_y_from_x_xk__);
end;%for nshuffle=0:n_shuffle-1;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

function f_z_ = affine_point_match_label_0(X_m_,Y_m_,ij_x_from_y_yk__,ij_y_from_x_xk__);
n_x = numel(X_m_); n_y = numel(Y_m_);
X_from_Y_m_yk__ = X_m_(ij_x_from_y_yk__);
Y_from_X_m_xk__ = Y_m_(ij_y_from_x_xk__);
X_from_Y_match_yx__ = (X_from_Y_m_yk__==repmat(Y_m_,[1,n_x]));
Y_from_X_match_xy__ = (Y_from_X_m_xk__==repmat(X_m_,[1,n_y]));
f_x_ = bsxfun(@rdivide,cumsum(sum(X_from_Y_match_yx__,1)),n_y*[1:n_x]);
f_y_ = bsxfun(@rdivide,cumsum(sum(Y_from_X_match_xy__,1)),n_x*[1:n_y]);
n_z = max(n_x,n_y);
f_z_ = 0.5*interp1(1:n_x,f_x_,linspace(1,n_x,n_z)) + 0.5*interp1(1:n_y,f_y_,linspace(1,n_y,n_z));
f_z_ = reshape(f_z_,[n_z,1]);
