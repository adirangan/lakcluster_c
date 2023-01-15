function ...
[ ...
 parameter ...
,a_est_ ...
,A_est__ ...
,tmp_X_dx__ ...
] = ...
affine_point_match_1( ...
 parameter ...
,X_dx__ ...
,Y_dy__ ...
,a_est_ ...
,A_est__ ...
,label_X_x_ ...
,label_Y_y_ ...
);

str_thisfunction = 'affine_point_match_1';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
rng(4);
nf=0;
n_d = 2;
n_x = 128; n_x = n_x + mod(n_x,4);
n_y = n_x;
a_tru_ = randn(n_d,1);
A_tru__ = eye(n_d,n_d) + 0.1*randn(n_d,n_d);
X_dx__ = [randn(n_d,3*n_x/4),rand(n_d,1*n_x/4)+[1.5;2.0]];
label_X_x_ = [1*ones(3*n_x/4,1);5*ones(1*n_x/4,1)];
Y_dy__ = a_tru_ + A_tru__*X_dx__;
label_Y_y_ = label_X_x_;
parameter = struct('type','parameter');
parameter.k_use = 32;
parameter.flag_disp = 1;
lambda = 0.15;
a_est_ = lambda*a_tru_ + (1.0-lambda)*randn(n_d,1);
A_est__ = lambda*A_tru__ + (1.0-lambda)*randn(n_d,n_d);
[parameter,a_est_,A_est__,tmp_X_dx__] = ...
affine_point_match_1(parameter,X_dx__,Y_dy__,a_est_,A_est__,label_X_x_,label_Y_y_);
disp('returning'); return;
end;%if nargin<1;

na=0;
if nargin<1+na; parameter=[]; end; na=na+1;
if nargin<1+na; X_dx__=[]; end; na=na+1;
if nargin<1+na; Y_dy__=[]; end; na=na+1;
if nargin<1+na; a_est_=[]; end; na=na+1;
if nargin<1+na; A_est__=[]; end; na=na+1;
if nargin<1+na; label_X_x_=[]; end; na=na+1;
if nargin<1+na; label_Y_y_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
flag_disp = parameter.flag_disp;
if ~isfield(parameter,'k_use'); parameter.k_use = 1; end;
k_use = parameter.k_use; k_use_orig = k_use;
if ~isfield(parameter,'k_gamma'); parameter.k_gamma = 0; end;
k_gamma = parameter.k_gamma;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_d = size(X_dx__,1); assert(size(Y_dy__,1)==n_d);
n_x = size(X_dx__,2); n_y = size(Y_dy__,2);
if isempty(a_est_); a_est_ = zeros(n_d,1); end;
if isempty(A_est__); A_est__ = eye(n_d,n_d); end;
if isempty(label_X_x_); label_X_x_ = zeros(n_x,1); end;
if isempty(label_Y_y_); label_Y_y_ = zeros(n_y,1); end;

tmp_opt = optimset('Display','off','MaxIter',1024);
aA_est__ = [a_est_(:);A_est__(:)];
flag_continue=1;
tmp_k_use = k_use;
while flag_continue;
if (flag_verbose>1); disp(sprintf(' %% k_use %d/%d',tmp_k_use,k_use)); end;
parameter.k_use = tmp_k_use;
tmp_f = @(aA__) affine_point_error_1(parameter,X_dx__,Y_dy__,aA__(1:n_d),aA__(n_d+(1:n_d^2)),label_X_x_,label_Y_y_);
aA_est__ = fminsearch(tmp_f,aA_est__,tmp_opt);
a_est_ = reshape(aA_est__(1:n_d),[n_d,1]);
A_est__ = reshape(aA_est__(n_d+(1:n_d.^2)),[n_d,n_d]);
tmp_X_dx__ = a_est_ + A_est__*X_dx__;
%%%%%%%%;
if flag_disp;
local_plot(parameter,X_dx__,Y_dy__,a_est_,A_est__,label_X_x_,label_Y_y_);
end;%if flag_disp;
%%%%%%%%;
if (tmp_k_use==1); flag_continue=0; end;
if (tmp_k_use> 1);
flag_continue=1;
tmp_k_use = max(1,min(tmp_k_use-1,floor(tmp_k_use*(1-k_gamma))));
end;%if (tmp_k_use> 1); 
end;%while flag_continue;

parameter.k_use = k_use_orig;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

function le = affine_point_error_1(parameter,X_dx__,Y_dy__,a_est_,A_est__,label_X_x_,label_Y_y_);
if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'k_use'); parameter.k_use = 1; end;
k_use = parameter.k_use;
n_d = size(X_dx__,1); assert(size(Y_dy__,1)==n_d);
n_x = size(X_dx__,2); n_y = size(Y_dy__,2);
a_est_ = reshape(a_est_,[n_d,1]);
A_est__ = reshape(A_est__,[n_d,n_d]);
u_label_ = union(unique(label_X_x_),unique(label_Y_y_)); n_u_label = numel(u_label_);
le = 0;
for nu_label=0:n_u_label-1;
u_label = u_label_(1+nu_label);
tmp_index_x_ = efind(label_X_x_==u_label); tmp_n_x = numel(tmp_index_x_);
tmp_index_y_ = efind(label_Y_y_==u_label); tmp_n_y = numel(tmp_index_y_);
if ( (tmp_n_x> 0) & (tmp_n_y> 0) ); %<-- we skip labels that do not match. ;
X_use_dx__ = X_dx__(:,1+tmp_index_x_);
Y_use_dy__ = Y_dy__(:,1+tmp_index_y_);
tmp_X_use_dx__ = a_est_ + A_est__*X_use_dx__;
[ij_x_from_y_yk__,d_x_from_y_yk__] = knnsearch(transpose(tmp_X_use_dx__),transpose(Y_use_dy__),'K',k_use);
[ij_y_from_x_xk__,d_y_from_x_xk__] = knnsearch(transpose(Y_use_dy__),transpose(tmp_X_use_dx__),'K',k_use);
tmp_le = sum(d_x_from_y_yk__.^2,'all')/max(1,tmp_n_y) + sum(d_y_from_x_xk__.^2,'all')/max(1,tmp_n_x);
le = le + tmp_le/max(1,k_use); %<-- note that we have normalized the contribution of each label by the number of points with that label. ;
clear X_use_dx__ Y_use_dy__ tmp_X_use_dx__ ij_x_from_y_yk__ d_x_from_y_yk__ ij_y_from_x_xk__ d_y_from_x_xk__ ;
end;%if ( (tmp_n_x> 0) & (tmp_n_y> 0) );
end;%for nu_label=0:n_u_label-1;

function local_plot(parameter,X_dx__,Y_dy__,a_est_,A_est__,label_X_x_,label_Y_y_);
if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'k_use'); parameter.k_use = 1; end;
k_use = parameter.k_use;
if ~isfield(parameter,'C_use__'); parameter.C_use__ = colormap('lines'); end;
C_use__ = parameter.C_use__;
n_d = size(X_dx__,1); assert(size(Y_dy__,1)==n_d);
n_x = size(X_dx__,2); n_y = size(Y_dy__,2);
a_est_ = reshape(a_est_,[n_d,1]);
A_est__ = reshape(A_est__,[n_d,n_d]);
n_C_use = size(C_use__,1);
figure(1);clf;figmed;
markersize_med = 8; markersize_sml = 6;
u_label_ = union(unique(label_X_x_),unique(label_Y_y_)); n_u_label = numel(u_label_);
for nu_label=0:n_u_label-1;
nc_use = max(0,min(n_C_use-1,floor(n_C_use*nu_label/max(1,n_u_label-1))));
u_label = u_label_(1+nu_label);
tmp_index_x_ = efind(label_X_x_==u_label); tmp_n_x = numel(tmp_index_x_);
tmp_index_y_ = efind(label_Y_y_==u_label); tmp_n_y = numel(tmp_index_y_);
if ( (tmp_n_x> 0) & (tmp_n_y> 0) ); %<-- we skip labels that do not match. ;
X_use_dx__ = X_dx__(:,1+tmp_index_x_);
Y_use_dy__ = Y_dy__(:,1+tmp_index_y_);
tmp_X_use_dx__ = a_est_ + A_est__*X_use_dx__;
subplot(1,2,1);
hold on; 
plot(Y_use_dy__(1,:),Y_use_dy__(2,:),'ko','MarkerSize',markersize_med,'MarkerFaceColor',C_use__(1+nc_use,:));
plot(X_use_dx__(1,:),X_use_dx__(2,:),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',C_use__(1+nc_use,:));
hold off;
subplot(1,2,2);
hold on; 
plot(Y_use_dy__(1,:),Y_use_dy__(2,:),'ko','MarkerSize',markersize_med,'MarkerFaceColor',C_use__(1+nc_use,:));
plot(tmp_X_use_dx__(1,:),tmp_X_use_dx__(2,:),'k^','MarkerSize',markersize_sml,'MarkerFaceColor',C_use__(1+nc_use,:));
hold off;
end;%if ( (tmp_n_x> 0) & (tmp_n_y> 0) );
end;%for nu_label=0:n_u_label-1;
subplot(1,2,1); grid on; title('original');
subplot(1,2,2); grid on; title('transformed');
sgtitle(sprintf('k_use %d',k_use),'Interpreter','none');
drawnow();





