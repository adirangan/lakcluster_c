function ...
[ ...
 parameter ...
,a_est_ ...
,A_est__ ...
,tmp_X_dx__ ...
] = ...
affine_point_match_0( ...
 parameter ...
,X_dx__ ...
,Y_dy__ ...
,a_est_ ...
,A_est__ ...
);

str_thisfunction = 'affine_point_match_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
rng(0);
nf=0;
n_d = 2;
n_x = 128; n_x = n_x + mod(n_x,2);
n_y = n_x;
a_tru_ = randn(n_d,1);
A_tru__ = eye(n_d,n_d) + 0.1*randn(n_d,n_d);
X_dx__ = [randn(n_d,n_x/2),rand(n_d,n_x/2)+[1.5;2.0]];
Y_dy__ = a_tru_ + A_tru__*X_dx__;
parameter = struct('type','parameter');
parameter.k_use = 32;
parameter.flag_disp = 1;
lambda = 0.15;
a_est_ = lambda*a_tru_ + (1.0-lambda)*randn(n_d,1);
A_est__ = lambda*A_tru__ + (1.0-lambda)*randn(n_d,n_d);
[parameter,a_est_,A_est__,tmp_X_dx__] = ...
affine_point_match_0(parameter,X_dx__,Y_dy__,a_est_,A_est__);
disp('returning'); return;
end;%if nargin<1;

na=0;
if nargin<1+na; parameter=[]; end; na=na+1;
if nargin<1+na; X_dx__=[]; end; na=na+1;
if nargin<1+na; Y_dy__=[]; end; na=na+1;
if nargin<1+na; a_est_=[]; end; na=na+1;
if nargin<1+na; A_est__=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
flag_disp = parameter.flag_disp;
if ~isfield(parameter,'k_use'); parameter.k_use = 1; end;
k_use = parameter.k_use;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_d = size(X_dx__,1); assert(size(Y_dy__,1)==n_d);
n_x = size(X_dx__,2); n_y = size(Y_dy__,2);
if isempty(a_est_); a_est_ = zeros(n_d,1); end;
if isempty(A_est__); A_est__ = eye(n_d,n_d); end;
tmp_opt = optimset('Display','off','MaxIter',1024);
aA_est__ = [a_est_(:);A_est__(:)];
for tmp_k_use=k_use:-1:1;
if (flag_verbose>1); disp(sprintf(' %% k_use %d/%d',tmp_k_use,k_use)); end;
parameter.k_use = tmp_k_use;
tmp_f = @(aA__) affine_point_error_0(parameter,X_dx__,Y_dy__,aA__(1:n_d),aA__(n_d+(1:n_d^2)));
aA_est__ = fminsearch(tmp_f,aA_est__,tmp_opt);
a_est_ = reshape(aA_est__(1:n_d),[n_d,1]);
A_est__ = reshape(aA_est__(n_d+(1:n_d.^2)),[n_d,n_d]);
tmp_X_dx__ = a_est_ + A_est__*X_dx__;
%%%%%%%%;
if flag_disp;
figure(1);clf;figmed; 
markersize_med = 8; markersize_sml = 4;
subplot(1,2,1);
hold on; 
plot(Y_dy__(1,:),Y_dy__(2,:),'ko','MarkerSize',markersize_med,'MarkerFaceColor','g');
plot(X_dx__(1,:),X_dx__(2,:),'k^','MarkerSize',markersize_sml,'MarkerFaceColor','r');
hold off;
title('original');
%%%%;
subplot(1,2,2);
hold on; 
plot(Y_dy__(1,:),Y_dy__(2,:),'ko','MarkerSize',markersize_med,'MarkerFaceColor','g');
plot(tmp_X_dx__(1,:),tmp_X_dx__(2,:),'k^','MarkerSize',markersize_sml,'MarkerFaceColor','r');
hold off;
title('transformed');
sgtitle(sprintf('k_use %d/%d',tmp_k_use,k_use),'Interpreter','none');
drawnow();
end;%if flag_disp;
%%%%%%%%;
end;%for tmp_k_use=k_use:-1:1;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

function le = affine_point_error_0(parameter,X_dx__,Y_dy__,a_est_,A_est__);
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
tmp_X_dx__ = a_est_ + A_est__*X_dx__;
[ij_x_from_y_yk__,d_x_from_y_yk__] = knnsearch(transpose(tmp_X_dx__),transpose(Y_dy__),'K',k_use);
[ij_y_from_x_xk__,d_y_from_x_xk__] = knnsearch(transpose(Y_dy__),transpose(tmp_X_dx__),'K',k_use);
le = mean(d_x_from_y_yk__.^2,'all') + mean(d_y_from_x_xk__.^2,'all');



