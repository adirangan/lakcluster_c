function ...
[ ...
 parameter ...
,vlXY_uxy___ ...
,vcXY_uxy___ ...
,u_L_u_ ...
] = ...
affine_point_match_vector_label_0( ...
 parameter ...
,X_dx__ ...
,Y_dy__ ...
,L_y_ ...
);
%%%%%%%%;
% Here we use use: ;
% 1. the labels (e.g., ctrls,case,bicl) stored in L_y_ ;
% 2. as well as the 'training' data-points Y_dy__ ;
% to provide approximate labels for the 'testing' data-points X_dx__. ;
%%%%%%%%;
% These approximate labels will be expressed in the form: ;
% vlXY_uxy___: a double-array of size numel(unique(L_y_))-by-size(X_dx__,2)-by-size(Y_dy__,2) ;
% where vlXY_uxy___(1+nu,1+nx,1+ny) is the strength of association of: ;
% label nu to point nx using ny nearest neighbors from Y. ;
% More specifically: ;
% Assume we define k = 1+ny ; %<-- this is the number of nearest-neighbors we use to define the vector_label. ;
% and set: vlXY_ux__ = vlXY_ux___(:,:,1+ny), ;
% and then define: u_L_u_ = unique(L_y_) ; %<-- this is a listing of the unique labels from Y. ;
% then vlXY_ux__(:,1+nx) ;
% will be the vector label (of length numel(u_L_u_)) associated with point nx defined by X_dx__(:,1+nx). ;
% The value vlXY_ux__(1+nu,1+nx) ;
% will be the value of this vector_label attributed to category u_L_u_(1+nu). ;
% Note that the putative labels of X are not involved in this process. ;
% After the vector_labels are assigned one can subselect (and renormalize) to obtain 'conditional' labels as necessary. ;
%%%%%%%%;
% In terms of actual calculation, ;
% vcXY_uxy___ holds the counts (i.e., numbers) of each label at each nx and k, ;
% while vlXY_uxy___ normalizes these counts by k. ;
%%%%%%%%;
str_thisfunction = 'affine_point_match_vector_label_0';

na=0;
if nargin<1+na; parameter=[]; end; na=na+1;
if nargin<1+na; X_dx__=[]; end; na=na+1;
if nargin<1+na; Y_dy__=[]; end; na=na+1;
if nargin<1+na; L_y_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_d = size(X_dx__,1); assert(size(Y_dy__,1)==n_d);
n_x = size(X_dx__,2); n_y = size(Y_dy__,2);
u_L_u_ = unique(L_y_);
n_u_L_u = numel(u_L_u_);
L_uy__ = sparse(n_u_L_u,n_y);
for nu_L_u=0:n_u_L_u-1;
u_L_u = u_L_u_(1+nu_L_u);
tmp_index_ = efind(L_y_==u_L_u);
L_uy__(1+nu_L_u,1+tmp_index_) = 1;
clear tmp_index_;
end;%for nu_L_u=0:n_u_L_u-1;

[ij_y_from_x_xk__,d_y_from_x_xk__] = knnsearch(transpose(Y_dy__),transpose(X_dx__),'K',n_y);
vlXY_uxy___ = zeros(n_u_L_u,n_x,n_y);
vcXY_uxy___ = zeros(n_u_L_u,n_x,n_y);
for nx=0:n_x-1;
ij_y_k_ = ij_y_from_x_xk__(1+nx,:);
tmp_L_uk__ = L_uy__(:,ij_y_k_);
tmp_L_cumsum_uk__ = cumsum(tmp_L_uk__,2);
vcXY_uk__ = tmp_L_cumsum_uk__;
vlXY_uk__ = bsxfun(@rdivide,vcXY_uk__,1:n_y);
vcXY_uxy___(:,1+nx,:) = vcXY_uk__;
vlXY_uxy___(:,1+nx,:) = vlXY_uk__;
clear ij_y_k_ tmp_L_uk__ tmp_L_cumsum_uk__ vcXY_uk__ vlXY_uk__ ;
end;%for nx=0:n_x-1;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
