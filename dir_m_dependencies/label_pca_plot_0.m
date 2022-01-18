function [AUl_n_,AUL_n_] = label_pca_plot_0(A_n_,label_);
% Applies labeled pca to A_n_ using numeric label_. ;
% Assumes A_n_ is M-by-N, with each column interpreted as a point in M-dimensional-space. ;
[M,N] = size(A_n_);
[Ul_,~,~] = svds(pca_l(A_n_,label_),3); AUl_n_ = A_n_*Ul_;
u_label_ = unique(label_);
n_label = numel(u_label_);
label__ = zeros(M,n_label);
for nlabel=0:n_label-1;
tmp_index_ = efind(label_==u_label_(1+nlabel));
label__(1+tmp_index_,1+nlabel);
end;%for nlabel=0:n_label-1;
[UL_,~,~] = svds(pca_L(A_n_,label__),3); AUL_n_ = A_n_*UL_;

if nargout<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

tmp_symbol_ = cell(n_label,1);
for nlabel=0:n_label-1;
if (mod(nlabel,5)==0); tmp_symbol = 'o'; end;
if (mod(nlabel,5)==1); tmp_symbol = '^'; end;
if (mod(nlabel,5)==2); tmp_symbol = 's'; end;
if (mod(nlabel,5)==3); tmp_symbol = 'p'; end;
if (mod(nlabel,5)==4); tmp_symbol = 'h'; end;
tmp_symbol_{1+nlabel} = tmp_symbol;
end;%for nlabel=0:n_label-1;

markersize_use = 64;
figure(1);clf; colormap(colormap_beach());
%%%%%%%%;
subplot(2,2,1);
hold on;
for nlabel=0:n_label-1;
tmp_index_ = efind(label_==u_label_(1+nlabel));
tmp_symbol = tmp_symbol_{1+nlabel};
scatter(AUl_n_(1+tmp_index_,1),AUl_n_(1+tmp_index_,2),markersize_use,label_(1+tmp_index_),tmp_symbol,'filled','MarkerEdgeColor','k');
end;%for nlabel=0:n_label-1;
hold off;
axis equal;
title('pca_l','Interpreter','none');
%%%%%%%%;
subplot(2,2,2);
hold on;
for nlabel=0:n_label-1;
tmp_index_ = efind(label_==u_label_(1+nlabel));
tmp_symbol = tmp_symbol_{1+nlabel};
scatter(AUL_n_(1+tmp_index_,1),AUL_n_(1+tmp_index_,2),markersize_use,label_(1+tmp_index_),tmp_symbol,'filled','MarkerEdgeColor','k');
end;%for nlabel=0:n_label-1;
hold off;
axis equal;
title('pca_L','Interpreter','none');
%%%%%%%%;
subplot(2,2,3);
hold on;
for nlabel=0:n_label-1;
tmp_index_ = efind(label_==u_label_(1+nlabel));
tmp_symbol = tmp_symbol_{1+nlabel};
scatter3(AUl_n_(1+tmp_index_,1),AUl_n_(1+tmp_index_,2),AUl_n_(1+tmp_index_,3),markersize_use,label_(1+tmp_index_),tmp_symbol,'filled','MarkerEdgeColor','k');
end;%for nlabel=0:n_label-1;
hold off;
axis vis3d;
title('pca_l','Interpreter','none');
%%%%%%%%;
subplot(2,2,4);
hold on;
for nlabel=0:n_label-1;
tmp_index_ = efind(label_==u_label_(1+nlabel));
tmp_symbol = tmp_symbol_{1+nlabel};
scatter3(AUL_n_(1+tmp_index_,1),AUL_n_(1+tmp_index_,2),AUL_n_(1+tmp_index_,3),markersize_use,label_(1+tmp_index_),tmp_symbol,'filled','MarkerEdgeColor','k');
end;%for nlabel=0:n_label-1;
hold off;
axis vis3d;
title('pca_L','Interpreter','none');
%%%%%%%%;
figbig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if nargout<1;

function B_ = pca_l(D_n_,k_);
D_t_ = transpose(D_n_);
[M,N] = size(D_n_); assert(numel(k_)==M); k_ = reshape(k_,[M,1]);
K_ = diag(k_);
K2_ = diag(k_.^2);
e_n_ = ones(M,1);
e_t_ = transpose(e_n_);
k1 = e_t_*k_;
k2 = e_t_.*(k_.^2);
B_ = ...
+ N * D_t_*K2_*D_n_ ...
;
return;
B_ = ...
+ N * D_t_*K2_*D_n_ ...
- 2*k1 * D_t_*K_*D_n_ ...
+ k2 * D_t_*D_n_ ...
- (D_t_*K2_*e_n_)*(e_t_*D_n_) ...
+ 2 * (D_t_*K_*e_n_)*(e_t_*K_*D_n_) ...
- (D_t_*e_n_)*(e_t_*K2_*D_n_) ...
;

function C_ = pca_L(D_n_,k__);
D_t_ = transpose(D_n_);
[M,N] = size(D_n_); assert(size(k__,1)==M);
C_ = zeros(N,N);
n_t = size(k__,2);
for nt=0:n_t-1;
C_ = C_ + pca_l(D_n_,k__(:,1+nt));
end;%for nt=0:n_t-1;


