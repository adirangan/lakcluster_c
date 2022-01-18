function [un_,vn_,niteration] = rrr_0(An_,Yn_,Bn_);
% solves for beta_ such that: ;
% An_*beta_ = Yn_ ;
% where beta_ = un_*transpose(vn_). ;
% if not provided, Bn_ = pinv(An_). ;
%%%%%%%%;
if nargin<2;
rng(1);
disp(sprintf(' %% First trying general An_.'));
disp(sprintf(' %% Result should be svd of Yn_.'));
An_ = randn(12,18); Yn_ = randn(12,4);
[un_,vn_,niteration] = rrr_0(An_,Yn_);
disp(sprintf(' %% ni %d ',niteration));
[tmp_U_,tmp_S_,tmp_V_] = svds(Yn_,4);
disp(sprintf(' %% norm(Yn_ - tmp_U_*tmp_S_*transpose(tmp_V_)) = %0.16f',norm(Yn_ - tmp_U_*tmp_S_*transpose(tmp_V_),'fro')));
disp(sprintf(' %% norm(An*un_ - tmp_U_*tmp_S_) = %0.16f',norm(An_*un_ - tmp_U_(:,1)*tmp_S_(1,1),'fro')));
disp(sprintf(' %% norm(vn_ - tmp_V_) = %0.16f',norm(vn_ - tmp_V_(:,1),'fro')));
disp(sprintf(' %% Now trying singular An_'));
disp(sprintf(' %% Result should not be svd of Yn_.'));
An_ = randn(12,18); 
[tmp_U_,tmp_S_,tmp_V_] = svds(An_,6);
An_ = tmp_U_*tmp_S_*transpose(tmp_V_); %<-- now An_ is rank 6 ;
Yn_ = randn(12,4);
[un_,vn_,niteration] = rrr_0(An_,Yn_);
disp(sprintf(' %% ni %d ',niteration));
[tmp_U_,tmp_S_,tmp_V_] = svds(Yn_,4);
disp(sprintf(' %% norm(Yn_ - tmp_U_*tmp_S_*transpose(tmp_V_)) = %0.16f',norm(Yn_ - tmp_U_*tmp_S_*transpose(tmp_V_),'fro')));
disp(sprintf(' %% norm(An*un_ - tmp_U_*tmp_S_) = %0.16f',norm(An_*un_ - tmp_U_(:,1)*tmp_S_(1,1),'fro')));
disp(sprintf(' %% norm(vn_ - tmp_V_) = %0.16f',norm(vn_ - tmp_V_(:,1),'fro')));
disp('returning');return;
end;%if nargin<2;
%%%%%%%%;
verbose=1;
assert(size(An_,1)==size(Yn_,1));
n_iteration = 64;
[n_M,n_A] = size(An_);
[n_M,n_Y] = size(Yn_);
if (nargin<3); Bn_ = pinv(An_); end;
un_ = ones(n_A,1);
vn_ = ones(n_Y,1);
for niteration=1:n_iteration;
vn_old_ = vn_;
un_old_ = un_;
tmp_Anun_ = An_*un_;
tmp_utAtAnun = transpose(tmp_Anun_)*tmp_Anun_;
tmp_utAtAnun = max(1e-12,tmp_utAtAnun);
tmp_utAtYn_ = transpose(transpose(tmp_Anun_)*Yn_);
vn_ = tmp_utAtYn_/tmp_utAtAnun;
vt_ = transpose(vn_);
vv = vt_*vn_;
vn_ = vn_/sqrt(vv);
vt_ = transpose(vn_);
%[un_,nb] = conjgrad_0(An_,Yn_*vn_,un_);
un_ = Bn_*(Yn_*vn_);
tmp_norm_Yn = norm(Yn_ - An_*un_*vt_,'fro');
tmp_norm_un = norm(un_ - un_old_,'fro');
tmp_norm_vn = norm(vn_ - vn_old_,'fro');
if (verbose); disp(sprintf(' %% ni %d/%d: tmp_norm_Yn = %0.16f tmp_norm_un = %0.16f tmp_norm_vn = %0.16f',niteration,n_iteration,tmp_norm_Yn,tmp_norm_un,tmp_norm_vn)); end;
if (tmp_norm_Yn<1e-12 | (tmp_norm_un<1e-12 & tmp_norm_vn<1e-12) ); break; end;
end;%for niteration=1:n_iteration;