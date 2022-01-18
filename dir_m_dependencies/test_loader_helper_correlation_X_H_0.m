function output = test_loader_helper_correlation_X_H_0(X_,H_,n_iteration);
% Assumes H_ = n_u by n_H and X_ = n_u by n_X. ;
%%%%%%%%;
if nargin<3; n_iteration = 64; end;
%%%%%%%%;
n_u = size(H_,1); n_H = size(H_,2); n_X = size(X_,2);
tmp_H_avg_ = mean(H_,1); tmp_H_std_ = max(1,std(H_,1,1));
tmp_H_ = (H_ - repmat(tmp_H_avg_,n_u,1))./repmat(tmp_H_std_,n_u,1);
tmp_X_avg_ = mean(X_,1); tmp_X_std_ = max(1,std(X_,1,1));
tmp_X_ = (X_ - repmat(tmp_X_avg_,n_u,1))./repmat(tmp_X_std_,n_u,1);
correlation_X_H__ = transpose(tmp_H_)*tmp_X_ / n_u;
[U_,S_,V_] = svds(correlation_X_H__,1);
[~,U_index_] = sort(U_,'descend'); U_index_ = U_index_ - 1;
[~,V_index_] = sort(V_,'descend'); V_index_ = V_index_ - 1;
correlation_X_H_ori__ = correlation_X_H__(1+U_index_,1+V_index_);
clear correlation_X_H__ U_ S_ V_ ;
%%%%%%%%;
correlation_X_H_avg__ = zeros(n_H,n_X);
correlation_X_H_std__ = zeros(n_H,n_X);
for niteration=1:n_iteration;
if (mod(niteration,8)==0); disp(sprintf(' %% niteration %d/%d',niteration,n_iteration)); end;
rng(1024*niteration); [tmp_Q_,~] = qr(randn(n_u)); tmp_correlation_X_H__ = transpose(tmp_H_)*tmp_Q_*tmp_X_ / n_u; clear tmp_Q_;
[tmp_U_,tmp_S_,tmp_V_] = svds(tmp_correlation_X_H__,1);
[~,tmp_U_index_] = sort(tmp_U_,'descend'); tmp_U_index_ = tmp_U_index_ - 1;
[~,tmp_V_index_] = sort(tmp_V_,'descend'); tmp_V_index_ = tmp_V_index_ - 1;
correlation_X_H_avg__ = correlation_X_H_avg__ + tmp_correlation_X_H__(1+tmp_U_index_,1+tmp_V_index_);
correlation_X_H_std__ = correlation_X_H_std__ + tmp_correlation_X_H__(1+tmp_U_index_,1+tmp_V_index_).^2;
clear tmp_correlation_X_H_ tmp_U_ tmp_S_ tmp_V_ tmp_U_index_ tmp_V_index_ ;
end;%for niteration=1:n_iteration;
correlation_X_H_avg__ = correlation_X_H_avg__/n_iteration;
correlation_X_H_std__ = sqrt(correlation_X_H_std__/n_iteration - correlation_X_H_avg__.^2);
correlation_X_H_Z__ = ( correlation_X_H_ori__ - correlation_X_H_avg__ ) ./ correlation_X_H_std__ ;
%correlation_X_H_p__ = 0.5 * erfc( correlation_X_H_Z__ / sqrt(2) ) ;
correlation_X_H_p__ = 1.0 * erfc( abs(correlation_X_H_Z__) / sqrt(2) ) ; %<-- two sided. ;
correlation_X_H_nlp__ = -erfcln( abs(correlation_X_H_Z__) / sqrt(2) ) ; %<-- two sided. ;
%%%%%%%%;
output = struct('type','correlation_X_H');
output.n_iteration = n_iteration;
output.U_index_ = U_index_;
output.V_index_ = V_index_;
output.correlation_X_H_ori__ = correlation_X_H_ori__;
output.correlation_X_H_avg__ = correlation_X_H_avg__;
output.correlation_X_H_std__ = correlation_X_H_std__;
output.correlation_X_H_Z__ = correlation_X_H_Z__;
output.correlation_X_H_p__ = correlation_X_H_p__;
output.correlation_X_H_nlp__ = correlation_X_H_nlp__;