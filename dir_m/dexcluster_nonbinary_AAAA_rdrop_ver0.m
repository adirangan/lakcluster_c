function [out_xdrop_,trace_] = dexcluster_nonbinary_AAAA_rdrop_ver0(A_n_,gamma);
% Simple differentially-expressed biclustering algorithm. ;
% Does not binarize. ;
% Only drops rows (not columns). ;

if nargin<2; gamma = 0.0; end;
verbose=1;
[n_r,n_c] = size(A_n_);
if min(n_r,n_c)<1024; verbose=0; end;

[rdrop_,cdrop_,lrij_,lcij_] = get_xdrop_(n_r,n_c,gamma);
cdrop_ = zeros(size(cdrop_)); %<-- we do not drop any columns. ;
r_ij_ = 1:n_r; c_ij_ = 1:n_c;
n_a = length(rdrop_);
e_t_ = ones(1,n_r);
%out_xdrop_ = zeros(n_r+n_c,2);
out_xdrop_ = zeros(n_r,2);
trace_ = zeros(n_a,6); %<-- [iteration , length(r_ij_) , length(c_ij_) , mean(QR_) , mean(QC_) , 1] ;
etAn_ = e_t_*A_n_; etAn2_ = e_t_*(A_n_.^2);
QR_pre_ = transpose((e_t_*A_n_)*transpose(A_n_)) - sum(A_n_.^2,2) ;
QC_pre_ = (e_t_*A_n_).^2 - e_t_*(A_n_.^2) ;
%%%%%%%%;
% A_n_ represented as: ;
% [ w_n_ y_n_ ] ;
% [ x_n_ z_n_ ] ;
% where x_n_ and y_n_ are removed. ;
%%%%%%%%;
nx=0;
tic;
t_s = toc;
for na=1:n_a;
if (verbose);
if (mod(na,10)==0); 
t_a = toc; t_e = t_a*n_a/na; 
disp(sprintf(' %% na %d/%d, time %0.2fs, estimated %0.2fs (%0.2fh)',na,n_a,t_a,t_e,t_e/3600)); 
end;%if (mod(na,10)==0); 
end;%if (verbose);
trace_(na,:) = [na , length(r_ij_) , length(c_ij_) , mean(QR_pre_) , mean(QC_pre_) , 1.0];
[~,tmp_r_ij_] = sort(QR_pre_(r_ij_));
r_rmv_ = r_ij_(tmp_r_ij_(1:rdrop_(na))); r_rtn_ = r_ij_(tmp_r_ij_(rdrop_(na)+1:end));
tmp_x_n_ = A_n_(r_rmv_,c_ij_); tmp_etxn_ = sum(tmp_x_n_,1); tmp_etxn2_ = sum(tmp_x_n_.^2,1);
tmp_etwn_ = etAn_(c_ij_) - tmp_etxn_; tmp_etwn2_ = etAn2_(c_ij_) - tmp_etxn2_;
etAn_(c_ij_) = etAn_(c_ij_) - tmp_etxn_; etAn2_(c_ij_) = etAn2_(c_ij_) - tmp_etxn2_;
QC_pos_ = zeros(1,n_c); QC_pos_(c_ij_) = QC_pre_(c_ij_) - (2*tmp_etwn_.*tmp_etxn_ + tmp_etxn_.^2 - tmp_etxn2_);
QC_pre_ = QC_pos_;
[~,tmp_c_ij_] = sort(QC_pre_(c_ij_));
c_rmv_ = c_ij_(tmp_c_ij_(1:cdrop_(na))); c_rtn_ = c_ij_(tmp_c_ij_(cdrop_(na)+1:end));
tmp_w_n_ = A_n_(r_rtn_,c_rtn_); tmp_y_n_ = A_n_(r_rtn_,c_rmv_);
tmp_x_n_ = A_n_(r_rmv_,c_rtn_); tmp_z_n_ = A_n_(r_rmv_,c_rmv_);
tmp_ynyten_ = tmp_y_n_*(transpose(sum(tmp_y_n_,1)));
tmp_wnxten_ = tmp_w_n_*(transpose(sum(tmp_x_n_,1)));
tmp_ynzten_ = tmp_y_n_*(transpose(sum(tmp_z_n_,1)));
tmp_yn2fn_ = sum(tmp_y_n_.^2,2);
QR_pos_ = zeros(n_r,1); QR_pos_(r_rtn_) = QR_pre_(r_rtn_) - (tmp_ynyten_ + tmp_wnxten_ + tmp_ynzten_ - tmp_yn2fn_);
QR_pre_ = QR_pos_;
etAn_(c_rmv_) = 0; etAn2_(c_rmv_) = 0;
QC_pre_(c_rmv_) = 0;
r_ij_(tmp_r_ij_(1:rdrop_(na))) = [];
c_ij_(tmp_c_ij_(1:cdrop_(na))) = [];
out_xdrop_(nx + (1:rdrop_(na)),:) = [r_rmv_(:)-1 , -ones(rdrop_(na),1)]; nx = nx+rdrop_(na);
%out_xdrop_(nx + (1:cdrop_(na)),:) = [-ones(cdrop_(na),1) , c_rmv_(:)-1]; nx = nx+cdrop_(na);
end;%for na=1:n_a;
