function [Q_n_] = dexcluster_nonbinary_AAAA_scorebox_ver0(A_n_);
% Simple differentially-expressed biclustering algorithm. ;
% Does not binarize. ;

%%%%%%%%;
if nargin<1;
M = 128; N = 256;
A_n_ = randn(M,N);
Q_n_ = dexcluster_nonbinary_AAAA_scorebox_ver0(A_n_);
m = 121; n = 237;
B_n_ = A_n_(1:m,1:n); e_t_ = ones(1,m);
QR_ = transpose((e_t_*B_n_)*transpose(B_n_)) - sum(B_n_.^2,2) ;
QC_ = (e_t_*B_n_).^2 - e_t_*(B_n_.^2) ;
QR = sum(QR_); QC = sum(QC_);
disp(sprintf(' %% Q_n_(m,n) %0.2f QR %0.2f QC %0.2f',Q_n_(m,n),QR,QC));
disp('returning');return;
end;%if nargin<1;
%%%%%%%%;

Q_n_ = cumsum(2*(cumsum(cumsum(A_n_,1).*A_n_,1)-cumsum(A_n_.^2,1)),2);
