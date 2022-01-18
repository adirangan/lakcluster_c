function [kld___,u_Sample_Label_] = kld_iso__0(C_,Sample_Label_);
% Returns matrix kld__ of (asymmetric) kullback-leibler divergences between pairs of clusters in C_. ;
% Clusters are defined by their Sample_Label_. ;
% kld___(nc1,nc2,j) contains the kld for clusters (nc1,nc2) calculated across columns [1:j]. ;
%%%%%%%%;
[n_Sample,n_Variable] = size(C_); d = n_Variable;
u_Sample_Label_ = unique(Sample_Label_); n_Sample_Label = length(u_Sample_Label_);
%%%%%%%%;
Sample_Label_ij_ = cell(n_Sample_Label,1);
n_Sample_Label_ = zeros(n_Sample_Label,1);
for nSample_Label=1:n_Sample_Label;
Sample_Label_ij_{nSample_Label} = find(Sample_Label_==u_Sample_Label_(nSample_Label));
n_Sample_Label_(nSample_Label) = length(Sample_Label_ij_{nSample_Label});
end;%for nSample_Label=1:n_Sample_Label;
%%%%%%%%;
mu__ = zeros(n_Sample_Label,n_Variable);
sg__ = zeros(n_Sample_Label,n_Variable);
for nSample_Label=1:n_Sample_Label;
[mu__(nSample_Label,:),sg__(nSample_Label,:)] = gaussian_iso__0(C_(Sample_Label_ij_{nSample_Label},:));
end;%for nSample_Label=1:n_Sample_Label;
%%%%%%%%;
for nSample_Label1=1:n_Sample_Label;
A_ = C_(Sample_Label_ij_{nSample_Label1},:);
n_A = n_Sample_Label_(nSample_Label1);
mu_A_ = mu__(nSample_Label1,:);
sg_A_ = sg__(nSample_Label1,:);
for nSample_Label2=nSample_Label1:n_Sample_Label;
B_ = C_(Sample_Label_ij_{nSample_Label2},:);
n_B = n_Sample_Label_(nSample_Label2);
mu_B_ = mu__(nSample_Label2,:);
sg_B_ = sg__(nSample_Label2,:);
%%%%;
A_mu_A_ = A_ - repmat(mu_A_,n_A,1);
A_mu_B_ = A_ - repmat(mu_B_,n_A,1);
B_mu_A_ = B_ - repmat(mu_A_,n_B,1);
B_mu_B_ = B_ - repmat(mu_B_,n_B,1);
kld_A_A_ = -[1:d].*log(sg_A_) - 0.5*[1:d]*log(2*pi) - mean(cumsum(abs(A_mu_A_).^2,2),1) ./(2*sg_A_.^2);
kld_A_B_ = -[1:d].*log(sg_B_) - 0.5*[1:d]*log(2*pi) - mean(cumsum(abs(A_mu_B_).^2,2),1) ./(2*sg_B_.^2);
kld_B_A_ = -[1:d].*log(sg_A_) - 0.5*[1:d]*log(2*pi) - mean(cumsum(abs(B_mu_A_).^2,2),1) ./(2*sg_A_.^2);
kld_B_B_ = -[1:d].*log(sg_B_) - 0.5*[1:d]*log(2*pi) - mean(cumsum(abs(B_mu_B_).^2,2),1) ./(2*sg_B_.^2);
kld___(nSample_Label1,nSample_Label2,:) = kld_A_B_ - kld_A_A_;
kld___(nSample_Label2,nSample_Label1,:) = kld_B_A_ - kld_B_B_;
%%%%;
end;%for nSample_Label2=1:n_Sample_Label;
end;%for nSample_Label1=1:n_Sample_Label;
%%%%%%%%;

