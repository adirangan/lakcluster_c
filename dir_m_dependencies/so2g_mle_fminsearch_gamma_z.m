function [p_threshold,gamma_z_out,param_out_,log_likelihood_ratio,label_,label_auc] = so2g_mle_fminsearch_gamma_z(d_);
% Uses so2g_ml3_fminsearch for each projection of 2-dimensional data d_. ;

verbose = 0;
if (verbose>0); disp(sprintf(' %% [entering so2g_mle_fminsearch_gamma_z]')); end;

n_gamma_z = 64; gamma_z_ = linspace(0,pi,n_gamma_z+1); gamma_z_ = gamma_z_(1:end-1);
llr_ = zeros(n_gamma_z,1);
for ngamma_z=1:n_gamma_z;
gamma_z = gamma_z_(ngamma_z);
if (verbose>0); disp(sprintf(' %% ngamma_z %d/%d',ngamma_z,n_gamma_z)); end;
llr_(ngamma_z) = -nllr(verbose,d_,gamma_z);
end;%for ngamma_z=1:n_gamma_z;
if (verbose>0); plot(gamma_z_,llr_,'k.-'); end;
[~,tmp_ij] = max(llr_);
gamma_z = gamma_z_(tmp_ij);
if (verbose<=0); gamma_z = fminsearch(@(gamma_z) nllr(verbose,d_,gamma_z),gamma_z); end;
if (verbose>=1); gamma_z = fminsearch(@(gamma_z) nllr(verbose,d_,gamma_z),gamma_z,optimset('Display','iter')); end;
gamma_z_out = gamma_z;
tmp_nu_ = [cos(gamma_z);sin(gamma_z)];
tmp_d_ = d_(:,1)*tmp_nu_(1) + d_(:,2)*tmp_nu_(2);
tmp_param_0in_(1) = prctile(tmp_d_,25);
tmp_param_0in_(2) = std(tmp_d_)/2;
tmp_param_0in_(3) = prctile(tmp_d_,75);
tmp_param_0in_(4) = std(tmp_d_)/2;
tmp_param_0in_(5) = 0.0;
[p_threshold,param_out_,log_likelihood_ratio,label_,label_auc] = so2g_mle_fminsearch(tmp_d_,tmp_param_0in_,0);

if (verbose>0); disp(sprintf(' %% [finished so2g_mle_fminsearch_gamma_z]')); end;

function output = nllr(verbose,d_,gamma_z);
if (verbose>0); disp(sprintf(' %% [entering nlrr] gamma_z %0.2f',gamma_z)); end;
tmp_nu_ = [cos(gamma_z);sin(gamma_z)];
tmp_d_ = d_(:,1)*tmp_nu_(1) + d_(:,2)*tmp_nu_(2);
tmp_param_0in_(1) = prctile(tmp_d_,25);
tmp_param_0in_(2) = std(tmp_d_)/2;
tmp_param_0in_(3) = prctile(tmp_d_,75);
tmp_param_0in_(4) = std(tmp_d_)/2;
tmp_param_0in_(5) = 0.0;
tmp_options = optimset('Display','iter');
if (verbose<=0); [~,~,output] = so2g_mle_fminsearch(tmp_d_,tmp_param_0in_,0); end;
if (verbose>=1); [~,~,output] = so2g_mle_fminsearch(tmp_d_,tmp_param_0in_,0,tmp_options); end;
output = -output;
if (verbose>0); disp(sprintf(' %% [finished nlrr] gamma_z %0.2f',gamma_z)); end;
