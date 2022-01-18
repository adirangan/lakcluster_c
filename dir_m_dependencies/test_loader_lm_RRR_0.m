function [beta_E_un_,beta_E_vn_,n_beta_E,beta_I_un_,beta_I_vn_,n_beta_I,zeta_E_un_,zeta_E_vn_,n_zeta_E,zeta_I_un_,zeta_I_vn_,n_zeta_I] = test_loader_lm_RRR_0(C_rank_,E_,I_,dir_trunk,infix,flag_load);
%%%%%%%%;
% Fit linear model linking covariates to data: ;
% E_ = [ones(n_u,1) , C_rank_] * zeta_E_un_(:,1:n_rank)*transpose(zeta_E_vn_(:,1:n_rank));
% I_ = [ones(n_u,1) , C_rank_] * zeta_I_un_(:,1:n_rank)*transpose(zeta_I_vn_(:,1:n_rank));
%%%%%%%%;
% Fit linear model linking data to covariates: ;
% C_rank_ = [ones(n_u,1) , E_] * beta_E_un_(:,1:n_rank)*transpose(beta_E_vn_(:,1:n_rank));
% C_rank_ = [ones(n_u,1) , I_] * beta_I_un_(:,1:n_rank)*transpose(beta_I_vn_(:,1:n_rank));
%%%%%%%%;

flag_E = 1; flag_I = 1;

[n_u,n_CCOV] = size(C_rank_);
[n_u,n_E_GENE] = size(E_);
[n_u,n_I_GENE] = size(I_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Fit linear model to log of data (xxxx) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

fname_beta = sprintf('%s/dir_mat/beta_%s_X_RRR_.mat',dir_trunk,infix);
fname_zeta = sprintf('%s/dir_mat/zeta_%s_X_RRR_.mat',dir_trunk,infix);
if (nargin<6); if (exist(fname_beta,'file') & exist(fname_zeta,'file')); flag_load=1; else flag_load=0; end; end;

if flag_load==0;

%%%%%%%%;
% Here we try and estimate the genes E_ or I_ using the covariates C_. ;
% For now we use all of C_ to estimate E_ and I_. ;
% Eventually, we might want to use only the C_ which are not associated with quality. ;
%%%%%%%%;
n_rank = 32;
%%%%%%%%;
% Can we estimate the link E_ = C_rank_ * zeta_E_  using a low-rank zeta_E_ ?
%%%%%%%%;
if flag_E;
zeta_E_un_ = zeros(1+n_CCOV,n_rank);
zeta_E_vn_ = zeros(n_E_GENE,n_rank);
tmp_Bn_ = pinv([ones(n_u,1) , C_rank_]);
tmp_Yn_ = E_;
for nrank=1:n_rank;
[tmp_un_,tmp_vn_] = rrr_0([ones(n_u,1) , C_rank_],tmp_Yn_,tmp_Bn_);
zeta_E_un_(:,nrank) = tmp_un_;
zeta_E_vn_(:,nrank) = tmp_vn_;
tmp_Yn_ = tmp_Yn_ - (ones(n_u,1)*zeta_E_un_(1,nrank) + C_rank_*zeta_E_un_(2:end,nrank))*transpose(zeta_E_vn_(:,nrank));
clear tmp_un_ tmp_vn_;
end;%for nrank=1:n_rank;
clear tmp_Yn_ tmp_Bn_;
end;%if flag_E;
%%%%%%%%
% Can we estimate the link I_ = C_rank_ * zeta_I_ using a low-rank zeta_I_ ?
%%%%%%%%
if flag_I;
zeta_I_un_ = zeros(1+n_CCOV,n_rank);
zeta_I_vn_ = zeros(n_I_GENE,n_rank);
tmp_Bn_ = pinv([ones(n_u,1) , C_rank_]);
tmp_Yn_ = I_;
for nrank=1:n_rank;
[tmp_un_,tmp_vn_] = rrr_0([ones(n_u,1) , C_rank_],tmp_Yn_,tmp_Bn_);
zeta_I_un_(:,nrank) = tmp_un_;
zeta_I_vn_(:,nrank) = tmp_vn_;
tmp_Yn_ = tmp_Yn_ - (ones(n_u,1)*zeta_I_un_(1,nrank) + C_rank_*zeta_I_un_(2:end,nrank))*transpose(zeta_I_vn_(:,nrank));
clear tmp_un_ tmp_vn_;
end;%for nrank=1:n_rank;
clear tmp_Yn_ tmp_Bn_;
end;%if flag_I;
%%%%%%%%;
% Appears to be rank 3 or so ;
%%%%%%%%;
flag_disp=1;
if flag_disp;
if flag_E;
subplot(2,2,1);
bar(1:n_rank,sqrt(sum((ones(n_u,1)*zeta_E_un_(1,:) + C_rank_*zeta_E_un_(2:end,:)).^2,1)),'r');
xlabel('pc #'); ylabel('l2-norm');
title(sprintf('numerical rank of zeta_%s_E_',infix),'Interpreter','none');
subplot(2,2,3);
bar(1:n_rank,log10(sqrt(sum((ones(n_u,1)*zeta_E_un_(1,:) + C_rank_*zeta_E_un_(2:end,:)).^2,1))),'r');
xlabel('pc #'); ylabel('log10(l2-norm)');
title(sprintf('numerical rank of zeta_%s_E_',infix),'Interpreter','none');
end;%if flag_E;
if flag_I;
subplot(2,2,2);
bar(1:n_rank,sqrt(sum((ones(n_u,1)*zeta_I_un_(1,:) + C_rank_*zeta_I_un_(2:end,:)).^2,1)),'b');
xlabel('pc #'); ylabel('l2-norm');
title(sprintf('numerical rank of zeta_%s_I_',infix),'Interpreter','none');
subplot(2,2,4);
bar(1:n_rank,log10(sqrt(sum((ones(n_u,1)*zeta_I_un_(1,:) + C_rank_*zeta_I_un_(2:end,:)).^2,1))),'b');
xlabel('pc #'); ylabel('log10(l2-norm)');
title(sprintf('numerical rank of zeta_%s_I_',infix),'Interpreter','none');
end;%if flag_I;
set(gcf,'Position',1+[0,0,512*2,512]);
fname_base = sprintf('%s/dir_jpg/zeta_%s_X_RRR_',dir_trunk,infix);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_disp;
%%%%%%%%;
n_zeta_E = 3;
n_zeta_I = 3;

flag_plot=1;
if flag_plot;
colormap(colormap_beach());
tmp_lim_ = mean(E_(:)) + std(E_(:))*1.5*[-1,+1];
subplot(2,2,1); imagesc( E_ , tmp_lim_); title(sprintf('E_%s_',infix),'Interpreter','none');
subplot(2,2,2); imagesc( ( ones(n_u,1)*zeta_E_un_(1,:) + C_rank_*zeta_E_un_(2:end,:) ) * transpose(zeta_E_vn_(:,:)) , tmp_lim_); title(sprintf('C_rank_*zeta_%s_E_',infix),'Interpreter','none');
subplot(2,2,3); imagesc( I_ , tmp_lim_); title(sprintf('I_%s_',infix),'Interpreter','none');
subplot(2,2,4); imagesc( ( ones(n_u,1)*zeta_I_un_(1,:) + C_rank_*zeta_I_un_(2:end,:) ) * transpose(zeta_I_vn_(:,:)) , tmp_lim_); title(sprintf('C_rank_*zeta_%s_I_',infix),'Interpreter','none');
figbig;
fname_base = sprintf('%s/dir_jpg/X_vs_C_zeta_%s_RRR_',dir_trunk,infix);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;

%%%%%%%%;
% Here we try and estimate the covariates C_ using the genes E_ or I_. ;
% For now we estimate all of C_. ;
% Eventually, we might want to only estimate the C_ which are not associated with a quality. ;
%%%%%%%%;
n_rank = 32;
%%%%%%%%;
% Can we estimate the link E_ * beta_E_ = C_rank_ using a low-rank beta_E_ ?
%%%%%%%%;
if flag_E;
beta_E_un_ = zeros(1+n_E_GENE,n_rank);
beta_E_vn_ = zeros(n_CCOV,n_rank);
tmp_Bn_ = pinv([ones(n_u,1) , E_]);
tmp_Yn_ = C_rank_;
for nrank=1:n_rank;
[tmp_un_,tmp_vn_] = rrr_0([ones(n_u,1) , E_],tmp_Yn_,tmp_Bn_);
beta_E_un_(:,nrank) = tmp_un_;
beta_E_vn_(:,nrank) = tmp_vn_;
tmp_Yn_ = tmp_Yn_ - (ones(n_u,1)*beta_E_un_(1,nrank) + E_*beta_E_un_(2:end,nrank))*transpose(beta_E_vn_(:,nrank));
clear tmp_un_ tmp_vn_;
end;%for nrank=1:n_rank;
clear tmp_Yn_ tmp_Bn_;
end;%if flag_E;
%%%%%%%%
% Can we estimate the link I_ * beta_I_ = C_rank_ using a low-rank beta_I_ ?
%%%%%%%%
if flag_I;
beta_I_un_ = zeros(1+n_I_GENE,n_rank);
beta_I_vn_ = zeros(n_CCOV,n_rank);
tmp_Bn_ = pinv([ones(n_u,1) , I_]);
tmp_Yn_ = C_rank_;
for nrank=1:n_rank;
[tmp_un_,tmp_vn_] = rrr_0([ones(n_u,1) , I_],tmp_Yn_,tmp_Bn_);
beta_I_un_(:,nrank) = tmp_un_;
beta_I_vn_(:,nrank) = tmp_vn_;
tmp_Yn_ = tmp_Yn_ - (ones(n_u,1)*beta_I_un_(1,nrank) + I_*beta_I_un_(2:end,nrank))*transpose(beta_I_vn_(:,nrank));
clear tmp_un_ tmp_vn_;
end;%for nrank=1:n_rank;
clear tmp_Yn_ tmp_Bn_;
end;%if flag_I;
%%%%%%%%;
% Appears to be rank 6 or so ;
%%%%%%%%;
flag_disp=1;
if flag_disp;
if flag_E;
subplot(2,2,1);
bar(1:n_rank,sqrt(sum((ones(n_u,1)*beta_E_un_(1,:) + E_*beta_E_un_(2:end,:)).^2,1)),'r');
xlabel('pc #'); ylabel('l2-norm');
title(sprintf('numerical rank of beta_%s_E_',infix),'Interpreter','none');
subplot(2,2,3);
bar(1:n_rank,log10(sqrt(sum((ones(n_u,1)*beta_E_un_(1,:) + E_*beta_E_un_(2:end,:)).^2,1))),'r');
xlabel('pc #'); ylabel('log10(l2-norm)');
title(sprintf('numerical rank of beta_%s_E_',infix),'Interpreter','none');
end;%if flag_E;
if flag_I;
subplot(2,2,2);
bar(1:n_rank,sqrt(sum((ones(n_u,1)*beta_I_un_(1,:) + I_*beta_I_un_(2:end,:)).^2,1)),'b');
xlabel('pc #'); ylabel('l2-norm');
title(sprintf('numerical rank of beta_%s_I_',infix),'Interpreter','none');
subplot(2,2,4);
bar(1:n_rank,log10(sqrt(sum((ones(n_u,1)*beta_I_un_(1,:) + I_*beta_I_un_(2:end,:)).^2,1))),'b');
xlabel('pc #'); ylabel('log10(l2-norm)');
title(sprintf('numerical rank of beta_%s_I_',infix),'Interpreter','none');
end;%if flag_I;
set(gcf,'Position',1+[0,0,512*2,512]);
fname_base = sprintf('%s/dir_jpg/beta_%s_X_RRR_',dir_trunk,infix);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_disp;
%%%%%%%%;
n_beta_E = 6;
n_beta_I = 6;

flag_plot=1;
if flag_plot;
colormap(colormap_beach());
subplot(1,3,1); imagesc(C_rank_,[1,n_u]); title('C_rank_','Interpreter','none');
subplot(1,3,2); imagesc( ( ones(n_u,1)*beta_E_un_(1,:) + E_*beta_E_un_(2:end,:) ) * transpose(beta_E_vn_(:,:)) , [1,n_u]); title(sprintf('E_%s_*beta_',infix),'Interpreter','none');
subplot(1,3,3); imagesc( ( ones(n_u,1)*beta_I_un_(1,:) + I_*beta_I_un_(2:end,:) ) * transpose(beta_I_vn_(:,:)) , [1,n_u]); title(sprintf('I_%s_*beta_',infix),'Interpreter','none');
figbig;
fname_base = sprintf('%s/dir_jpg/C_vs_X_beta_%s_RRR_',dir_trunk,infix);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;

save(fname_beta,'beta_E_un_','beta_E_vn_','n_beta_E','beta_I_un_','beta_I_vn_','n_beta_I');
save(fname_zeta,'zeta_E_un_','zeta_E_vn_','n_zeta_E','zeta_I_un_','zeta_I_vn_','n_zeta_I');

end;%if flag_load==0;

if flag_load==1;

load(fname_beta,'beta_E_un_','beta_E_vn_','n_beta_E','beta_I_un_','beta_I_vn_','n_beta_I');
load(fname_zeta,'zeta_E_un_','zeta_E_vn_','n_zeta_E','zeta_I_un_','zeta_I_vn_','n_zeta_I');

end;%if flag_load==1;

%%%%%%%%;
if flag_E;
tmp_ = zeta_E_un_; fname_tsv = sprintf('%s/dir_mat/E_%s_zeta_un_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = zeta_E_vn_; fname_tsv = sprintf('%s/dir_mat/E_%s_zeta_vn_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = beta_E_un_; fname_tsv = sprintf('%s/dir_mat/E_%s_beta_un_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = beta_E_vn_; fname_tsv = sprintf('%s/dir_mat/E_%s_beta_vn_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
end;%if flag_E;
%%%%%%%%;
if flag_I;
tmp_ = zeta_I_un_; fname_tsv = sprintf('%s/dir_mat/I_%s_zeta_un_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = zeta_I_vn_; fname_tsv = sprintf('%s/dir_mat/I_%s_zeta_vn_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = beta_I_un_; fname_tsv = sprintf('%s/dir_mat/I_%s_beta_un_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = beta_I_vn_; fname_tsv = sprintf('%s/dir_mat/I_%s_beta_vn_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
end;%if flag_I;
%%%%%%%%;
