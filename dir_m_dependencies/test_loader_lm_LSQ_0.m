function [beta_E_,beta_I_,zeta_E_,zeta_I_] = test_loader_lm_LSQ_0(C_rank_,E_,I_,dir_trunk,infix,flag_load);

flag_E = 1; flag_I = 1;

[n_u,n_CCOV] = size(C_rank_);
[n_u,n_E_GENE] = size(E_);
[n_u,n_I_GENE] = size(I_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Fit linear model to log of data (xxxx) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

fname_beta = sprintf('%s/dir_mat/beta_%s_X_LSQ_.mat',dir_trunk,infix);
fname_zeta = sprintf('%s/dir_mat/zeta_%s_X_LSQ_.mat',dir_trunk,infix);
if (nargin<6); if (exist(fname_beta,'file') & exist(fname_zeta,'file')); flag_load=1; else flag_load=0; end; end;

if flag_load==0;

%%%%%%%%;
% Here we try and estimate the genes E_ or I_ using the covariates C_. ;
% e.g., : E_ = C_rank_ * zeta_E_. ;
%%%%%%%%;
if flag_E;
zeta_E_ = [ones(n_u,1) , C_rank_]\E_;
end;%if flag_E;
if flag_I;
zeta_I_ = [ones(n_u,1) , C_rank_]\I_;
end;%if flag_I;
%%%%%%%%;
flag_plot=1;
if flag_plot;
colormap(colormap_beach());
tmp_lim_ = mean(E_(:)) + std(E_(:))*1.5*[-1,+1];
subplot(2,2,1); imagesc( E_ , tmp_lim_); title(sprintf('E_%s_',infix),'Interpreter','none');
subplot(2,2,2); imagesc( [ones(n_u,1) , C_rank_]*zeta_E_ , tmp_lim_); title(sprintf('C_rank_*zeta_%s_E_',infix),'Interpreter','none');
subplot(2,2,3); imagesc( I_ , tmp_lim_); title(sprintf('I_%s_',infix),'Interpreter','none');
subplot(2,2,4); imagesc( [ones(n_u,1) , C_rank_]*zeta_I_ , tmp_lim_); title(sprintf('C_rank_*zeta_%s_I_',infix),'Interpreter','none');
figbig;
fname_base = sprintf('%s/dir_jpg/X_vs_C_zeta_%s_LSQ_',dir_trunk,infix);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;

%%%%%%%%;
% Here we try and estimate the covariates C_ using the genes E_ or I_. ;
% e.g., E_ * beta_E_ = C_rank_. ;
%%%%%%%%;
if flag_E;
beta_E_ = [ones(n_u,1) , E_]\C_rank_;
end;%if flag_E;
if flag_I;
beta_I_ = [ones(n_u,1) , I_]\C_rank_;
end;%if flag_I;
%%%%%%%%;
flag_plot=1;
if flag_plot;
colormap(colormap_beach());
subplot(1,3,1); imagesc(C_rank_,[1,n_u]); title('C_rank_','Interpreter','none');
subplot(1,3,2); imagesc( [ones(n_u,1) , E_]*beta_E_ , [1,n_u]); title(sprintf('E_%s_*beta_',infix),'Interpreter','none');
subplot(1,3,3); imagesc( [ones(n_u,1) , I_]*beta_I_ , [1,n_u]); title(sprintf('I_%s_*beta_',infix),'Interpreter','none');
figbig;
fname_base = sprintf('%s/dir_jpg/C_vs_X_beta_%s_LSQ_',dir_trunk,infix);
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot;

save(fname_beta,'beta_E_','beta_I_');
save(fname_zeta,'zeta_E_','zeta_I_');

end;%if flag_load==0;

if flag_load==1;

load(fname_beta,'beta_E_','beta_I_');
load(fname_zeta,'zeta_E_','zeta_I_');

%%%%%%%%;
if flag_E;
tmp_ = zeta_E_; fname_tsv = sprintf('%s/dir_mat/E_%s_zeta_LSQ_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = beta_E_; fname_tsv = sprintf('%s/dir_mat/E_%s_beta_LSQ_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
end;%if flag_E;
%%%%%%%%;
if flag_I;
tmp_ = zeta_I_; fname_tsv = sprintf('%s/dir_mat/I_%s_zeta_LSQ_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
tmp_ = beta_I_; fname_tsv = sprintf('%s/dir_mat/I_%s_beta_LSQ_.tsv',dir_trunk,infix); save(fname_tsv,'tmp_','-ascii','-tabs');
end;%if flag_I;
%%%%%%%%;

end;%if flag_load==1;
