function [C_,D_,niteration] = svd_impute_fit_3(B_,ij_missed_,n_rank,str_prefix,rseed,n_iteration,tolerance);
% Assuming simple svd-based imputation. ;
% initializes based on median of each column. ;
% C_ is imputed, D_ is imputed and fitted. ;

if nargin<1;
n_r = 128*4; n_c = 10*n_r+1; n_dim = floor(n_r/4); s_ = transpose(exp(-linspace(0,1,n_dim)/0.0625));
A_ = randn(n_r,n_c);
[tmp_U_,tmp_S_,tmp_V_] = svds(A_,n_dim); tmp_S_ = diag(diag(tmp_S_).*s_);
A_ = tmp_U_*tmp_S_*transpose(tmp_V_);
ij_missed_ = randperm(numel(A_)); ij_missed_ = ij_missed_(1:floor(numel(A_)/6));
ij_filled_ = setdiff(1:numel(A_),ij_missed_);
B_ = A_; B_(ij_missed_) = 0;
[C_,D_,niteration] = svd_impute_fit_3(B_,ij_missed_,4);
figure(1);clf;
colormap(colormap_pm());
l_ = mean(A_(:)) + 1.5*std(A_(:))*[-1,+1];
subplot(2,3,1); imagesc(A_,l_); set(gca,'XTick',[],'YTick',[]); title('original');
subplot(2,3,4); imagesc(B_,l_); set(gca,'XTick',[],'YTick',[]); title('perforated');
subplot(2,3,2); imagesc(C_,l_); set(gca,'XTick',[],'YTick',[]); title('imputed');
tmp_C_ = C_; tmp_C_(ij_filled_)=0;
subplot(2,3,5); imagesc(tmp_C_,l_); set(gca,'XTick',[],'YTick',[]); title('imputed (missing only)');
subplot(2,3,3); imagesc(D_,l_); set(gca,'XTick',[],'YTick',[]); title('fitted');
tmp_D_ = D_; tmp_D_(ij_filled_)=0;
subplot(2,3,6); imagesc(tmp_D_,l_); set(gca,'XTick',[],'YTick',[]); title('fitted (missing only)');
figbig;
figure(2);clf;
hold on;
plot(1:n_dim,svds(A_,n_dim),'ko-');
plot(1:n_dim,svds(B_,n_dim),'bo-');
plot(1:n_dim,svds(C_,n_dim),'go-');
plot(1:n_dim,svds(D_,n_dim),'ro-');
legend({'original','perforated','imputed','fitted'},'Location','NorthEast');
disp('returning'); return;
end;%if nargin<1;

na=2;
if nargin<na; ij_missed_ = find(B_==0); end; na=na+1;
if nargin<na; n_rank = 1; end; na=na+1;
if nargin<na; str_prefix = 'test'; end; na=na+1;
if nargin<na; rseed = 1; end; na=na+1;
if nargin<na; n_iteration = 256; end; na=na+1;
if nargin<na; tolerance = 1e-6; end; na=na+1;
% n_rank = 1; str_prefix = 'test'; rseed = 1; n_iteration = 256; tolerance = 1e-6; 
rng(rseed);

verbose_flag=1;

if (verbose_flag); disp(sprintf(' %% [entering svd_impute_fit_3] str_prefix %s',str_prefix)); end;
p_missed = numel(ij_missed_)/numel(B_);
ij_filled_ = setdiff(1:numel(B_),ij_missed_);
p_filled = numel(ij_filled_)/numel(B_);
e_ = zeros(n_iteration,1);
[n_row,n_col] = size(B_); n_dim = min(n_row,n_col);
C_ = impute_randperm(B_,ij_missed_);
sC_median_ = svd(impute_median(B_,ij_missed_));
sC_randperm_ = svd(impute_randperm(B_,ij_missed_));

if (verbose_flag); disp(sprintf(' %% initial imputation')); end;
%%%%%%%%;
% svd imputation. ;
%%%%%%%%;
continue_flag=1;
niteration=0;
while continue_flag;
niteration = niteration+1;
if (verbose_flag); if (mod(niteration,10)==0); disp(sprintf(' %% niteration %d/%d',niteration,n_iteration)); end; end;
[tmp_U_,tmp_S_,tmp_V_] = svds(C_,min(min(size(C_)),n_rank));
tmp_B_ = tmp_U_*tmp_S_*transpose(tmp_V_);
tmp_val_ = C_(ij_missed_);
C_(ij_missed_) = tmp_B_(ij_missed_);
e_(niteration) = fnorm(C_(ij_missed_)-tmp_val_)./fnorm(C_(ij_missed_));
if e_(niteration) < tolerance; continue_flag=0; end;
if niteration>=n_iteration; continue_flag=0; end;
end;%for niteration=1:n_iteration;
clear tmp_U_ tmp_S_ tmp_V_ tmp_B_ ;

if (verbose_flag & niteration>1 & n_rank>0);
figure;clf;
subplot(1,2,1);
plot(1:niteration,e_(1:niteration),'k.-');
xlabel('iteration');ylabel('error');title('error in svd imputation');
subplot(1,2,2);
plot(1:niteration,log10(e_(1:niteration)),'k.-');
xlabel('iteration');ylabel('log10(error)');title('log10(error) in svd imputation');
disp(sprintf(' %% writing %s_FIGA.jpg',str_prefix));
print('-djpeg',sprintf('%s_FIGA.jpg',str_prefix));
print('-depsc',sprintf('%s_FIGA.eps',str_prefix));
end;%if (verbose_flag);

sC_impute_ = svd(C_);
if (verbose_flag); disp(sprintf(' %% fitting spectrum')); end;
%%%%%%%%;
% fit spectrum with noise. ;
%%%%%%%%;
[tmp_U_,tmp_S_,tmp_V_] = svd(C_);
[tmp_Q_,~] = qr(randn(size(C_,1)-n_rank));
tmp_UQ_ = orth(tmp_U_(:,n_rank+1:end)*tmp_Q_);
tmp_Q_ = randn(n_col,n_dim-n_rank);
tmp_VQ_ = orth(tmp_Q_ - tmp_V_(:,1:n_rank) * ( transpose(tmp_V_(:,1:n_rank)) * tmp_Q_ ));
%tmp_UQ_ = orth(randn(n_row,n_dim-n_rank));
%tmp_VQ_ = orth(randn(n_col,n_dim-n_rank));
%%%%%%%%;
%[tmp_M_,~] = qr(randn(size(C_,1))); %<-- random rotation. ;
[tmp_M_] = sparse(1:n_row,randperm(n_row),1); %<-- random permutation. ;
tmp_S_scale_ = diag(sC_randperm_(1:n_rank));
tmp_C_top_ = tmp_U_(:,1:n_rank)*tmp_S_scale_*transpose(tmp_V_(:,1:n_rank));
if (verbose_flag); disp(sprintf(' %% [%d,%d] , [%d,%d] , [%d,%d]',size(tmp_UQ_),size(tmp_S_(n_rank+1:n_dim,n_rank+1:n_dim)/sqrt(p_filled)),size(transpose(tmp_VQ_)))); end;
tmp_C_bot_ = tmp_UQ_ * tmp_S_(n_rank+1:n_dim,n_rank+1:n_dim)/sqrt(p_filled) * transpose(tmp_VQ_);
fh_ = @(f_scale_) fnorm(sC_randperm_ - svd(impute_randperm(tmp_M_*svd_impute_fit_helper_0(f_scale_,C_,tmp_C_top_,tmp_C_bot_,ij_missed_),ij_missed_)));
f_scale_opt_ = fminsearch(fh_,[1;1],optimset('TolFun',tolerance,'TolX',tolerance));
if (verbose_flag); disp(sprintf(' %% f_scale_opt_ [%0.2f %0.2f]',f_scale_opt_)); end;
D_ = svd_impute_fit_helper_0(f_scale_opt_,C_,tmp_C_top_,tmp_C_bot_,ij_missed_);

if (verbose_flag);
figure;clf;
plot(1:n_dim,sC_randperm_,'ko',1:n_dim,svd(impute_randperm(tmp_M_*D_,ij_missed_)),'rx');
xlabel('rank'); ylabel('singular value'); title('fitting perforated spectrum'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthEast');
xlim([1,n_dim]);
disp(sprintf(' %% writing %s_FIGB.jpg',str_prefix));
print('-djpeg',sprintf('%s_FIGB.jpg',str_prefix));
print('-depsc',sprintf('%s_FIGB.eps',str_prefix));
figure;clf;
plot(1:n_dim,log(sC_randperm_),'ko',1:n_dim,log(svd(impute_randperm(tmp_M_*D_,ij_missed_))),'rx');
xlabel('rank'); ylabel('log(singular value)'); title('fitting perforated spectrum'); legend({'perforated as given','imputed + spun +perforated'},'location','NorthEast');
xlim([1,n_dim]);
disp(sprintf(' %% writing %s_FIGC.jpg',str_prefix));
print('-djpeg',sprintf('%s_FIGC.jpg',str_prefix));
print('-depsc',sprintf('%s_FIGC.eps',str_prefix));
end;%if (verbose_flag);


if (verbose_flag); disp(sprintf(' %% [finished svd_impute_fit_3] str_prefix %s',str_prefix)); end;
