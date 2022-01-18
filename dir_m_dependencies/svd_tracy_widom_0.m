function [svd_up,lambda_tw,avg_tw,std_tw] = svd_tracy_widom_0(n_row,n_col,p_val);
% Returns the lambda_tw above which an eigenvalue is more significant than p_val. ;
%%%%%%%%;
% based on: ;
% Folkmar Bornemann table, ;
% as well as:  ;
% https://www.sciencedirect.com/science/article/pii/S0047259X14000761?via%3Dihub ;
% Journal of Multivariate Analysis, Volume 129, August 2014, Pages 69-81, ;
% Distribution of the largest eigenvalue for real Wishart and Gaussian random matrices and a simple approximation for the Tracy–Widom distribution ;
% Marco Chiani ;
%%%%%%%%;
% gammainc(x,a)*gamma(a) = \int_{0}^{x} t^(a-1) * exp(-t) dt ;
% test with: ;
% x=2;a=4; disp(integral(@(t) t.^(a-1).*exp(-t),0,x) - gammainc(x,a)*gamma(a));
%%%%%%%%;
% test out pdf and cdf of tracy-widom distribution. ;
%{
n_row = 127; n_col = 1781;
gamma_tw = n_row/n_col;
n_sample = 1024*1;
T_svd_ = zeros(n_sample,1);
a1 = -0.5; a2 = -0.5; %<-- beta = 1 for real. ;
avg_tw = ( sqrt(n_row + a1) + sqrt(n_col + a2) ).^2 ;
std_tw = sqrt(avg_tw) * ( 1/sqrt(n_row + a1) + 1/sqrt(n_col + a2) ).^(1/3);
for nsample=0:n_sample-1;
if (mod(nsample,128)==0); disp(sprintf('%% nsample %d/%d',nsample,n_sample)); end;
A_ = randn(n_row,n_col); %<-- note real so that beta = 1;
% H_ = A_*transpose(A_); H_svd = svds(H_,1); %<-- note that H_svd = A_svd^2. ;
A_svd = svds(A_,1); T_svd_(1+nsample) = (A_svd^2 - avg_tw)/std_tw;
end;%for nsample=0:n_sample-1;
%%%%%%%%;
tw_bornemann = load('TW_beta1.mat'); %<-- bornemann. ;
n_h = 32;
h_x_ = linspace(min(tw_bornemann.x),max(tw_bornemann.x),n_h);
h_ = hist(T_svd_,h_x_); h_ = h_./sum(h_)/mean(diff(h_x_));
figure(1);clf;
subplot(1,2,1);
hold on;
stairs(h_x_,h_,'r-');
plot(tw_bornemann.x,tw_bornemann.TW_s_tag,'k-');
plot(tw_bornemann.x,tw_chiani_pdf(tw_bornemann.x),'g-');
hold off;
xlim([min(tw_bornemann.x),max(tw_bornemann.x)]);
xlabel('x'); ylabel('#');
title('pdf');
subplot(1,2,2);
hold on;
stairs(h_x_,cumsum(h_)/sum(h_),'r-');
plot(tw_bornemann.x,tw_bornemann.TW_s,'k-');
plot(tw_bornemann.x,tw_chiani_cdf(tw_bornemann.x),'g-');
hold off;
xlabel('x'); ylabel('#');
title('cdf');
figbig;
 %}

if (nargin<1);
disp(sprintf(' %% testing svd_tracy_widom_0'));
n_row = 127; n_col = 1781;
[svd_up,lambda_tw,avg_tw,std_tw] = svd_tracy_widom_0(n_row,n_col,0.01); disp(sprintf(' %% p = 0.01: %.6f %.6f %.6f %.6f',svd_up,lambda_tw,avg_tw,std_tw));
svd_up_01 = svd_up;
[svd_up,lambda_tw,avg_tw,std_tw] = svd_tracy_widom_0(n_row,n_col,0.05); disp(sprintf(' %% p = 0.05: %.6f %.6f %.6f %.6f',svd_up,lambda_tw,avg_tw,std_tw));
svd_up_05 = svd_up;
[svd_up,lambda_tw,avg_tw,std_tw] = svd_tracy_widom_0(n_row,n_col,0.15); disp(sprintf(' %% p = 0.15: %.6f %.6f %.6f %.6f',svd_up,lambda_tw,avg_tw,std_tw));
svd_up_15 = svd_up;
n_sample = 128;
svd_ = zeros(n_sample,1);
b = 0.010;
for nsample=0:n_sample-1; svd_(1+nsample) = svds(b*randn(n_row,n_col),1); end;
figure(1);clf;hold on;
plot((1:n_sample)/n_sample,sort(svd_),'k-');
plot((1:n_sample)/n_sample,ones(1,n_sample)*svd_up_01*b,'r-','LineWidth',2.0);
plot((1:n_sample)/n_sample,ones(1,n_sample)*svd_up_05*b,'r-','LineWidth',1.0);
plot((1:n_sample)/n_sample,ones(1,n_sample)*svd_up_15*b,'r-','LineWidth',0.5);
grid on;
hold off;
xlim([0,1]);
title('cdf');
disp('returning');return;
end;%if (nargin<1);

a1 = -0.5; a2 = -0.5; %<-- beta = 1 for real. ;
avg_tw = ( sqrt(n_row + a1) + sqrt(n_col + a2) ).^2 ;
std_tw = sqrt(avg_tw) * ( 1/sqrt(n_row + a1) + 1/sqrt(n_col + a2) ).^(1/3);
lambda_tw = fzero(@(x) tw_chiani_cdf(x) - (1-p_val),0);
svd_up = sqrt(avg_tw + std_tw*lambda_tw) ;







