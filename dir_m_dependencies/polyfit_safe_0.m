function [p_,s,mu_] = polyfit_safe_0(x_,y_,n);

if nargin<1;
%%%%%%%%;
% check polyfit. ;
%%%%%%%%;
rng(0);
tmp_n = 1024*1;
tmp_p_ = [0.3;1;2];
tmp_x_ = randn(tmp_n,1);
tmp_y_ = polyval(tmp_p_,tmp_x_) + 0.010*randn(tmp_n,1);
tmp_x0 = sum(tmp_x_.^0);
tmp_x1 = sum(tmp_x_.^1);
tmp_x2 = sum(tmp_x_.^2);
tmp_x3 = sum(tmp_x_.^3);
tmp_x4 = sum(tmp_x_.^4);
tmp_x0y = sum(tmp_x_.^0.*tmp_y_);
tmp_x1y = sum(tmp_x_.^1.*tmp_y_);
tmp_x2y = sum(tmp_x_.^2.*tmp_y_);
tmp_A__ = [ ...
     tmp_x4 , tmp_x3 , tmp_x2 ; ...
     tmp_x3 , tmp_x2 , tmp_x1 ; ...
     tmp_x2 , tmp_x1 , tmp_x0 ; ...
	    ];
tmp_b_ = [ ...
      tmp_x2y ; ...
      tmp_x1y ; ...
      tmp_x0y ; ...
	   ];
tmp_q_ = tmp_A__\tmp_b_;
tmp_r_ = polyfit_safe_0(tmp_x_,tmp_y_,2); tmp_r_ = tmp_r_(:);
disp(sprintf(' %% tmp_p_ vs tmp_q_: %0.16f',fnorm(tmp_p_-tmp_q_)/fnorm(tmp_p_)));
disp(sprintf(' %% tmp_p_ vs tmp_r_: %0.16f',fnorm(tmp_p_-tmp_r_)/fnorm(tmp_p_)));
disp(sprintf(' %% tmp_q_ vs tmp_r_: %0.16f',fnorm(tmp_q_-tmp_r_)/fnorm(tmp_q_)));
flag_plot=1;
if flag_plot;
figure(1);clf;figmed;
subplot(1,2,1);
hold on;
plot(tmp_x_,tmp_y_,'k.','MarkerSize',16);
plot(sort(tmp_x_),polyval(tmp_p_,sort(tmp_x_)),'r-','LineWidth',3);
plot(sort(tmp_x_),polyval(tmp_q_,sort(tmp_x_)),'g-','LineWidth',3);
n_degree = 16;
[ ...
 c_obj ...
,x_use_ ...
,y_use_ ...
,bic_ ...      
,c_ ...
] = ...
chebfit_safe_1( ...
 tmp_x_ ...
,tmp_y_ ...
,n_degree ...
);
plot(sort(tmp_x_),c_obj(sort(tmp_x_)),'m-','LineWidth',2);
hold off;
xlabel('x');
ylabel('y');
subplot(1,2,2);
hold on;
plot(0:n_degree,bic_,'ro');
hold off;
xlabel('degree');
ylabel('BIC');
end;%if flag_plot;
disp('returning'); return;
end;%if nargin<1;

index_use_ = intersect(efind(isfinite(x_)),efind(isfinite(y_)));
n_use = min(n,numel(index_use_)-1);
if nargout<3;
[p_,s] = polyfit(x_(1+index_use_),y_(1+index_use_),n_use);
end;%if nargout<3;
if nargout==3;
[p_,s,mu_] = polyfit(x_(1+index_use_),y_(1+index_use_),n_use);
end;%if nargout==3;
