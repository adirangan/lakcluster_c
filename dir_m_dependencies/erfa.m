function output = erfa(input);
% approximates erfc(input). ;
% This is quite accurate (errors < 1e-6) ;
% when the input is > 3 or so. ;
if nargin<1;
x_=linspace(-5,5,1024);
tmp_erfc_ = erfc(x_);
tmp_erfa_ = erfa(x_);
subplot(1,3,1); hold on;
plot(x_,(tmp_erfc_),'k-',x_,(tmp_erfa_),'r-'); xlabel('x'); ylabel('erfc');
subplot(1,3,2); hold on;
plot(x_,log(tmp_erfc_),'k-',x_,log(tmp_erfa_),'r-'); xlabel('x'); ylabel('log(erfc)');
subplot(1,3,3); hold on;
plot(x_,log(tmp_erfc_)-log(tmp_erfa_),'k-'); xlabel('x'); ylabel('log(erfc)-log(erfa)');
disp('returning');return;
end;%if nargin<1;

%%%%%%%%;

output = zeros(size(input));
tmp_ij_ = find(input>=0);
x1_ = input(tmp_ij_);
x2_ = x1_.*x1_;
x3_ = x1_.*x2_;
x5_ = x3_.*x2_;
x7_ = x5_.*x2_;
x9_ = x7_.*x2_;
erft_ = erfc(0) + 2/sqrt(pi) * ( ...
	  + (1/factorial(1))*(    -1)*x1_ ...
	  + (1/factorial(3))*(    +2)*x3_ ...
	  + (1/factorial(5))*(   -12)*x5_ ...
	  + (1/factorial(7))*(  +120)*x7_ ...
	  + (1/factorial(9))*( -1680)*x9_ ...
	  );
erff_ = exp(-x2_) .* x1_ ./sqrt(pi) ...
  ./ ( x2_ + 0.5 ...
       ./ (1 + 1 ...
	   ./ (x2_ + 1.5 ...
	       ./ ( 1 + 2.0 ...
		    ./ ( x2_ + 2.5 ...
			 ) ) ) ) ) ;
erfa_ = max(erft_,erff_);
output(tmp_ij_) = erfa_;
tmp_ij_ = find(input<0);
if (length(tmp_ij_)>0); output(tmp_ij_) = 2-erfa(-input(tmp_ij_)); end;
