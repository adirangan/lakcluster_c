function imagedisp_0(input_,lim_,n_x,n_y);
% displays an image. uses limits lim_. ;
% try: ;
%{
  x_ = linspace(0,2*pi,128);
  y_ = linspace(0,2*pi,128*2);
  A_ = sin(transpose(x_))*cos(y_);
  imagedisp_0(A_,[-1,+1]);
 %}

if nargin<2; lim_ = mean(input_(:)) + 1.5*std(input_(:))*[-1,+1]; end;
if nargin<3; n_x = min(size(input_,1),36*1.25); end;
if nargin<4; n_y = min(size(input_,2),72*1.25); end;
[n_r,n_c] = size(input_);
r_sum_ = sparse(max(1,min(n_x,round(linspace(1,n_x,n_r)))),1:n_r,1,n_x,n_r);
tmp_ = max(1,sum(r_sum_,2));
r_sum_ = sparse(1:n_x,1:n_x,1./tmp_,n_x,n_x)*r_sum_;
c_sum_ = sparse(max(1,min(n_y,round(linspace(1,n_y,n_c)))),1:n_c,1,n_y,n_c);
tmp_ = max(1,sum(c_sum_,2));
c_sum_ = sparse(1:n_y,1:n_y,1./tmp_,n_y,n_y)*c_sum_;
input_avg_ = r_sum_*input_*transpose(c_sum_);
%symbol_ = [32,183,164,215,43,35];
%symbol_ = [32,183,43,35];
symbol_ = [32,183,'x',35];
n_symbol = length(symbol_);
input_tab_ = max(1,min(n_symbol,1 + floor( n_symbol * (input_avg_ - lim_(1))/diff(lim_) )));
for nx=1:n_x;
disp(char(symbol_(input_tab_(nx,:))));
end;%for nx=1:n_x;


