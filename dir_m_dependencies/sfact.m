function output = sfact(d_);
%%%%%%%%;
if (nargin<1);
d_ = 25:-1:0;
%ln_1_ = log(factorial(d_));
ln_1_ = gammaln(1+d_);
ln_2_ = sfact(d_);
subplot(1,2,1);
plot(d_,ln_1_,'k.-',d_,ln_2_,'ro-');
xlabel('d');ylabel('log(d!)');
subplot(1,2,2);
plot(d_,log10(abs(ln_1_-ln_2_)./ln_1_),'kx-');
xlabel('d');ylabel('log10(relative error in log(d!))');
disp('returning');return;
end;%if (nargin<1);
%%%%%%%%;
%{
tmp_ij_ = find(d_>10);
%tmp_ij_ = find(d_>=0);
d1_ij_ = d_(tmp_ij_);
e1_ij_ = 1./d1_ij_;
e2_ij_ = e1_ij_.^2;
e3_ij_ = e1_ij_.*e2_ij_;
e5_ij_ = e3_ij_.*e2_ij_;
e7_ij_ = e5_ij_.*e2_ij_;
output(tmp_ij_) = ...
  + d1_ij_.*log(d1_ij_) ...
  - d1_ij_ ...
  + 0.5.*log(2*pi*d1_ij_) ...
  + 0.0833333333333333.*e1_ij_ ... %<-- 1/12;
  - 0.0027777777777778.*e3_ij_ ... %<-- 1/360;
  + 0.0007936507936508.*e5_ij_ ... %<-- 1/1260;
  - 0.0005952380952381.*e7_ij_ ; %<-- 1/1680;
%tmp_ = [0.0000000000000000 0.0000000000000000 0.6931471805599453 1.7917594692280550 3.1780538303479458 4.7874917427820458 6.5792512120101012 8.5251613610654147 10.6046029027452509 12.8018274800814691 15.1044125730755159]; %<-- log(factorial([0:10]));
%tmp_ij_ = find(d_<=10 & d_>=0);
%output(tmp_ij_) = tmp_(1+d_(tmp_ij_));
tmp_ij_ = find(d_<=10 & d_>=0);
output(tmp_ij_) = gammaln(1+d_(tmp_ij_));
tmp_ij_ = find(d_<0);
output(tmp_ij_) = -Inf;
 %}
output = zeros(size(d_));
tmp_ij_ = find(d_>=0);
output(tmp_ij_) = gammaln(1+d_(tmp_ij_));
tmp_ij_ = find(d_<0);
output(tmp_ij_) = -Inf;
