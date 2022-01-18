function output = igammaln(s,x);
% simple approximation to log(igamma(s,x));
if (nargin<1);
x_ = 0:.5:60; s = 44; plot(x_,log(igamma(s,x_.^2/2)),'.',x_,igammaln(s,x_.^2/2),'o');
xlabel('x'); ylabel('igammaln(44,x)');
disp('returning');return;
end;%if (nargin<1);

%output = log(igamma(s,x));
output = log(gammainc(x,s,'upper')) + gammaln(s);
tmp_index_ = efind(x>0 & ~isfinite(output));
if (~isempty(tmp_index_));
output(1+tmp_index_) = (s-1).*log(x(1+tmp_index_)) - x(1+tmp_index_);
end;%if (~isempty(tmp_index_));



