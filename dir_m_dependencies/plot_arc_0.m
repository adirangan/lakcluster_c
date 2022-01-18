function l = plot_arc_0(x1,x2,r,plot_flag);
% test with;
%{

  cla;hold on;
  for iteration = 1:16;
  x1 = randn(2,1);
  x2 = randn(2,1);
  r_ = -0.5:0.1:+0.5;
  cra = colormap('jet');
  plot(x1(1),x1(2),'b.');
  plot(x2(1),x2(2),'r.');
  for nr=1:length(r_);
  cb = max(1,min(size(cra,1),floor(size(cra,1)*(nr+0.5)/length(r_))));
  l = plot_arc_0(x1,x2,r_(nr)); set(l,'Color',cra(cb,:));
  end;% for nr=1:length(r_);
  end%  for iteration = 1:16;
  hold off;
  xlim([-2,2]); ylim([-2,2]); axis square;

  %}
if nargin<4; plot_flag=1; end;
x1=x1(:); x2=x2(:);
xm = 0.5*(x1+x2);
%plot(xm(1),xm(2),'g.'); axis equal;
dx = x1-x2; m = norm(dx/2); r = sign(r)*min(0.75*m,abs(r));
xt = [-dx(2);+dx(1)]; xt = xt/norm(xt);
alpha = (m^2 - r^2)/(2*r); d = sqrt(alpha^2+m^2); 
xc = xm + alpha*xt;
%[m, alpha, r, abs(alpha)+abs(r) , d , norm(xc-x1)],;
%plot(xc(1),xc(2),'k.'); axis equal;
x1c = x1-xc; x2c = x2-xc;
theta1c = atan2(x1c(2),x1c(1));theta2c = atan2(x2c(2),x2c(1));
%plot(xc(1)+d*cos(theta1c),xc(2)+d*sin(theta1c),'bo'); axis equal;
%plot(xc(1)+d*cos(theta2c),xc(2)+d*sin(theta2c),'ro'); axis equal;
dtheta = theta2c-theta1c;
if (dtheta>0 & dtheta<pi); theta_first = theta1c; theta_secnd = theta2c; theta_span = dtheta; else; theta_first = theta2c; theta_secnd = theta1c; theta_span = theta1c-theta2c; end;
while (theta_span>2*pi); theta_span = theta_span - 2*pi; end;
while (theta_span<0); theta_span = theta_span + 2*pi; end;
if (theta_span>pi); tmp = theta_first; theta_first = theta_secnd; theta_secnd = theta_first; theta_span = 2*pi - theta_span; end;
npts = 32;
theta_ = linspace(theta_first,theta_first + theta_span,npts); 
x_ = d*cos(theta_) + xc(1);
y_ = d*sin(theta_) + xc(2);
xord_ = [x_(1:end-1) ; x_(2:end)];
yord_ = [y_(1:end-1) ; y_(2:end)];
if plot_flag; l = line(xord_,yord_);
 else; l = {xord_,yord_}; end;

