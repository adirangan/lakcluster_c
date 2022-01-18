function dolphin_ellipse_from_DD_inv_0(DD_inv__,c_,linewidth_use,zradius);
% plots 2d-ellipse from DD_inv__ = dt*BB_inv__+CC_inv__. ;
% Assumes DD_inv__ is 2-by-2. ;

na=0;
if (nargin<1+na); DD_inv__=[]; end; na=na+1;
if (nargin<1+na); c_=[]; end; na=na+1;
if (nargin<1+na); linewidth_use=[]; end; na=na+1;
if (nargin<1+na); zradius=[]; end; na=na+1;

if isempty(c_); c_ = [0,0,0]; end;
if isempty(linewidth_use); linewidth_use = 2; end;
if isempty(zradius); zradius = 1; end;

D_inv__ = sqrtm(DD_inv__);
[tmp_U__,tmp_D__] = svds(D_inv__,2);
n_psi = 64;
psi_ = linspace(0,2*pi,n_psi+1);
x_ = cos(psi_);
y_ = sin(psi_);
xy__ = [ x_ ; y_ ];
D_inv_xy__ = D_inv__*xy__;
hold on;
p=patch(zradius*transpose(D_inv_xy__(1,:)),zradius*transpose(D_inv_xy__(2,:)),-1e-3*zradius*ones(1+n_psi,1),c_);
set(p,'EdgeColor','none');
%p=patch(1.0*transpose(D_inv_xy__(1,:)),1.0*transpose(D_inv_xy__(2,:)),-1e-6*ones(1+n_psi,1),c_.^1.00);
%p=patch(2.0*transpose(D_inv_xy__(1,:)),2.0*transpose(D_inv_xy__(2,:)),-1e-3*ones(1+n_psi,1),c_.^0.50);
%set(p,'EdgeColor','none');
%p=patch(1.0*transpose(D_inv_xy__(1,:)),1.0*transpose(D_inv_xy__(2,:)),-1e-6*ones(1+n_psi,1),c_.^1.00);
%set(p,'EdgeColor','none');
%plot([0;+tmp_D__(1,1)*tmp_U__(1,1)],[0;+tmp_D__(1,1)*tmp_U__(2,1)],'-','Color',c_,'LineWidth',linewidth_use);
%plot([0;-tmp_D__(1,1)*tmp_U__(1,1)],[0;-tmp_D__(1,1)*tmp_U__(2,1)],'-','Color',c_,'LineWidth',linewidth_use);
%plot([0;+tmp_D__(2,2)*tmp_U__(1,2)],[0;+tmp_D__(2,2)*tmp_U__(2,2)],'-','Color',c_,'LineWidth',linewidth_use);
%plot([0;-tmp_D__(2,2)*tmp_U__(1,2)],[0;-tmp_D__(2,2)*tmp_U__(2,2)],'-','Color',c_,'LineWidth',linewidth_use);
%plot(D_inv_xy__(1,:),D_inv_xy__(2,:),'-','Color',c_,'LineWidth',linewidth_use);
hold off;



