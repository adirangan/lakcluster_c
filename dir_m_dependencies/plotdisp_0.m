function plotdisp_0(x_,y_,xlim_,ylim_,n_x,n_y);
% displays an image. uses limits lim_. ;
% try: ;
%{
  x_ = linspace(0,2*pi,128);
  y_ = sin(x_).*x_;
  plotdisp_0(x_,y_);
 %}

ni = 2;
if nargin<ni; y_ = x_; x_ = 1:length(y_(:)); end; ni=ni+1;
if nargin<ni; xlim_ = [min(x_(:)),max(x_(:))]; end; ni=ni+1;
if nargin<ni; ylim_ = [min(y_(:)),max(y_(:))]; end; ni=ni+1;
if nargin<ni; n_x = min(length(x_),72*1.25); end; ni=ni+1;
if nargin<ni; n_y = min(length(y_),36*1.25); end; ni=ni+1;
nx_ = max(1,min(n_x,1+floor(n_x*(x_-xlim_(1))/diff(xlim_))));
ny_ = max(1,min(n_y,1+floor(n_y*(y_-ylim_(1))/diff(ylim_))));
input_ = sparse(1+n_y-ny_,nx_,1,n_y,n_x);
imagedisp_0(input_,[0,1],n_y,n_x);


