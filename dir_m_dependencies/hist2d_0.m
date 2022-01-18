function [output_,xl_,yl_] = hist2d_0(x_,y_,nxbins,nybins,xl_,yl_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% function output_ = hist2d_0(x_,y_,nxbins,nybins,xl_,yl_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% 
% Sets up a 2d histogram;
%
% Inputs: 
% x_ double array of length N storing x-coordinates of data. ;
% y_ double array of length N storing y-coordinates of data. ;
% xl_ double array of length 2 storing min and max values for x-coordinates. ;
% yl_ double array of length 2 storing min and max values for y-coordinates. ;
% nxbins integer number of bins in x-direction. ;
% nybins integer number of bins in y-direction. ;
%
% Outputs:
% output_ matrix of size nybins-x-nxbins storing 2d-histogram of data. ;
% Note that output_ is in sparse format. ;
%
% Test with: hist2d_0();
% 

if (nargin<1);
G1_ = 0.15*randn(1024*4,2) + repmat([-0.15,-0.25],1024*4,1);
G2_ = 0.15*randn(1024*1,2) + repmat([+0.15,+0.25],1024*1,1);
G_ = [G1_;G2_];
colormap('hot');
imagesc(hist2d_0(G_(:,1),G_(:,2),32,33,[-1,+1],[-1,+1])); colorbar;
disp('returning');return;
end;%if (nargin<1);

if (nargin<3);
nxbins = 32; nybins = 32;
end;%if (nargin<3);

if (nargin<5);
xl_=[min(x_),max(x_)]; xl_ = mean(xl_) + diff(xl_)/2*1.0625*[-1,1];
yl_=[min(y_),max(y_)]; yl_ = mean(yl_) + diff(yl_)/2*1.0625*[-1,1];
end;%if (nargin<5);

x_ = x_(:);
y_ = y_(:);
len = length(x_);
bx = 1+min(nxbins-1,max(0,floor(nxbins*(x_-min(xl_))/(max(xl_)-min(xl_)))));
by = 1+min(nybins-1,max(0,floor(nybins*(y_-min(yl_))/(max(yl_)-min(yl_)))));
output_ = sparse(by,bx,ones(len,1),nybins,nxbins);
