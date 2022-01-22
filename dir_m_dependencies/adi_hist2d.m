function output = adi_hist2d(xra,yra,xl,nxbins,yl,nybins);
% sets up a 2d histogram;
if (nargin<3);
xl=[min(xra),max(xra)]; xl = mean(xl) + diff(xl)/2*1.0625*[-1,1];
yl=[min(yra),max(yra)]; yl = mean(yl) + diff(yl)/2*1.0625*[-1,1];
nxbins = 32;
nybins = 32;
end;%if (nargin<3);

xra = xra(:);
yra = yra(:);
len = length(xra);
bx = 1+min(nxbins-1,max(0,floor(nxbins*(xra-min(xl))/(max(xl)-min(xl)))));
by = 1+min(nybins-1,max(0,floor(nybins*(yra-min(yl))/(max(yl)-min(yl)))));
output = sparse(by,bx,ones(len,1),nybins,nxbins);
