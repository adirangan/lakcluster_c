function s = surfline_0(x_,y_,z_,c_);
n_n = numel(x_);
if nargin<4; c_ = z_; z_ = zeros(size(x_)); end;
x_ = reshape(x_,[1,n_n]);
y_ = reshape(y_,[1,n_n]);
z_ = reshape(z_,[1,n_n]);
c_ = reshape(c_,[1,n_n]);
s = surface([x_;x_],[y_;y_],[z_;z_],[c_;c_],...
        'FaceColor','none',...
        'EdgeColor','interp',...
        'LineWidth',2);
