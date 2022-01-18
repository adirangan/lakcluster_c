function DDinv_tvv__ = chebout_dolphin_estimate_DD_0(dt_,BB__,CC__,eig_tolerance);
if (nargin<4); eig_tolerance = 1e-6; end;
n_var = size(BB__,1);
DDinv_tvv__ = zeros(numel(dt_),n_var^2);
for ndt=0:numel(dt_)-1;
dt = dt_(1+ndt);
DDinv_tvv__(1+ndt,:) = reshape(dolphin_estimate_DD_0(dt,BB__,CC__,eig_tolerance),[1,n_var^2]);
end;%for ndt=0:numel(dt_)-1;
