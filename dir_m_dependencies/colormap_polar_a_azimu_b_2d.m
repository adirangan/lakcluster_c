function c_ = colormap_polar_a_azimu_b_2d(polar_a_,azimu_b_,gamma);
if nargin<2;
n_x = 129; dx1 = pi/(n_x-1); dx2 = dx1/2;
polar_a_ = linspace(0,1*pi,1*n_x);
azimu_b_ = linspace(0,2*pi,2*n_x);
[Azimu_b__,Polar_a__] = ndgrid(azimu_b_,polar_a_);
c_ = colormap_polar_a_azimu_b_2d(Polar_a__(:),Azimu_b__(:));
c__ = reshape(c_,[1,2*n_x^2,3]);
azimu_b__ = 0*pi+transpose(repmat(Azimu_b__(:),[1,4]) + dx2*repmat([-1,+1,+1,-1],[2*n_x^2,1]));
polar_a__ = 1*pi-transpose(repmat(Polar_a__(:),[1,4]) + dx2*repmat([-1,-1,+1,+1],[2*n_x^2,1]));
patch(azimu_b__,polar_a__,c__,'EdgeColor','none');
axisnotick;
axis image;
figbig;
disp('returning');return;
end;%if nargin<2;

if nargin<3; gamma = 0.35; end;
polar_a_ = polar_a_(:);
azimu_b_ = azimu_b_(:);
c_ = colorpsphere_8points_gamma(polar_a_,azimu_b_,gamma);

  
