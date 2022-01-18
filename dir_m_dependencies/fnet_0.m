function z_out_ = fnet_0(M__,rseed,z_0in_);
% runs energy minimizer (see fnet_helper_0.m) ;
% test with: ;
%{

  n_pt = 8; M__ = rand(n_pt); 
  %rij = 1:floor(n_pt/2); M__(rij,rij) = M__(rij,rij) + 1.3;
  %rij = floor(n_pt/2) + (1:floor(n_pt/2)); M__(rij,rij) = M__(rij,rij) + 1.8;
  rij = 1:2:n_pt; M__(rij,rij) = M__(rij,rij) + 1.3;
  rij = 2:2:n_pt; M__(rij,rij) = M__(rij,rij) + 1.8;
  M__ = (M__+transpose(M__))/2; for npt1=1:n_pt; M__(npt1,npt1)=0; end;
  z_ = zeros(n_pt,2); z_(:,1) = cos(2*pi*(0:n_pt-1)/n_pt); z_(:,2) = sin(2*pi*(0:n_pt-1)/n_pt);
  z_out_ = fnet_0(z_,M__);
  cra = colormap(jet);
  cla; 
  subplot(1,2,1); imagesc(M__); axis equal;
  subplot(1,2,2); hold on;
  M_thresh = 0;
  for npt1 = 1:n_pt;
  for npt2 = 1:n_pt;
  if (M__(npt1,npt2)>M_thresh); 
  l = line([z_out_(npt1,1);z_out_(npt2,1)],[z_out_(npt1,2);z_out_(npt2,2)]);
  set(l,'LineWidth',M__(npt1,npt2)-M_thresh,'Color','k');
  end;%if (M__(npt1,npt2)>M_thresh); 
  end;%for npt2 = 1:n_pt;
  end;%for npt1 = 1:n_pt;
  for npt1 = 1:n_pt;
  cb = max(1,min(size(cra,1),floor(size(cra,1)*(npt1+0.5)/n_pt)));
  plot(z_out_(npt1,1),z_out_(npt1,2),'.','markersize',50,'Color',cra(cb,:)); 
  end;%for npt1 = 1:n_pt;
  hold off;
  axis equal;

  %}
if (nargin<2); rseed=[]; end;
if (isempty(rseed)); rseed = 0; end;
if (nargin<3); z_0in_=[]; end;
if (isempty(z_0in_)); z_0in_ = 2*pi*(transpose(0:size(M__,1)-1))/size(M__,1); z_0in_ = [cos(z_0in_),sin(z_0in_)]; end;
rng(rseed);
n_z = size(z_0in_,1); for nz=1:n_z; M__(nz,nz)=0; end;
n_iteration = 16; output_ = zeros(2*n_z,n_iteration); fval_ = zeros(n_iteration,1);
for niteration = 1:n_iteration;
if niteration==1; tmp_z_ = z_0in_(randperm(n_z),:); else tmp_z_ = pre_z_(randperm(n_z),:); end;
[output,fval,fsolve_exitflag,fsolve_output] = fminsearch(@(z_0in_) fnet_helper_0(z_0in_,M__),tmp_z_(:),optimset('Display','off','Jacobian','off'));
pre_z_ = zeros(n_z,2); pre_z_(:,1) = output(1:n_z); pre_z_(:,2) = output(n_z+(1:n_z));
output_(:,niteration) = output; fval_(niteration) = norm(fval);
end;%for niteration = 1:n_iteration;
[tmp,tmp_ij] = min(fval_); z_out_(:,1) = output_(1:n_z,tmp_ij); z_out_(:,2) = output_(n_z+(1:n_z),tmp_ij);
