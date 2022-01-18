function output = fnet_helper_0(z_,M__);
% sets up energy equation for positions at z_ given interaction matrix M__ ;
% We expect the force to be something like: ;
% F_{i} = \sum_{j} n_{ij}*[1/d_{ij} - M_{ij}] ;
% where n_{ij} = unit vector from j to i ;
% d_{ij} = distance from j to i ;
% M_{ij} = interaction between j and i ;
% The energy associated with this force would be something like: ;
% E_{i} = \sum_{j} 1/d_{ij}^2 + M_{ij}d_{ij}^2 ;
% test with: ;
%{

  n_pt = 8; M__ = rand(n_pt); 
  %rij = 1:floor(n_pt/2); M__(rij,rij) = M__(rij,rij) + 1.3;
  %rij = floor(n_pt/2) + (1:floor(n_pt/2)); M__(rij,rij) = M__(rij,rij) + 1.8;
  rij = 1:2:n_pt; M__(rij,rij) = M__(rij,rij) + 1.3;
  rij = 2:2:n_pt; M__(rij,rij) = M__(rij,rij) + 1.8;
  M__ = (M__+transpose(M__))/2; for npt1=1:n_pt; M__(npt1,npt1)=0; end;
  z_ = zeros(n_pt,2); z_(:,1) = cos(2*pi*(0:n_pt-1)/n_pt); z_(:,2) = sin(2*pi*(0:n_pt-1)/n_pt);
  iteration_max = 16; output_ = zeros(2*n_pt,iteration_max); fval_ = zeros(iteration_max,1);
  for niteration = 1:iteration_max;
  z_tmp_ = z_(randperm(n_pt),:);
  [output,fval,fsolve_exitflag,fsolve_output] = fminsearch(@(z_) fnet_helper_0(z_,M__),z_tmp_(:),optimset('Display','off','Jacobian','off'));
  output_(:,niteration) = output; fval_(niteration) = norm(fval);
  end;%for niteration = 1:iteration_max;
  [tmp,tmp_ij] = min(fval_); M_z__(:,1) = output_(1:n_pt,tmp_ij); M_z__(:,2) = output_(n_pt+(1:n_pt),tmp_ij);
  cra = colormap(jet);
  cla; 
  subplot(1,2,1); imagesc(M__); axis equal;
  subplot(1,2,2); hold on;
  M_thresh = 0;
  for npt1 = 1:n_pt;
  for npt2 = 1:n_pt;
  if (M__(npt1,npt2)>M_thresh); 
  l = line([M_z__(npt1,1);M_z__(npt2,1)],[M_z__(npt1,2);M_z__(npt2,2)]);
  set(l,'LineWidth',M__(npt1,npt2)-M_thresh,'Color','k');
  end;%if (M__(npt1,npt2)>M_thresh); 
  end;%for npt2 = 1:n_pt;
  end;%for npt1 = 1:n_pt;
  for npt1 = 1:n_pt;
  cb = max(1,min(size(cra,1),floor(size(cra,1)*(npt1+0.5)/n_pt)));
  plot(M_z__(npt1,1),M_z__(npt1,2),'.','markersize',50,'Color',cra(cb,:)); 
  end;%for npt1 = 1:n_pt;
  hold off;
  axis equal;

  %}

n_pt = length(z_)/2;
x_ = reshape(z_(1:n_pt),n_pt,1);
y_ = reshape(z_(n_pt+1:end),n_pt,1);
dx__ = repmat(x_,1,n_pt) - repmat(transpose(x_),n_pt,1);
dy__ = repmat(y_,1,n_pt) - repmat(transpose(y_),n_pt,1);
n2__ = max(0.0001,dx__.^2 + dy__.^2) ;
gamma = 2.0 ; beta = 2 ;
term1__ = 1./n2__.^gamma; term2__ = (1+beta*M__).*n2__.^2;
term1__ = term1__ - diag(diag(term1__)); term2__ = term2__ - diag(diag(term2__));
E_ = sum(term1__,2) + sum(term2__,2);
output = sum(E_);
