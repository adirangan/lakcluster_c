function ...
[ ...
 parameter...
] = ...
vector_field_0( ...
 parameter ...
,A__ ...
);
% This function plots a simple vector field for a 2-dimensional ODE ;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); A__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'x0_lim_'); parameter.x0_lim_ = [-1,+1]; end;
if ~isfield(parameter,'x1_lim_'); parameter.x1_lim_ = [-1,+1]; end;
if ~isfield(parameter,'n_bin'); parameter.n_bin = 16; end;
if ~isfield(parameter,'markersize_use'); parameter.markersize_use = 4; end;
if ~isfield(parameter,'linewidth_use'); parameter.linewidth_use = 2; end;
x0_lim_ = parameter.x0_lim_;
x1_lim_ = parameter.x1_lim_;
n_bin = parameter.n_bin;
markersize_use = parameter.markersize_use;
linewidth_use = parameter.linewidth_use;

x0_val_ = linspace(x0_lim_(1),x0_lim_(2),n_bin); x0_val_ = reshape(x0_val_,1,length(x0_val_));
x1_val_ = linspace(x1_lim_(1),x1_lim_(2),n_bin); x1_val_ = reshape(x1_val_,1,length(x1_val_));
x_ = zeros(2,length(x0_val_)*length(x1_val_));
rhs_ = zeros(2,length(x0_val_)*length(x1_val_));
tmp_ = repmat(x0_val_,length(x1_val_),1);
tmp_=reshape(tmp_,1,length(x0_val_)*length(x1_val_));
x_(1,:) = tmp_;
tmp_ = repmat(x1_val_,1,length(x0_val_));
x_(2,:) = tmp_;
clear tmp_;
rhs_ = A__*x_;
dx0 = diff(x0_lim_)/n_bin;
dx1 = diff(x1_lim_)/n_bin;
theta_ = atan2(rhs_(2,:)/dx1,rhs_(1,:)/dx0);
cos_ = cos(theta_); sin_ = sin(theta_);
v_ = sqrt((rhs_(2,:)/dx1).^2 + (rhs_(1,:)/dx0).^2);
v_mean = mean(v_); v_std = std(v_); v_lim = v_mean+[-1,+1]*1.5*v_std; v_lim(1)=max(0,v_lim(1)); v_lim(2) = min(max(v_),v_lim(2));
c_map__ = colormap(1-colormap('gray'));
n_c_map_bin = length(c_map__);
c_map_nbin_ = max(0,min(n_c_map_bin-1,floor(n_c_map_bin*(v_-v_lim(1))/diff(v_lim))));
c_map_load__ = c_map__(1+c_map_nbin_,:);
dl = 0.45;
lx__ = [x_(1,:) + dl*dx0*cos_ ; x_(1,:) - dl*dx0*cos_];
ly__ = [x_(2,:) + dl*dx1*sin_ ; x_(2,:) - dl*dx1*sin_];
for nl=0:size(lx__,2)-1;
hold on;
l = line(lx__(:,1+nl),ly__(:,1+nl),'Color',c_map_load__(1+nl,:));
set(l,'LineWidth',linewidth_use);
plot( ...
 lx__(1,1+nl) ...
,ly__(1,1+nl) ...
,'o' ...
,'MarkerSize',markersize_use ...
,'MarkerEdgeColor','none' ...
,'MarkerFaceColor',c_map_load__(1+nl,:) ...
);
hold off;
end;%for nl=0:size(lx__,2)-1;
