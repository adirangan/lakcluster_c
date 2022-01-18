function D_ = test_loader_helper_CLabel_Distance_0(dir_trunk,E_,prefix,CLabel_);
% Generates distance matrix D_, where D_(na,nb) = average distance between unique CLabels na and nb. ;

[n_u,N_E] = size(E_);

[u_CLabel_,~,tmp_ij_] = unique(CLabel_);
n_CLabel = length(u_CLabel_);
n_CLabel_ = zeros(n_CLabel,1);
for nCLabel=1:n_CLabel;
n_CLabel_(nCLabel) = length(find(CLabel_==u_CLabel_(nCLabel)));
end;%for nCLabel=1:n_CLabel;

z_ = zeros(n_u,n_CLabel);
na=0;
for nCLabel=1:n_CLabel;
tmp_l = n_CLabel_(nCLabel);
z_(na+(1:tmp_l),nCLabel) = 1.0./tmp_l;
na = na+tmp_l;
end;%for nCLabel=1:n_CLabel;
assert(na==n_u);

tmp_E_ = E_(tmp_ij_,:);
tmp_S_ = sum(tmp_E_.^2,2);
tmp_D_ = ones(n_u,1)*transpose(tmp_S_) + tmp_S_*ones(1,n_u) - 2*tmp_E_*transpose(tmp_E_);
D_ = transpose(z_) * tmp_D_ * z_;
tmp_0in = mean(D_(find(eye(n_CLabel)==1)));
tmp_out = mean(D_(find(eye(n_CLabel)==0)));
tmp_v = (tmp_out - tmp_0in)/tmp_0in;
disp(sprintf(' %% %s: D_ %0.4f - %0.4f = %0.4f --> %0.4f',prefix,tmp_out,tmp_0in,tmp_out-tmp_0in,tmp_v));

clf;imagesc(D_,[0,max(D_(:))]); colormap(colormap_beach());
set(gca,'XTick',[],'YTick',[]);
title(sprintf('%s: %0.4f',prefix,tmp_v),'Interpreter','none');
figbig;
fname_jpg = sprintf('%s/dir_jpg/%sD.jpg',dir_trunk,prefix);
disp(sprintf(' %% writing %s',fname_jpg));
print('-djpeg',fname_jpg);
fname_eps = sprintf('%s/dir_jpg/%sD.eps',dir_trunk,prefix);
disp(sprintf(' %% writing %s',fname_eps));
print('-depsc',fname_eps);



