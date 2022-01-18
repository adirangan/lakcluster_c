function [xdrop_condensed_,Z_avg_,Z_avg_condensed_,xdrop_,xrmv_] = tensor_auc_ver3(A);
% rank-score for tensor A;
% test with: tensor_auc_ver3();

if nargin<1;
%rng(1);
N=96; X = ceil(N^0.55);
A = randn(N,N+1); B=randn(X,X+1)+0.5; 
for nj=1:3; 
pf_{nj} = randperm(size(A,nj)); 
%pf_{nj} = 1:size(A,nj);
[~,pi_{nj}] = sort(pf_{nj}); [~,pi_{nj}] = sort(pi_{nj}); 
end;% for nj=1:3;
A(pf_{1}(1:X),pf_{2}(1:X+1)) = B;
tic; [xdrop_,Z_avg_,~,~,~] = tensor_auc_ver3(A); t_tot = toc; 
disp(sprintf(' time: %f',t_tot));
for nj=1:3; xdrop_{nj} = xdrop_{nj}(end:-1:1); end;
figure;
subplot(1,3,1); imagesc(sum(A(pi_{1},pi_{2},pi_{3}),3));
subplot(1,3,2); imagesc(sum(A(xdrop_{1},xdrop_{2},xdrop_{3}),3));
subplot(1,3,3); hold on; plot(Z_avg_{1});plot(Z_avg_{2});plot(Z_avg_{3});hold off;
for nj=1:3; 
[~,tmp_Bn_{nj}] = intersect(xdrop_{nj},pi_{nj}(1:size(B,nj)),'stable');
[~,tmp_Bc_{nj}] = intersect(xdrop_{nj},pi_{nj}(size(B,nj)+1:end),'stable');
Auc_(nj) = auc_0(tmp_Bn_{nj},tmp_Bc_{nj});
end;%for nj=1:3; 
Auc_avg = mean(Auc_); disp(sprintf(' %% Auc %0.2f',Auc_avg));
return;
end;% if nargin<1;

if nargin<1;
rng(1);
N=96; X = ceil(N^0.65);
A = randn(N,N+1,N+2); B=randn(X,X+1,X+2)+0.95; 
for nj=1:3; 
%pf_{nj} = randperm(size(A,nj)); 
pf_{nj} = 1:size(A,nj);
[~,pi_{nj}] = sort(pf_{nj}); [~,pi_{nj}] = sort(pi_{nj}); 
end;% for nj=1:3;
A(pf_{1}(1:X),pf_{2}(1:X+1),pf_{3}(1:X+2)) = B;
tic; [xdrop_,Z_avg_,~,~,~] = tensor_auc_ver3(A); t_tot = toc; 
disp(sprintf(' time: %f',t_tot));
for nj=1:3; xdrop_{nj} = xdrop_{nj}(end:-1:1); end;
figure;
subplot(1,3,1); imagesc(sum(A(pi_{1},pi_{2},pi_{3}),3));
subplot(1,3,2); imagesc(sum(A(xdrop_{1},xdrop_{2},xdrop_{3}),3));
subplot(1,3,3); hold on; plot(Z_avg_{1});plot(Z_avg_{2});plot(Z_avg_{3});hold off;
for nj=1:3; 
[~,tmp_Bn_{nj}] = intersect(xdrop_{nj},pi_{nj}(1:size(B,nj)),'stable');
[~,tmp_Bc_{nj}] = intersect(xdrop_{nj},pi_{nj}(size(B,nj)+1:end),'stable');
Auc_(nj) = auc_0(tmp_Bn_{nj},tmp_Bc_{nj});
end;%for nj=1:3; 
Auc_avg = mean(Auc_); disp(sprintf(' %% Auc %0.2f',Auc_avg));
return;
end;% if nargin<1;

verbose=0;

[xol_(1),xol_(2),xol_(3)] = size(A);
xtot = sum(xol_);
xrmv_ = zeros(1,xtot);
xtmp_ = xol_;

for nx=1:xtot; 
tmp_ij_ = find(xtmp_>1);
if (length(tmp_ij_)>0);
xrat = xtmp_(tmp_ij_)./(1+xol_(tmp_ij_));
[~,xij] = max(xrat);
xrmv_(nx) = tmp_ij_(xij);
xtmp_(tmp_ij_(xij)) = xtmp_(tmp_ij_(xij))-1; 
end;%if (length(tmp_ij_)>0);
if (length(tmp_ij_)==0);
xrat = xtmp_./(1+xol_);
[~,xij] = max(xrat);
xrmv_(nx) = xij;
xtmp_(xij) = xtmp_(xij)-1; 
end;%if (length(tmp_ij_)==0);
end;%for nx=1:xtot;

%{
xrmv_ = zeros(1,xtot);
xrmv_(1:xol_(1)) = 1;
xrmv_(xol_(1) + (1:xol_(2))) = 2;
xrmv_(end) = 3;
 %}

A_{1} = A;
for n2 = 1:xol_(2); for n3 = 1:xol_(3);
A_tmp = A_{1}(:,n2,n3); [~,tmp_s] = sort(A_tmp,'ascend'); [~,tmp_i] = sort(tmp_s,'ascend'); A_{1}(:,n2,n3) = tmp_i/xol_(1);
end;end;%for n2 = 1:xol_(2); for n3 = 1:xol_(3);

A_{2} = A;
for n1 = 1:xol_(1); for n3 = 1:xol_(3);
A_tmp = A_{2}(n1,:,n3); [~,tmp_s] = sort(A_tmp,'ascend'); [~,tmp_i] = sort(tmp_s,'ascend'); A_{2}(n1,:,n3) = tmp_i/xol_(2);
end;end;%for n1 = 1:xol_(1); for n3 = 1:xol_(3);

A_{3} = A;
for n1 = 1:xol_(1); for n2 = 1:xol_(2);
A_tmp = A_{3}(n1,n2,:); [~,tmp_s] = sort(A_tmp,'ascend'); [~,tmp_i] = sort(tmp_s,'ascend'); A_{3}(n1,n2,:) = tmp_i/xol_(3);
end;end;%for n1 = 1:xol_(1); for n2 = 1:xol_(2);

clear Z_ xdrop_ ;

Z_{1} = sum(squeeze(sum(A_{1},3)),2); Z_{1} = Z_{1}(:);
Z_{2} = sum(squeeze(sum(A_{2},3)),1); Z_{2} = Z_{2}(:);
Z_{3} = sum(squeeze(sum(A_{3},2)),1); Z_{3} = Z_{3}(:);

for nj=1:3; x_rem_{nj} = 1:xol_(nj); end;%for nj=1:3;
Z_denominator_ = zeros(xtot,3);

for nx=1:xtot;

for nj=1:3; 
x_l = zeros(3,1);
x_l(1) = length(x_rem_{1}); x_l(2) = length(x_rem_{2}); x_l(3) = length(x_rem_{3}); x_l(nj) = 1;
Z_denominator_(nx,nj) = x_l(1)*x_l(2)*x_l(3);
tm2 = mean(Z_{nj}(x_rem_{nj})); tmp = tm2/Z_denominator_(nx,nj);
Z_av2_{nj}(nx) = tm2;
if (verbose); disp(sprintf(' %% Z_av2_{%d}(%d) = %f',nj,nx,tm2)); end; %if (verbose);
Z_avg_{nj}(nx) = tmp;
if (verbose); disp(sprintf(' %% Z_avg_{%d}(%d) = %f',nj,nx,tmp)); end; %if (verbose);
Z_avg_condensed_{nj}(1+xol_(nj)-length(x_rem_{nj})) = tmp;
if (verbose); disp(sprintf(' %% Z_avg_condensed_{%d}(%d) = %f',nj,1+xol_(nj)-length(x_rem_{nj}),tmp)); end; %if (verbose);
end;%for nj=1:3; 

if xrmv_(nx)==1;
%[~,ij_rmv] = min(Z_{1}(x_rem_{1})); 
[~,ij_rmv] = min(abs(Z_{1}(x_rem_{1})/Z_denominator_(nx,1)-Z_av2_{1}(1)/Z_denominator_(1,1))); 
xdrop_{1}(nx) = x_rem_{1}(ij_rmv); xdrop_{2}(nx) = 0; xdrop_{3}(nx) = 0;
tmp = squeeze(sum(A_{2}(xdrop_{1}(nx),:,:),3)); Z_{2} = Z_{2} - tmp(:);
tmp = squeeze(sum(A_{3}(xdrop_{1}(nx),:,:),2)); Z_{3} = Z_{3} - tmp(:);
A_{1}(xdrop_{1}(nx),:,:) = 0; A_{2}(xdrop_{1}(nx),:,:) = 0; A_{3}(xdrop_{1}(nx),:,:) = 0;
x_rem_{1} = setdiff(x_rem_{1},xdrop_{1}(nx));
end;%if xrmv_(nx)==1;

if xrmv_(nx)==2;
%[~,ij_rmv] = min(Z_{2}(x_rem_{2})); 
[~,ij_rmv] = min(abs(Z_{2}(x_rem_{2})/Z_denominator_(nx,2)-Z_av2_{2}(1)/Z_denominator_(1,2))); 
xdrop_{2}(nx) = x_rem_{2}(ij_rmv); xdrop_{1}(nx) = 0; xdrop_{3}(nx) = 0;
tmp = squeeze(sum(A_{1}(:,xdrop_{2}(nx),:),3)); Z_{1} = Z_{1} - tmp(:);
tmp = squeeze(sum(A_{3}(:,xdrop_{2}(nx),:),1)); Z_{3} = Z_{3} - tmp(:);
A_{1}(:,xdrop_{2}(nx),:) = 0; A_{2}(:,xdrop_{2}(nx),:) = 0; A_{3}(:,xdrop_{2}(nx),:) = 0;
x_rem_{2} = setdiff(x_rem_{2},xdrop_{2}(nx));
end;%if xrmv_(nx)==2;

if xrmv_(nx)==3;
%[~,ij_rmv] = min(Z_{3}(x_rem_{3})); 
[~,ij_rmv] = min(abs(Z_{3}(x_rem_{3})/Z_denominator_(nx,3)-Z_av2_{3}(1)/Z_denominator_(1,3))); 
xdrop_{3}(nx) = x_rem_{3}(ij_rmv); xdrop_{1}(nx) = 0; xdrop_{2}(nx) = 0;
tmp = squeeze(sum(A_{1}(:,:,xdrop_{3}(nx)),2)); Z_{1} = Z_{1} - tmp(:);
tmp = squeeze(sum(A_{2}(:,:,xdrop_{3}(nx)),1)); Z_{2} = Z_{2} - tmp(:);
A_{1}(:,:,xdrop_{3}(nx)) = 0; A_{2}(:,:,xdrop_{3}(nx)) = 0; A_{3}(:,:,xdrop_{3}(nx)) = 0;
x_rem_{3} = setdiff(x_rem_{3},xdrop_{3}(nx));
end;%if xrmv_(nx)==3;

end;%for nx=1:xtot;

if (verbose);
%disp(transpose(Z_{1}/(length(x_rem_{2})*length(x_rem_{3}))));
%disp(transpose(Z_{2}/(length(x_rem_{1})*length(x_rem_{3}))));
%disp(transpose(Z_{3}/(length(x_rem_{1})*length(x_rem_{2}))));
disp(transpose(Z_{1}(x_rem_{1})/(length(x_rem_{2})*length(x_rem_{3}))));
disp(transpose(Z_{2}(x_rem_{2})/(length(x_rem_{1})*length(x_rem_{3}))));
disp(transpose(Z_{3}(x_rem_{3})/(length(x_rem_{1})*length(x_rem_{2}))));
end;%if (verbose);

for nj=1:3; xdrop_condensed_{nj} = xdrop_{nj}(find(xdrop_{nj}>0)); end;%for nj=1:3;
