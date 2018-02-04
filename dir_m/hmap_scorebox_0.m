function [area_out,ik_out,ij_out] = hmap_scorebox_0(LR2_ra,Lx2_ra,rthr,ar_min,ar_max,plot_flag,LR2_lim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [area_out,ik_out,ij_out] = hmap_scorebox_0(LR2_ra,Lx2_ra,rthr,ar_min,ar_max,plot_flag,LR2_lim);
%
% This function finds and returns the corner-indices (ij_out,ik_out) of the array LR2_ra
% that satisfy the following constraints: 
% 1. the product area_out = ik_out*ij_out should be as large as possible (i.e., maximum area) 
% 2. the entry Lx2_ra(ij_out,ik_out) should be at a maximum (i.e., the full number of covariate-categories)
% 3. the aspect ratio (ij_out/nrows)/(ik_out/ncols) should lie in between ar_min and ar_max
% 4. the value LR2_ra(ij_out,ik_out) should be at least a threshold-value LR2_thr.
% Within the above constraints, the threshold-value LR2_thr is determined using the relative threshold rthr (an input).
% specifically, LR2_thr = LR2_min + rthr*(LR2_max - LR2_min)
% where LR2_min and LR2_max are either: 
% (i) the smallest and largest values of LR2_ra, respectively, or
% (ii) the smallest and largest values of LR2_lim, if LR2_lim has two different values.
% The inputs to this function are:
% LR2_ra: an array of values (e.g., average row-scores) to threshold in constraint #4.
% Lx2_ra: an array of values (e.g., number of covariate-categories remaining) to use in constraint #2.
% rthr: a real number between 0 and 1 to use in constraint #4.
% ar_min, ar_max; positive real numbers to use in constraint #3.
% plot_flag: an integer (0 or 1); if this is set to 1 the area_out is plotted.
% LR2_lim: forced maximum and minimum for LR2_min and LR2_max (optional).
%
% test by running with no arguments:
% i.e., >> hmap_scorebox_0();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<1);
disp(sprintf(' '));
disp(sprintf(' testing hmap_scorebox_0: '));
disp(sprintf(' First we set up a simple array ''LR2_ra'' with large values in the upper left corner (to simulate biclustered data). '));
LR2_ra = sort(rand(156,1),'descend')*sort(rand(1,212),'descend');
disp(sprintf(' Then we set up another array ''Lx2_ra'' of the same size, with various integer values ranging from 1-3. '));
Lx2_ra = [1*ones(2,1);2*ones(1,1);3*ones(153,1)]*ones(1,212);
disp(sprintf(' This second array will simulate the number of remaining covariates associated with each pixel of LR2_ra. '));
disp(sprintf(' Then we pick a relative threshold (in this case rthr = 0.7).'));
rthr = 0.7;
disp(sprintf(' Then we pick a minimum and maximum acceptable aspect-ratio (in this case ar_min = 0.5, ar_max = 1.0/ar_min). '));
ar_min = 0.5; ar_max = 1.0/ar_min;
disp(sprintf(' Finally, we use hmap_scorebox_0 to search for the corner submatrix of LR2_ra which has the largest '));
disp(sprintf(' possible area while remaining larger than the relative threshold (and maintaining an acceptable aspect-ratio): '));
disp(sprintf(' >> hmap_scorebox_0(LR2_ra,Lx2_ra,rthr,ar_min,ar_max,1); '));
hmap_scorebox_0(LR2_ra,Lx2_ra,rthr,ar_min,ar_max,1);
disp(sprintf(' '));
disp(sprintf(' This output is plotted in the figure. '));
disp(sprintf(' The white dots indicate pixels that correspond to corner-submatrices above the relative threshold with '));
disp(sprintf(' an acceptable aspect-ratio. The black circle-x shows the maximum area point. '));
return;
end;%if (nargin<1);

if (nargin<6); plot_flag=0; end;

[nrows,ncols] = size(LR2_ra);
area_out = zeros(1,1);
D = LR2_ra; 

[ij,ik] = find(Lx2_ra==max(Lx2_ra(:)));
if (nargin<7 | length(unique(LR2_lim))<2);
LR2_max = 0;
for nl=1:length(ij);
LR2_max = max(LR2_max,LR2_ra(ij(nl),ik(nl)));
LR2_min = max(0,min(LR2_max,LR2_ra(ij(nl),ik(nl))));
end;%for nl=1:length(ij);
end;%if (nargin<7 | length(unique(LR2_lim))<2);
if (nargin>=7 & length(unique(LR2_lim))>=2);
LR2_min = min(LR2_lim(:)); LR2_max = max(LR2_lim(:));
end;%if (nargin>=7 & length(unique(LR2_lim))>=2);
LR2_thr = LR2_min + rthr*(LR2_max - LR2_min);
[ij,ik] = find(LR2_ra>LR2_thr & Lx2_ra==max(Lx2_ra(:)));
loc_valid = find((ij/nrows)./(ik/ncols) > ar_min & (ij/nrows)./(ik/ncols) < ar_max);
if isempty(loc_valid); loc_valid=[1]; ij=1;ik=1; end;
area_out=0; if (length(loc_valid)>0); [area_out,loc] = max(ij(loc_valid).*ik(loc_valid)); end;%if (length(loc_valid)>0);

if plot_flag;
figure;cla; colormap('hot'); % colormap('jet');
if (length(loc_valid)>0);
hold on;
imagesc(LR2_ra(end:-1:1,:),[LR2_min,LR2_max]); 
if (nargin<7 | length(unique(LR2_lim))<2); colorbar; end;
if plot_flag>2;
for nl=1:length(loc_valid);
plot(ik(loc_valid(nl)),size(LR2_ra,1)-ij(loc_valid(nl)),'w.');
end;%for nl=1:length(loc_valid);
end;%if plot_flag>2;
if plot_flag==1;
plot(ik(loc_valid(loc)),size(LR2_ra,1)-ij(loc_valid(loc)),'kx','Markersize',10);
end;%if plot_flag==1;
if plot_flag>1;
plot(ik(loc_valid(loc)),size(LR2_ra,1)-ij(loc_valid(loc)),'ko','Markersize',35);
plot(ik(loc_valid(loc)),size(LR2_ra,1)-ij(loc_valid(loc)),'kx','Markersize',35);
end;%if plot_flag>1;
hold off;
xlim([1,size(LR2_ra,2)]); xlabel('columns n');
ylim([1,size(LR2_ra,1)]); ylabel('rows m');
set(gca,'Xtick',[ik(loc_valid(loc))],'Xticklabel',[ik(loc_valid(loc))],'Ytick',[size(LR2_ra,1)-ij(loc_valid(loc))],'Yticklabel',[ij(loc_valid(loc))]);
if plot_flag>1;
title(sprintf('valid areas in white; maximum area at ''x'': relative threshold %0.2f, abs threshold %0.2f',rthr,LR2_thr));
end; %if plot_flag>1;
if plot_flag==1;
title(sprintf('maximum area at ''x'': relative threshold %0.2f, abs threshold %0.2f',rthr,LR2_thr));
end; %if plot_flag==1;
end;%if (length(loc_valid)>0);
end;%  if plot_flag;

ik_out = ik(loc_valid(loc));
ij_out = ij(loc_valid(loc));
%disp(sprintf(' %% ncov max %d, using rthr %0.3f = %0.2f [%0.2f,%0.2f] --> [%d-x-%d]=%d',max(Lx2_ra(:)),LR2_thr,rthr,LR2_min,LR2_max,ij_out,ik_out,area_out));
