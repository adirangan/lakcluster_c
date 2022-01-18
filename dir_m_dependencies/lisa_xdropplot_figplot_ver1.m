%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run lisa_xdropplot_initialize and lisa_xdropplot_loadtrace first. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

dir_out_s0000_jpg = sprintf('%s/dir_jpg',dir_out_s0000);
if (~exist(dir_out_s0000_jpg,'dir')); disp(sprintf(' %% making %s',dir_out_s0000_jpg)); mkdir(dir_out_s0000_jpg); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% plot rdrop_a_ and cdrop_a_ ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
figure;clf;
subplot(2,2,1); plot(rdrop_a_,'.'); ll=size(mds_sort_,1); xlim([1,length(rdrop_a_)]);ylim([1,ll]); title('rdrop_a_ scatter','Interpreter','none');
subplot(2,2,3); ll=size(mds_sort_,1); imagesc(log2(1+adi_hist2d(1:length(rdrop_a_),rdrop_a_,[1,length(rdrop_a_)],256,[1,ll],256))); title('rdrop_a_ histogram','Interpreter','none');
subplot(2,2,2);plot(cdrop_a_,'.'); ll=n_snp; xlim([1,length(cdrop_a_)]);ylim([1,ll]); title('cdrop_a_ scatter','Interpreter','none');
subplot(2,2,4); ll=n_snp; imagesc(log2(1+adi_hist2d(1:length(cdrop_a_),cdrop_a_,[1,length(cdrop_a_)],256,[1,ll],256))); title('cdrop_a_ histogram','Interpreter','none');
set(gcf,'Position',1+[0,0,1024*1.5,768]);
print('-djpeg',sprintf('%s/%s_xdrop_A.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_xdrop_A.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);

if n_study>1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% plot study stack ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
rem_r_ = interp1(r_rem_,1:length(r_rem_),length(rdrop_a_):-1:1); rem_r_(find(~isfinite(rem_r_)))=iteration_max;
rdrop_study_ = study_index_(rdrop_a_);
length_study_rdrop_ = zeros(length(rdrop_a_),n_study);
fraction_study_rdrop_ = zeros(length(rdrop_a_),n_study);
for ns=1:n_study; length_study_rdrop_(1,ns) = length(find(rdrop_study_==ns)); end;%for ns=1:n_study;
fraction_study_rdrop_(1,:) = length_study_rdrop_(1,:)/length(rdrop_a_); 
for nr=1:length(rdrop_a_)-1;
length_study_rdrop_(nr+1,:) = length_study_rdrop_(nr+0,:);
length_study_rdrop_(nr+1,rdrop_study_(nr)) = length_study_rdrop_(nr+1,rdrop_study_(nr))-1;
fraction_study_rdrop_(nr+1,:) = length_study_rdrop_(nr+1,:)/(length(rdrop_a_)-nr); 
end;%for nr=1:length(rdrop_a_)-1;
cra_study_rdrop_ = [zeros(length(rdrop_a_),1) , cumsum(fraction_study_rdrop_,2) , ones(length(rdrop_a_),1)];
rem_r_edge_ = [0 , 0.5*(rem_r_(1:end-1) + rem_r_(2:end+0)) , rem_r_(end)+0.5 ] ;
cra_hsv = colormap('hsv'); clen_hsv = size(cra_hsv,1);
cra = zeros(n_study,3); clen = size(cra,1);
tab=1; dtab = floor(clen_hsv * (1-1/n_study)/2);
for ns=1:n_study;
cra(ns,:) = cra_hsv(tab,:); tab = tab + dtab; if tab>clen_hsv; tab = tab-clen_hsv; end;
end;%for ns=1:n_study;
%figure;
clf;
subplot(2,1,1);
colormap(gca,cra);
ygap = 0.0125;
for ns=1:n_study;
x_ra_{ns} = zeros(5,length(rdrop_a_));
x_ra_{ns}(1,:) = rem_r_edge_(1:end-1);
x_ra_{ns}(2,:) = rem_r_edge_(2:end+0);
x_ra_{ns}(3,:) = rem_r_edge_(2:end+0);
x_ra_{ns}(4,:) = rem_r_edge_(1:end-1);
x_ra_{ns}(5,:) = rem_r_edge_(1:end-1);
y_ra_{ns} = zeros(5,length(rdrop_a_));
y_ra_{ns}(1,:) = cra_study_rdrop_(:,ns+0);
y_ra_{ns}(2,:) = cra_study_rdrop_(:,ns+0);
y_ra_{ns}(3,:) = cra_study_rdrop_(:,ns+1);
y_ra_{ns}(4,:) = cra_study_rdrop_(:,ns+1);
y_ra_{ns}(5,:) = cra_study_rdrop_(:,ns+0);
c_ra_{ns} = max(1,min(clen,round(clen*(ns+0.0)/n_study)))*ones(1,length(rdrop_a_));
hold on; p=patch(x_ra_{ns},y_ra_{ns} + (ns-1)*ygap,c_ra_{ns});set(p,'LineStyle','none','MarkerEdgeColor','none'); hold off;
end;%for ns=1:n_study;
xlim([0,rem_r_(end)+0.5]); ylim([0,1+(n_study-1)*ygap]); set(gca,'Ytick',[]);
xlabel('number of iterations remaining'); title('study distribution');
subplot(2,1,2);
colormap(gca,cra);
ygap = 0.0125;
for ns=1:n_study;
x_ra_{ns} = zeros(5,length(rdrop_a_));
x_ra_{ns}(1,:) = (1:length(rdrop_a_)) - 0.5;
x_ra_{ns}(2,:) = (1:length(rdrop_a_)) + 0.5;
x_ra_{ns}(3,:) = (1:length(rdrop_a_)) + 0.5;
x_ra_{ns}(4,:) = (1:length(rdrop_a_)) - 0.5;
x_ra_{ns}(5,:) = (1:length(rdrop_a_)) - 0.5;
y_ra_{ns} = zeros(5,length(rdrop_a_));
y_ra_{ns}(1,:) = cra_study_rdrop_(:,ns+0);
y_ra_{ns}(2,:) = cra_study_rdrop_(:,ns+0);
y_ra_{ns}(3,:) = cra_study_rdrop_(:,ns+1);
y_ra_{ns}(4,:) = cra_study_rdrop_(:,ns+1);
y_ra_{ns}(5,:) = cra_study_rdrop_(:,ns+0);
c_ra_{ns} = max(1,min(clen,round(clen*(ns+0.0)/n_study)))*ones(1,length(rdrop_a_));
hold on; p=patch(x_ra_{ns},y_ra_{ns} + (ns-1)*ygap,c_ra_{ns});set(p,'LineStyle','none','MarkerEdgeColor','none'); hold off;
end;%for ns=1:n_study;
xlim([0.5,length(rdrop_a_)+0.5]); ylim([0,1+(n_study-1)*ygap]); set(gca,'Ytick',[]);
xlabel('number of patients eliminated'); title('study distribution');
set(gcf,'Position',1+[0,0,1024,768]);
print('-djpeg',sprintf('%s/%s_study_stack_A.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_study_stack_A.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);
end;%if n_study>1;

if n_study>1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% plot study distribution ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
clear cra_hsv clen_hsv cra clen tab dtab;
cra_hsv = colormap('hsv'); clen_hsv = size(cra_hsv,1);
cra = zeros(n_study,3); clen = size(cra,1);
tab=1; dtab = floor(clen_hsv * (1-1/n_study)/2);
for ns=1:n_study;
cra(ns,:) = cra_hsv(tab,:); tab = tab + dtab; if tab>clen_hsv; tab = tab-clen_hsv; end;
end;%for ns=1:n_study;
rcut_ = max(1,min(length(rdrop_a_),floor(length(rdrop_a_)./1.5.^((1:8)-1))));
figure;clf;
for np=1:8;
subplot(2,4,np);
clear iXX_ iXX
rcut = rcut_(np);
rij = rdrop_a_(max(1,length(rdrop_a_)-rcut):end); 
for ns=1:n_study; iXX(ns) = sum(rdrop_study_(max(1,length(rdrop_study_)-rcut):end)==ns); end;
n_prm = 1024*1.0;
iXX_ = zeros(n_study,n_prm);
for ni=1:n_prm;
prm_ = randperm(length(rdrop_study_)); rdrop_study_tmp_ = rdrop_study_(prm_(1:rcut));
for ns=1:n_study;iXX_(ns,ni) = sum(rdrop_study_tmp_==ns); end;
end;%for ni=1:n_prm;
iXX_avg = mean(iXX_,2);
iXX_50 = prctile(iXX_,50,2);
iXX_85 = prctile(iXX_,85,2);
iXX_15 = prctile(iXX_,15,2);
iXX_95 = prctile(iXX_,95,2);
iXX_05 = prctile(iXX_,05,2);
hold on;
for ns=1:n_study;
xvals = ns + [0;1;1;0] - 0.5; yvals = [0;0;1;1]*iXX(ns); p=patch(xvals,yvals,cra(ns,:)); set(p,'LineStyle','none','MarkerEdgeColor','none');
xvals = ns + [0.5;0.5] - 0.5; yvals = [iXX_15(ns),iXX_85(ns)]; l = line(xvals,yvals,'LineWidth',3,'LineStyle','-','Color',[0,0,0]);
xvals = ns + [0.25;0.75] - 0.5; yvals = [iXX_15(ns),iXX_15(ns)]; l = line(xvals,yvals,'LineWidth',3,'LineStyle','-','Color',[0,0,0]);
xvals = ns + [0.25;0.75] - 0.5; yvals = [iXX_85(ns),iXX_85(ns)]; l = line(xvals,yvals,'LineWidth',3,'LineStyle','-','Color',[0,0,0]);
xvals = ns + [0.5;0.5] - 0.5; yvals = [iXX_05(ns),iXX_95(ns)]; l = line(xvals,yvals,'LineWidth',1,'LineStyle','-','Color',[0,0,0]);
xvals = ns + [0.35;0.65] - 0.5; yvals = [iXX_05(ns),iXX_05(ns)]; l = line(xvals,yvals,'LineWidth',1,'LineStyle','-','Color',[0,0,0]);
xvals = ns + [0.35;0.65] - 0.5; yvals = [iXX_95(ns),iXX_95(ns)]; l = line(xvals,yvals,'LineWidth',1,'LineStyle','-','Color',[0,0,0]);
end;%for ns=1:n_study;
hold off;
set(gca,'Xtick',1:n_study); xlim([0,1+n_study]);
xlabel('study index'); ylabel('number of patients');
title(sprintf('patients remaining %d',length(rij)));
end;%for np=1:8;
set(gcf,'Position',1 + [0,0,1024*2,1024*1]);
print('-djpeg',sprintf('%s/%s_study_dist_A.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_study_dist_A.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);
end;%if n_study>1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% plot BD phenotype stack ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
rem_r_ = interp1(r_rem_,1:length(r_rem_),length(rdrop_a_):-1:1); rem_r_(find(~isfinite(rem_r_)))=iteration_max;
rdrop_pheno_ = pheno_index_(rdrop_a_);
length_pheno_rdrop_ = zeros(length(rdrop_a_),n_pheno);
fraction_pheno_rdrop_ = zeros(length(rdrop_a_),n_pheno);
for ns=1:n_pheno; length_pheno_rdrop_(1,ns) = length(find(rdrop_pheno_==ns)); end;%for ns=1:n_pheno;
fraction_pheno_rdrop_(1,:) = length_pheno_rdrop_(1,:)/length(rdrop_a_); 
for nr=1:length(rdrop_a_)-1;
length_pheno_rdrop_(nr+1,:) = length_pheno_rdrop_(nr+0,:);
length_pheno_rdrop_(nr+1,rdrop_pheno_(nr)) = length_pheno_rdrop_(nr+1,rdrop_pheno_(nr))-1;
fraction_pheno_rdrop_(nr+1,:) = length_pheno_rdrop_(nr+1,:)/(length(rdrop_a_)-nr); 
end;%for nr=1:length(rdrop_a_)-1;
cra_pheno_rdrop_ = [zeros(length(rdrop_a_),1) , cumsum(fraction_pheno_rdrop_,2) , ones(length(rdrop_a_),1)];
rem_r_edge_ = [0 , 0.5*(rem_r_(1:end-1) + rem_r_(2:end+0)) , rem_r_(end)+0.5 ] ;
cra_hsv = colormap('hsv'); clen_hsv = size(cra_hsv,1);
cra = zeros(n_pheno,3); clen = size(cra,1);
tab=1; dtab = floor(clen_hsv * (1-1/n_pheno)/2);
for ns=1:n_pheno;
cra(ns,:) = cra_hsv(tab,:); tab = tab + dtab; if tab>clen_hsv; tab = tab-clen_hsv; end;
end;%for ns=1:n_pheno;
%figure;
clf;
subplot(2,1,1);
colormap(cra);
ygap = 0.0125;
for ns=1:n_pheno;
x_ra_{ns} = zeros(5,length(rdrop_a_));
x_ra_{ns}(1,:) = rem_r_edge_(1:end-1);
x_ra_{ns}(2,:) = rem_r_edge_(2:end+0);
x_ra_{ns}(3,:) = rem_r_edge_(2:end+0);
x_ra_{ns}(4,:) = rem_r_edge_(1:end-1);
x_ra_{ns}(5,:) = rem_r_edge_(1:end-1);
y_ra_{ns} = zeros(5,length(rdrop_a_));
y_ra_{ns}(1,:) = cra_pheno_rdrop_(:,ns+0);
y_ra_{ns}(2,:) = cra_pheno_rdrop_(:,ns+0);
y_ra_{ns}(3,:) = cra_pheno_rdrop_(:,ns+1);
y_ra_{ns}(4,:) = cra_pheno_rdrop_(:,ns+1);
y_ra_{ns}(5,:) = cra_pheno_rdrop_(:,ns+0);
c_ra_{ns} = max(1,min(clen,round(clen*(ns+0.0)/n_pheno)))*ones(1,length(rdrop_a_));
hold on; p=patch(x_ra_{ns},y_ra_{ns} + (ns-1)*ygap,c_ra_{ns});set(p,'LineStyle','none','MarkerEdgeColor','none'); hold off;
end;%for ns=1:n_pheno;
xlim([0,rem_r_(end)+0.5]); ylim([0,1+(n_pheno-1)*ygap]); set(gca,'Ytick',[]);
xlabel('number of iterations remaining'); title('pheno distribution');
subplot(2,1,2);
colormap(cra);
ygap = 0.0125;
for ns=1:n_pheno;
x_ra_{ns} = zeros(5,length(rdrop_a_));
x_ra_{ns}(1,:) = (1:length(rdrop_a_)) - 0.5;
x_ra_{ns}(2,:) = (1:length(rdrop_a_)) + 0.5;
x_ra_{ns}(3,:) = (1:length(rdrop_a_)) + 0.5;
x_ra_{ns}(4,:) = (1:length(rdrop_a_)) - 0.5;
x_ra_{ns}(5,:) = (1:length(rdrop_a_)) - 0.5;
y_ra_{ns} = zeros(5,length(rdrop_a_));
y_ra_{ns}(1,:) = cra_pheno_rdrop_(:,ns+0);
y_ra_{ns}(2,:) = cra_pheno_rdrop_(:,ns+0);
y_ra_{ns}(3,:) = cra_pheno_rdrop_(:,ns+1);
y_ra_{ns}(4,:) = cra_pheno_rdrop_(:,ns+1);
y_ra_{ns}(5,:) = cra_pheno_rdrop_(:,ns+0);
c_ra_{ns} = max(1,min(clen,round(clen*(ns+0.0)/n_pheno)))*ones(1,length(rdrop_a_));
hold on; p=patch(x_ra_{ns},y_ra_{ns} + (ns-1)*ygap,c_ra_{ns});set(p,'LineStyle','none','MarkerEdgeColor','none'); hold off;
end;%for ns=1:n_pheno;
xlim([0.5,length(rdrop_a_)+0.5]); ylim([0,1+(n_pheno-1)*ygap]); set(gca,'Ytick',[]);
xlabel('number of patients remaining'); title('pheno distribution');
set(gcf,'Position',1+[0,0,1024,768]);
print('-djpeg',sprintf('%s/%s_pheno_stack_A.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_pheno_stack_A.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% plot mds distribution ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
xl = 0.1*[-0.5,1]; yl = 0.065*[-1,1];
v1=1;v2=2;
xm = median(mds_sort_(rdrop_a_,v1)); ym = median(mds_sort_(rdrop_a_,v2));
rcut_ = max(1,min(length(rdrop_a_),floor(length(rdrop_a_)./1.5.^((1:8)-1))));
figure;clf;
for np=1:8;
rcut = rcut_(np);
rij = rdrop_a_(max(1,length(rdrop_a_)-rcut):end); 
subplot(2,4,np); colormap(1-colormap('gray')); hold on;
ah = adi_hist2d(mds_sort_(rij,v1),mds_sort_(rij,v2),xl,96,yl,96);
imagesc(log2(1+log2(1+ah))); 
xlg = xlim(); ylg = ylim(); 
xlm = xlg(1) + diff(xlg)*(xm-xl(1))/diff(xl);
ylm = ylg(1) + diff(ylg)*(ym-yl(1))/diff(yl);
l=line(xlg,ylm*[1,1],'LineWidth',0.25,'Color',0.50*[1,1,1]);
l=line(xlm*[1,1],ylg,'LineWidth',0.25,'Color',0.50*[1,1,1]);
xlim(xlg);ylim(ylg);
set(gca,'Xtick',[],'Ytick',[]); axis square; hold off;
title(sprintf('patients %d',length(rij)));
end;%for np=1:8;
set(gcf,'Position',1+[0,0,1024,470]);
print('-djpeg',sprintf('%s/%s_mds_dist_A.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_mds_dist_A.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);
figure;clf;
for np=1:8;
rcut = rcut_(np);
rij = rdrop_a_(max(1,length(rdrop_a_)-rcut):end); 
subplot(2,4,np); cra = colormap('hot'); colormap(cra(end:-1:1,:)); hold on;
contourf(log2(1+adi_hist2d(mds_sort_(rij,v1),mds_sort_(rij,v2),xl,32,yl,32)),4); 
xlg = xlim(); ylg = ylim(); 
xlm = xlg(1) + diff(xlg)*(xm-xl(1))/diff(xl);
ylm = ylg(1) + diff(ylg)*(ym-yl(1))/diff(yl);
l=line(xlg,ylm*[1,1],'LineWidth',0.25,'Color',0.50*[1,1,1]);
l=line(xlm*[1,1],ylg,'LineWidth',0.25,'Color',0.50*[1,1,1]);
xlim(xlg);ylim(ylg);
set(gca,'Xtick',[],'Ytick',[]);  axis square; hold off;
title(sprintf('patients %d',length(rij)));
end;%for np=1:8;
set(gcf,'Position',1+[0,0,1024,470]);
print('-djpeg',sprintf('%s/%s_mds_dist_B.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_mds_dist_B.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% plot mds distribution for a particular iteration ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
xl = 0.1*[-0.5,1]; yl = 0.065*[-1,1];
v1=1;v2=2;
xm = median(mds_sort_(rdrop_a_,v1)); ym = median(mds_sort_(rdrop_a_,v2));
% put the particular iteration into rcut_ ;
tmp_iteration = 275;
rcut_ = [length(rdrop_a_) , trace_{1}(1+tmp_iteration,2)];
figure;clf;
for np=1:2;
rcut = rcut_(np);
rij = rdrop_a_(max(1,length(rdrop_a_)-rcut):end); 
subplot(1,2,np); colormap(1-colormap('gray')); hold on;
ah = adi_hist2d(mds_sort_(rij,v1),mds_sort_(rij,v2),xl,96,yl,96);
imagesc(log2(1+log2(1+ah))); 
xlg = xlim(); ylg = ylim(); 
xlm = xlg(1) + diff(xlg)*(xm-xl(1))/diff(xl);
ylm = ylg(1) + diff(ylg)*(ym-yl(1))/diff(yl);
l=line(xlg,ylm*[1,1],'LineWidth',0.25,'Color',0.50*[1,1,1]);
l=line(xlm*[1,1],ylg,'LineWidth',0.25,'Color',0.50*[1,1,1]);
xlim(xlg);ylim(ylg);
set(gca,'Xtick',[],'Ytick',[]); axis square; hold off;
title(sprintf('patients %d',length(rij)));
end;%for np=1:2;
set(gcf,'Position',1+[0,0,1024,470]);
print('-djpeg',sprintf('%s/%s_mds_dist_A_ni%.3d.jpg',dir_out_s0000_jpg,string_name_s0000,tmp_iteration));
print('-depsc',sprintf('%s/%s_mds_dist_A_ni%.3d.eps',dir_out_s0000_jpg,string_name_s0000,tmp_iteration));
drawnow();pause(1);
figure;clf;
for np=1:2;
rcut = rcut_(np);
rij = rdrop_a_(max(1,length(rdrop_a_)-rcut):end); 
subplot(1,2,np); cra = colormap('hot'); colormap(cra(end:-1:1,:)); hold on;
contourf(log2(1+adi_hist2d(mds_sort_(rij,v1),mds_sort_(rij,v2),xl,32,yl,32)),4); 
xlg = xlim(); ylg = ylim(); 
xlm = xlg(1) + diff(xlg)*(xm-xl(1))/diff(xl);
ylm = ylg(1) + diff(ylg)*(ym-yl(1))/diff(yl);
l=line(xlg,ylm*[1,1],'LineWidth',0.25,'Color',0.50*[1,1,1]);
l=line(xlm*[1,1],ylg,'LineWidth',0.25,'Color',0.50*[1,1,1]);
xlim(xlg);ylim(ylg);
set(gca,'Xtick',[],'Ytick',[]);  axis square; hold off;
title(sprintf('patients %d',length(rij)));
end;%for np=1:2;
set(gcf,'Position',1+[0,0,1024,470]);
print('-djpeg',sprintf('%s/%s_mds_dist_B_ni%.3d.jpg',dir_out_s0000_jpg,string_name_s0000,tmp_iteration));
print('-depsc',sprintf('%s/%s_mds_dist_B_ni%.3d.eps',dir_out_s0000_jpg,string_name_s0000,tmp_iteration));
drawnow();pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% plot trace distribution ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if (n_finished>2);
QR_ = zeros(n_finished,iteration_max);
QC_ = zeros(n_finished,iteration_max);
for nf=1:n_finished;
ns = trace_finished_ij_(nf); trace_tmp_ = trace_{ns};
QR_(nf,:) = transpose(trace_tmp_(:,4));
QC_(nf,:) = transpose(trace_tmp_(:,5));
end;%for nf=1:n_finished;
%%%%%%%%%%%%%%%% ;
PR_ = zeros(n_finished,iteration_max);
PC_ = zeros(n_finished,iteration_max);
for ni=1:iteration_max;
[~,tmp_ij] = sort(QR_(:,ni),'ascend'); [~,PR_(:,ni)] = sort(tmp_ij,'ascend');
[~,tmp_ij] = sort(QC_(:,ni),'ascend'); [~,PC_(:,ni)] = sort(tmp_ij,'ascend');
end;%for ni=1:iteration_max;
%%%%%%%%%%%%%%%% ;
QR_ = transpose(QR_); QC_ = transpose(QC_);
PR_ = transpose(PR_); PC_ = transpose(PC_);
QR_avg_ = mean(QR_,2); QR_std_ = std(QR_,[],2); QC_avg_ = mean(QC_,2); QC_std_ = std(QC_,[],2);
ZR_ = (QR_ - repmat(QR_avg_,1,n_finished))./repmat(QR_std_,1,n_finished); ZC_ = (QC_ - repmat(QC_avg_,1,n_finished))./repmat(QC_std_,1,n_finished);
%%%%%%%%%%%%%%%% ;
figure;clf;
subplot(2,1,1); plot(max(t_ret_)-t_ret_(end:-1:1),ZR_(t_ret_,1),'ro-',max(t_ret_)-t_ret_(end:-1:1),ZR_(t_ret_,2:end),'b-'); 
t_tick_ = 0:25:max(t_ret_); set(gca,'XTick',t_tick_,'XTickLabel',t_tick_); hold on; for nt=1:length(t_tick_); plot(t_tick_(nt)*[1,1],[min(ZR_(:)),max(ZR_(:))],'k-'); end; ylim([min(ZR_(:)),max(ZR_(:))]);
xlabel('row iteration number'); ylabel('ZR_','Interpreter','none'); xlim([0,max((t_ret_))+1]);
subplot(2,1,2); plot(max(t_ret_)-t_ret_(end:-1:1),ZC_(t_ret_,1),'ro-',max(t_ret_)-t_ret_(end:-1:1),ZC_(t_ret_,2:end),'b-'); 
t_tick_ = 0:25:max(t_ret_); set(gca,'XTick',t_tick_,'XTickLabel',t_tick_); hold on; for nt=1:length(t_tick_); plot(t_tick_(nt)*[1,1],[min(ZC_(:)),max(ZC_(:))],'k-'); end; ylim([min(ZC_(:)),max(ZC_(:))]);
xlabel('col iteration number'); ylabel('ZC_','Interpreter','none'); xlim([0,max((t_ret_))+1]);
set(gcf,'Position',1+[0,0,1024*2,768]);
print('-djpeg',sprintf('%s/%s_trace_C.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_trace_C.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);
%%%%%%%%%%%%%%%%;
figure;clf;
disp(sprintf(' %% out_trace results; %d pats -x- %d snps ; %d iterations',max(rdrop_b_),max(cdrop_b_),max(t_ret_)));
subplot(2,4,1); plot(max(r_rem_)-r_rem_(t_ret_),ZR_(t_ret_,1),'ro-',max(r_rem_)-r_rem_(t_ret_),ZR_(t_ret_,2:end),'b-'); xlabel('rows eliminated'); ylabel('ZR_','Interpreter','none'); xlim([0,max(max(r_rem_)-r_rem_(t_ret_))+1]);
subplot(2,4,2); plot(max(r_rem_)-r_rem_(t_ret_),PR_(t_ret_,1),'ro-',max(r_rem_)-r_rem_(t_ret_),PR_(t_ret_,2:end),'b-');  xlabel('rows eliminated'); ylabel('PR_','Interpreter','none'); xlim([0,max(max(r_rem_)-r_rem_(t_ret_))+1]);ylim([0,n_finished+1]);
subplot(2,4,3); plot(max(t_ret_)-t_ret_(end:-1:1),ZR_(t_ret_,1),'ro-',max(t_ret_)-t_ret_(end:-1:1),ZR_(t_ret_,2:end),'b-'); xlabel('row iteration number'); ylabel('ZR_','Interpreter','none'); xlim([0,max((t_ret_))+1]);
subplot(2,4,4); plot(max(t_ret_)-t_ret_(end:-1:1),PR_(t_ret_,1),'ro-',max(t_ret_)-t_ret_(end:-1:1),PR_(t_ret_,2:end),'b-'); xlabel('row iteration number'); ylabel('PR_','Interpreter','none'); xlim([0,max((t_ret_))+1]);ylim([0,n_finished+1]);
subplot(2,4,5); plot(max(c_rem_)-c_rem_(t_ret_),ZC_(t_ret_,1),'ro-',max(c_rem_)-c_rem_(t_ret_),ZC_(t_ret_,2:end),'b-'); xlabel('cols eliminated'); ylabel('ZC_','Interpreter','none'); xlim([0,max(max(c_rem_)-c_rem_(t_ret_))+1]);
subplot(2,4,6); plot(max(c_rem_)-c_rem_(t_ret_),PC_(t_ret_,1),'ro-',max(c_rem_)-c_rem_(t_ret_),PC_(t_ret_,2:end),'b-'); xlabel('cols eliminated'); ylabel('PC_','Interpreter','none'); xlim([0,max(max(c_rem_)-c_rem_(t_ret_))+1]);ylim([0,n_finished+1]);
subplot(2,4,7); plot(max(t_ret_)-t_ret_(end:-1:1),ZC_(t_ret_,1),'ro-',max(t_ret_)-t_ret_(end:-1:1),ZC_(t_ret_,2:end),'b-'); xlabel('col iteration number'); ylabel('ZC_','Interpreter','none'); xlim([0,max((t_ret_))+1]);
subplot(2,4,8); plot(max(t_ret_)-t_ret_(end:-1:1),PC_(t_ret_,1),'ro-',max(t_ret_)-t_ret_(end:-1:1),PC_(t_ret_,2:end),'b-'); xlabel('col iteration number'); ylabel('PC_','Interpreter','none'); xlim([0,max((t_ret_))+1]);ylim([0,n_finished+1]);
set(gcf,'Position',1+[0,0,1024*2,768]);
print('-djpeg',sprintf('%s/%s_trace_A.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_trace_A.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);
%%%%%%%%%%%%%%%%;
figure;clf;
iteration_alo = max(1,floor(iteration_max*0.05));
iteration_ahi = min(iteration_max-1,floor(iteration_max*0.95));
iteration_avg_ = 1:iteration_ahi; 
dr_ = diff(r_rem_); dc_ = diff(c_rem_); 
dr2_ = transpose(dr_(iteration_avg_))/sum(dr_(iteration_avg_)); 
dc2_ = transpose(dc_(iteration_avg_))/sum(dc_(iteration_avg_));
iteration_top_ = iteration_alo:iteration_ahi;
subplot(2,2,1); title('top-vs-avg ZR_','Interpreter','none');
hold on;
plot(max(ZR_(iteration_top_,2:end),[],1),dr2_*ZR_(iteration_avg_,2:end),'bx');
plot(max(ZR_(iteration_top_,1    ),[],1),dr2_*ZR_(iteration_avg_,1    ),'rx');
hold off; xlabel('top'); ylabel('avg');
subplot(2,2,2); title('top-vs-avg ZC_','Interpreter','none');
hold on;
plot(max(ZC_(iteration_top_,2:end),[],1),dc2_*ZC_(iteration_avg_,2:end),'bx');
plot(max(ZC_(iteration_top_,1    ),[],1),dc2_*ZC_(iteration_avg_,1    ),'rx');
hold off; xlabel('top'); ylabel('avg');
subplot(2,2,3); title('top-vs-avg PR_','Interpreter','none');
hold on;
plot(max(PR_(iteration_top_,2:end),[],1),dr2_*PR_(iteration_avg_,2:end),'bx');
plot(max(PR_(iteration_top_,1    ),[],1),dr2_*PR_(iteration_avg_,1    ),'rx');
hold off; xlabel('top'); ylabel('avg');
subplot(2,2,4); title('top-vs-avg PC_','Interpreter','none');
hold on;
plot(max(PC_(iteration_top_,2:end),[],1),dc2_*PC_(iteration_avg_,2:end),'bx');
plot(max(PC_(iteration_top_,1    ),[],1),dc2_*PC_(iteration_avg_,1    ),'rx');
hold off; xlabel('top'); ylabel('avg');
set(gcf,'Position',1+[0,0,1024*1.5,768]);
print('-djpeg',sprintf('%s/%s_trace_B.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_trace_B.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);
end;%if (n_finished>2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% plot colored scatterplot using ZR_ ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if (n_finished>2);
iteration_alo = max(1,floor(iteration_max*0.05));
iteration_ahi = min(iteration_max-1,floor(iteration_max*0.95));
iteration_avg_ = 1:iteration_ahi; 
dr_ = diff(r_rem_); dc_ = diff(c_rem_); 
dr2_ = transpose(dr_(iteration_avg_))/sum(dr_(iteration_avg_)); 
dc2_ = transpose(dc_(iteration_avg_))/sum(dc_(iteration_avg_));
iteration_top_ = iteration_alo:iteration_ahi;
tr_avg_s0000 = dr2_*ZR_(iteration_avg_,1);
tr_avg_ = dr2_*ZR_(iteration_avg_,2:end);
tr_avg_p_ = zeros(n_finished-1,1);
tr_top_s0000 = max(ZR_(iteration_top_,1),[],1); disp(tr_top_s0000); 
tr_top_ = max(ZR_(iteration_top_,2:end),[],1);
tr_top_p_ = zeros(n_finished-1,1);
p_ra_ = zeros(n_finished-1,1);
for np=(1:n_finished-1);
tr_top_p_(np) = (length(find(tr_top_>tr_top_(np))) + 0.5*length(find(tr_top_==tr_top_(np))))/(n_finished-1);
tr_avg_p_(np) = (length(find(tr_avg_>tr_avg_(np))) + 0.5*length(find(tr_avg_==tr_avg_(np))))/(n_finished-1);
tau_tmp = min([tr_top_p_(np),tr_avg_p_(np)]); tau_tmp = max(0.5/(n_finished-1),tau_tmp);
ls_rm = find(tr_top_(:)>=prctile(tr_top_(:),100*(1-tau_tmp)));
ls_ra = find(tr_avg_(:)>=prctile(tr_avg_(:),100*(1-tau_tmp)));
p_ra_(np) = length(unionall({ls_rm,ls_ra}))/(n_finished-1);
end;%for np=(1:n_finished-1);
tr_top_s0000_p = (length(find(tr_top_>tr_top_s0000)) + 0.5*length(find(tr_top_==tr_top_s0000)))/(n_finished-1);
tr_avg_s0000_p = (length(find(tr_avg_>tr_avg_s0000)) + 0.5*length(find(tr_avg_==tr_avg_s0000)))/(n_finished-1);
tau_tmp = min([tr_top_s0000_p,tr_avg_s0000_p]); tau_tmp = max(0.5/(n_finished-1),tau_tmp);
ls_rm = find(tr_top_(:)>=prctile(tr_top_(:),100*(1-tau_tmp)));
ls_ra = find(tr_avg_(:)>=prctile(tr_avg_(:),100*(1-tau_tmp)));
p_s0000 = length(unionall({ls_rm,ls_ra}))/(n_finished-1);
disp(sprintf('p_s0000 %0.3f',p_s0000));
figure;
hold on;
l10plim=[0,2];
cmap = colormap('jet'); clen = size(cmap,1); %cmap = cmap(end:-1:1,:);
Msize1 = 25; Msize2 = 35;
for np=1:(n_finished-1);
nc = max(1,min(clen,floor(clen*((-log10(p_ra_(np))-min(l10plim))/diff(l10plim)))));
plot(tr_top_(np),tr_avg_(np),'.','Color',cmap(nc,:),'MarkerSize',Msize1);
end;%for np=1:(n_finished-1);
for no_flag=0:0;%for no_flag=0:NO-1;
plot(tr_top_s0000,tr_avg_s0000,'x','Color',[0,0,0],'MarkerSize',Msize2);
plot(tr_top_s0000,tr_avg_s0000,'o','Color',[0,0,0],'MarkerSize',Msize2);
end;%for no_flag=0:NO-1;
hold off;
%xlim([-1.50,+3.25]);ylim([-1.75,+2.75]);
disp(sprintf('xlim %0.2f %0.2f ylim %0.2f %0.2f',xlim(),ylim()));
xlabel('row Z-score (top)'); ylabel('row Z-score (avg)'); 
title('scatterplot of distribution of top vs avg: bicluster (x) and permutations (dots)');
%set(gca,'Xtick',[],'Ytick',[]);
axis equal;
set(gcf,'Position',1+[0,0,1024,1024]);
print('-djpeg',sprintf('%s/%s_pval_scatter_A.jpg',dir_out_s0000_jpg,string_name_s0000));
print('-depsc',sprintf('%s/%s_pval_scatter_A.eps',dir_out_s0000_jpg,string_name_s0000));
drawnow();pause(1);
end;%if (n_finished>2);

%disp('returning');return;