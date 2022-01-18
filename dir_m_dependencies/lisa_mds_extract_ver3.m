% 20180304: recovering mds file ;
% also generating plot indicating which clusters lie where in covariate-space. ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% loading fam file for each study. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

dir_base = '/data/rangan/dir_bcc/dir_PGC_20180304';

fam_fname_ = {...
'v1/bip_bonn_eur_sr-qc.fam',...
'v1/bip_dub1_eur_sr-qc.fam',...
'v1/bip_edi1_eur_sr-qc.fam',...
'v1/bip_fat2_eur_sr-qc.fam',...
'v1/bip_fran_eur_sr-qc.fam',...
'v1/bip_gsk1_eur_sr-qc.fam',...
'v1/bip_icuk_eur_sr-qc.fam',...
'v1/bip_may1_eur_sr-qc.fam',...
'v1/bip_mich_eur_sr-qc.fam',...
'v1/bip_st2c_eur_sr-qc.fam',...
'v1/bip_stp1_eur_sr-qc.fam',...
'v1/bip_swa2_eur_sr-qc.fam',...
'v1/bip_swei_eur_sr-qc.fam',...
'v1/bip_top7_eur_sr-qc.fam',...
'v1/bip_top8_eur_sr-qc2.fam',...
'v1/bip_ucl2_eur_sr-qc2.fam',...
'v1/bip_uclo_eur_sr-qc.fam',...
'v1/bip_ume4_eur_sr-qc.fam',...
'v1/bip_usc2_eur_sr-qc.fam',...
'v1/bip_wtcc_eur_sr-qc.fam',...
'v1_boma/bip_bmau_eur_sr-qc.fam',...
'v1_boma/bip_bmg2_eur_sr-qc.fam',...
'v1_boma/bip_bmg3_eur_sr-qc.fam',...
'v1_boma/bip_bmpo_eur_sr-qc.fam',...
'v1_boma/bip_bmsp_eur_sr-qc2.fam',...
'v1_boma/bip_hal2_eur_sr-qc.fam',...
'v1_boma/bip_rom3_eur_sr-qc.fam',...
'v1_gain/bip_gain_eur_sr-qc.fam',...
'v1/ms.bip_butr_eur-qc-pcx.fam',...
'v1/ms.bip_uktr_eur-qc-pcx.fam',...
};

fam_label_ = {...
'bonn',...
'dub1',...
'edi1',...
'fat2',...
'fran',...
'gsk1',...
'icuk',...
'may1',...
'mich',...
'st2c',...
'stp1',...
'swa2',...
'swei',...
'top7',...
'top8',...
'ucl2',...
'uclo',...
'ume4',...
'usc2',...
'wtcc',...
'bmau',...
'bmg2',...
'bmg3',...
'bmpo',...
'bmsp',...
'hal2',...
'rom3',...
'gain',...
'butr',...
'uktr',...
};

n_study = length(fam_fname_);

fam_name_ = cell(n_study,1);
n_patient_ = zeros(n_study,1);
for nstudy=1:n_study;
fam_fname = fam_fname_{nstudy};
disp(sprintf(' %% reading %s',fam_fname));
fp = fopen(sprintf('%s/%s',dir_base,fam_fname),'r');
fam_{nstudy} = textscan(fp,'%s %s %s %s %d %d','headerlines',0);
n_patient_(nstudy) = length(fam_{nstudy}{1});
for np=1:n_patient_(nstudy);
fam_name_{nstudy}{np} = sprintf('%s%s%s',fam_{nstudy}{1}{np},'&',fam_{nstudy}{2}{np});
end;%for np=1:n_patient_(nstudy);
fclose(fp);
end;%for nstudy=1:n_study;

figure;
bar(1:n_study,n_patient_);
set(gca,'XTick',1:n_study,'XTickLabel',fam_label_,'XTickLabelRotation',90);
set(gcf,'Position',1+[0,0,1024,1024]);
ylabel('number');
title('number of patients in each study');
print('-depsc','./lisa_mds_extract_A.eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% loading newer mds_w2.cov, ; 
% naming convention should match fam files. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

fid = fopen('/data/rangan/dir_bcc/dir_PGC_20180304/mds_w2.cov');
A_ = textscan(fid,'%s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d','headerlines',1);
fclose(fid);

n_study_A = 32;

study_A_ = zeros(size(A_{24},1),1);
for np=1:size(A_{24},1);
study_A_(np) = 0;
for ns=1:n_study_A-1;
study_A_(np) = study_A_(np) + ns*cast(A_{23+ns}(np),'double');
end;%for ns=1:n_study_A-1;
end;%for np=1:size(A_{24},1);

n_1_ = zeros(n_study_A,1);
tmp_sum = zeros(size(A_{24},1),1);
for ns=1:n_study_A-1;
n_0 = sum(A_{23+ns}==0); n_1 = sum(A_{23+ns}==1);
disp(sprintf(' %% study index %.2d: number of (0) vs (1): (%.5d,%.5d)',ns,n_0,n_1));
tmp_sum = tmp_sum + cast(A_{23+ns},'double');
n_1_(1+ns) = n_1;
end;%for ns=1:n_study_A-1;
ns=0;
n_0 = sum(tmp_sum==1); n_1 = sum(tmp_sum==0);
disp(sprintf(' %% study index %.2d: number of (0) vs (1): (%.5d,%.5d)',ns,n_0,n_1));
n_1_(1+ns) = n_1;
figure;
bar(0:n_study_A-1,n_1_); xlim([-1,n_study_A]); xlabel('study index'); ylabel('number'); title('number per study');
set(gcf,'Position',1+[1,1,1024,512]);
print('-depsc','./lisa_mds_extract_study_A_histogram.eps');

n_patient_A = length(A_{1});
A_name_ = cell(n_patient_A,1); for np=1:n_patient_A; A_name_{np} = sprintf('%s%s%s',A_{1}{np},'&',A_{2}{np}); end;%for np=1:n_patient_A;

study_A_overlap_ = zeros(n_study,n_study_A);
for nstudy=1:n_study;
for nstudy_A=1:n_study_A;
study_A_overlap_(nstudy,nstudy_A) = length(intersect(fam_name_{nstudy},A_name_(find(study_A_==nstudy_A-1))));
end;%for nstudy_A=1:n_study_A;
end;%for nstudy=1:n_study;
imagesc(study_A_overlap_);
set(gca,'XTick',0:n_study_A-1,'XTickLabel',0:nstudy_A-1);
set(gca,'YTick',1:n_study,'YTickLabel',fam_label_);
set(gcf,'Position',1+[1,1,1024,524]);
print('-depsc','./lisa_mds_extract_study_A_overlap.eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% loading famex file. ;
% former name will be excluded based on relatedness to latter. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

fp = fopen('/data/rangan/dir_bcc/dir_PGC_20180304/prune.bfile.cobg.PGC_BIP32b.mepr.famex','r');
famex_ = textscan(fp,'%s %s %s %s','headerlines',0);
fclose(fp);
n_famex = length(famex_{1});
famex_name_former_ = cell(n_famex,1);
famex_name_latter_ = cell(n_famex,1);
for np=1:n_famex;
famex_name_former_{np} = sprintf('%s%s%s',famex_{1}{np},'&',famex_{2}{np});
famex_name_latter_{np} = sprintf('%s%s%s',famex_{3}{np},'&',famex_{4}{np});
end;%for np=1:n_famex;

famex_A_overlap_ = zeros(n_study,1);
for nstudy=1:n_study;
famex_A_overlap_(nstudy) = length(intersect(fam_name_{nstudy},famex_name_former_));
end;%for nstudy=1:n_study;
bar(1:n_study,famex_A_overlap_);
set(gca,'XTick',1:n_study,'XTickLabel',fam_label_,'XTickLabelRotation',90);
set(gcf,'Position',1+[1,1,1024,524]);
print('-depsc','./lisa_mds_extract_famex_A_overlap.eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% loading older mds_for_pract.cov, ; 
% naming convention should not match fam files, ;
% should not be used. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

fid = fopen('/data/rangan/dir_bcc/dir_PGC_20180304/mds_for_pract.cov');
B_ = textscan(fid,'%s %s %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f','headerlines',1);
fclose(fid);

n_study_B = 11;

study_B_ = zeros(size(B_{3},1),1);
for np=1:size(B_{3},1);
study_B_(np) = 0;
for ns=1:n_study_B-1;
study_B_(np) = study_B_(np) + ns*cast(B_{2+ns}(np),'double');
end;%for ns=1:n_study_B-1;
end;%for np=1:size(B_{3},1);

n_1_ = zeros(n_study_B,1);
tmp_sum = zeros(size(B_{3},1),1);
for ns=1:n_study_B-1;
n_0 = sum(B_{2+ns}==0); n_1 = sum(B_{2+ns}==1);
disp(sprintf(' %% study index %.2d: number of (0) vs (1): (%.5d,%.5d)',ns,n_0,n_1));
tmp_sum = tmp_sum + cast(B_{2+ns},'double');
n_1_(1+ns) = n_1;
end;%for ns=1:n_study_B-1;
ns=0;
n_0 = sum(tmp_sum==1); n_1 = sum(tmp_sum==0);
disp(sprintf(' %% study index %.2d: number of (0) vs (1): (%.5d,%.5d)',ns,n_0,n_1));
n_1_(1+ns) = n_1;
figure;
bar(0:n_study_B-1,n_1_); xlim([-1,n_study_B]); xlabel('study index'); ylabel('number'); title('number per study');
set(gcf,'Position',1+[1,1,1024,512]);
print('-depsc','./lisa_mds_extract_study_B_histogram.eps');

n_patient_B = length(B_{1}); 
B_name_ = cell(n_patient_B,1); for np=1:n_patient_B; B_name_{np} = sprintf('%s%s%s',B_{1}{np},'&',B_{2}{np}); end;%for np=1:n_patient_B;

study_B_overlap_ = zeros(n_study,n_study_B);
for nstudy=1:n_study;
for nstudy_B=1:n_study_B;
study_B_overlap_(nstudy,nstudy_B) = length(intersect(fam_name_{nstudy},B_name_(find(study_B_==nstudy_B-1))));
end;%for nstudy_B=1:n_study_B;
end;%for nstudy=1:n_study;
imagesc(study_B_overlap_);
set(gca,'XTick',0:n_study_B-1,'XTickLabel',0:nstudy_B-1);
set(gca,'YTick',1:n_study,'YTickLabel',fam_label_);
set(gcf,'Position',1+[1,1,1024,524]);
print('-depsc','./lisa_mds_extract_study_B_overlap.eps');

assert(length(intersect(A_name_,B_name_))==0);
n_patient_mds = n_patient_A + n_patient_B;
mds_ = zeros(n_patient_mds,20);

rij_cas = zeros(n_patient_mds,1); rij_con = zeros(n_patient_mds,1); rij_fam = zeros(n_patient_mds,1);
for np=1:n_patient_A; 
rij_cas(np) = strcmp(A_name_{np}(1:3),'cas'); 
rij_con(np) = strcmp(A_name_{np}(1:3),'con'); 
rij_fam(np) = strcmp(A_name_{np}(1:3),'fam'); 
end;%for np=1:n_patient_A;
for np=1:n_patient_B; 
rij_cas(n_patient_A+np) = strcmp(B_name_{np}(1:3),'cas'); 
rij_con(n_patient_A+np) = strcmp(B_name_{np}(1:3),'con'); 
rij_fam(n_patient_A+np) = strcmp(B_name_{np}(1:3),'fam'); 
end;%for np=1:n_patient_B;

rng(0);
for np=1:n_patient_A;
for ip=1:20; mds_(np,ip) = A_{3+ip}(np); end;% for ip=1:20;
end;%for np=1:n_patient_A;
mds_avg_ = mean(mds_(1:n_patient_A,:));
mds_std_ = std(mds_(1:n_patient_A,:));
for np=1:n_patient_B;
for ip=1:5; mds_(n_patient_A+np,ip) = B_{12+ip}(np); end;% for ip=1:5;
mds_(n_patient_A+np,6:20) = mds_avg_(6:20) + mds_std_(6:20).*randn(1,15);
end;%for np=1:n_patient_B;

%{
cra = colormap('lines');
figure;
v1=1;v2=2;v3=3; %v1=4;v2=5;v3=6; 
hold on;
for nstudy=1:n_study;
[~,tmp_max] = max(study_A_overlap_(nstudy,:));
tmp_ij = find(study_A_== tmp_max-1);
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'.','MarkerSize',15,'Color',cra(nstudy,:))
end;%for nstudy=1:n_study;
hold off;
axis vis3d;
 %}

%v1=1;v2=2;plot(mds_(:,v1),mds_(:,v2),'.','MarkerSize',5);
cra = colormap('lines');
figure;
v1=1;v2=2;
prows = 5; pcols = 6;
for nstudy=1:n_study;
[~,tmp_max] = max(study_A_overlap_(nstudy,:));
tmp_ij = find(study_A_== tmp_max-1);
subplot(prows,pcols,nstudy);
hold on;
plot(mds_(tmp_ij,v1),mds_(tmp_ij,v2),'.','MarkerSize',5,'Color',cra(nstudy,:));
grid;
hold off;
xlim([-0.04,0.08]);ylim([-0.04,0.08]);
title(fam_label_{nstudy});
end;%for nstudy=1:n_study;
set(gcf,'Position',1+[0,0,512*3,512*2]);
print('-depsc','./lisa_mds_extract_mds12.eps');

%v1=3;v2=4;plot(mds_(:,v1),mds_(:,v2),'.','MarkerSize',5);
cra = colormap('lines');
figure;
v1=3;v3=4;
prows = 5; pcols = 6;
for nstudy=1:n_study;
[~,tmp_max] = max(study_A_overlap_(nstudy,:));
tmp_ij = find(study_A_== tmp_max-1);
subplot(prows,pcols,nstudy);
hold on;
plot(mds_(tmp_ij,v1),mds_(tmp_ij,v2),'.','MarkerSize',5,'Color',cra(nstudy,:));
grid;
hold off;
xlim([-0.06,0.08]);ylim([-0.05,0.06]);
title(fam_label_{nstudy});
end;%for nstudy=1:n_study;
set(gcf,'Position',1+[0,0,512*3,512*2]);
print('-depsc','./lisa_mds_extract_mds34.eps');


mds_name_ = [A_name_;B_name_];
save('/data/rangan/dir_bcc/dir_PGC_20180304/mds.mat','n_patient_','n_patient_A','n_patient_B','n_patient_mds','mds_','rij_cas','rij_con','rij_fam','mds_name_','n_famex','famex_name_former_','famex_name_latter_');

%{
figure;
v1=1;v2=2;v3=3; %v1=4;v2=5;v3=6; 
subplot(1,2,1);
hold on;
tmp_ij = intersect(find(rij_cas),1:n_patient_A);
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'r.')
tmp_ij = intersect(find(rij_con),1:n_patient_A);
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'b.'); 
tmp_ij = intersect(find(rij_fam),1:n_patient_A);
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'k.'); 
axis vis3d;
hold off;
subplot(1,2,2);
hold on;
tmp_ij = intersect(find(rij_cas),n_patient_A + [1:n_patient_B]);
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'r.')
tmp_ij = intersect(find(rij_con),n_patient_A + [1:n_patient_B]);
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'b.'); 
tmp_ij = intersect(find(rij_fam),n_patient_A + [1:n_patient_B]);
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'k.'); 
axis vis3d;
hold off;
 %}

%{
figure;
v1=1;v2=2;v3=3; %v1=4;v2=5;v3=6; 
hold on;
tmp_ij = 1:n_patient_A;
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'r.')
tmp_ij = n_patient_A + [1:n_patient_B];
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'b.'); 
hold off;
axis vis3d;
 %}

%{
figure;
v1=1;v2=2;v3=3; %v1=4;v2=5;v3=6; 
hold on;
tmp_ij = 1:n_patient_A;
plot3(mds_(tmp_ij,v1),mds_(tmp_ij,v2),mds_(tmp_ij,v3),'k.')
hold off;
axis vis3d;
 %}
