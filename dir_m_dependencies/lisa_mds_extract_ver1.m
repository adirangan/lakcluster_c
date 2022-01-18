% 20180304: recovering mds file ;
% run from lisa: ;
% diff v1/chr23/prune.bfile.cobg.PGC_BIP32b.menv.mds_cov v1/prune.bfile.cobg.PGC_BIP32b.menv.mds_cov
% scp -p v1/chr23/prune.bfile.cobg.PGC_BIP32b.menv.mds_cov rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_PGC_20180304/mds_w2.cov

fid = fopen('/data/rangan/dir_bcc/dir_PGC_20180304/mds_w2.cov');
A_ = textscan(fid,'%s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d','headerlines',1);
fclose(fid);

n_1_ = zeros(32,1);
tmp_sum = zeros(size(A_{24},1),1);
for ns=1:31;
n_0 = sum(A_{23+ns}==0); n_1 = sum(A_{23+ns}==1);
disp(sprintf(' %% study index %.2d: number of (0) vs (1): (%.5d,%.5d)',ns,n_0,n_1));
tmp_sum = tmp_sum + cast(A_{23+ns},'double');
n_1_(1+ns) = n_1;
end;%for ns=1:31;
ns=0;
n_0 = sum(tmp_sum==1); n_1 = sum(tmp_sum==0);
disp(sprintf(' %% study index %.2d: number of (0) vs (1): (%.5d,%.5d)',ns,n_0,n_1));
n_1_(1+ns) = n_1;
figure;
bar(0:31,n_1_); xlim([-1,32]); xlabel('study index'); ylabel('number'); title('number per study');
set(gcf,'Position',1+[1,1,1024,512]);
print('-depsc','./mds_A_study_index_histogram.eps');

fid = fopen('/data/rangan/dir_bcc/dir_PGC_20180304/mds_for_pract.cov');
B_ = textscan(fid,'%s %s %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f','headerlines',1);
fclose(fid);

n_1_ = zeros(11,1);
tmp_sum = zeros(size(B_{3},1),1);
for ns=1:10;
n_0 = sum(B_{2+ns}==0); n_1 = sum(B_{2+ns}==1);
disp(sprintf(' %% study index %.2d: number of (0) vs (1): (%.5d,%.5d)',ns,n_0,n_1));
tmp_sum = tmp_sum + cast(B_{2+ns},'double');
n_1_(1+ns) = n_1;
end;%for ns=1:10;
ns=0;
n_0 = sum(tmp_sum==1); n_1 = sum(tmp_sum==0);
disp(sprintf(' %% study index %.2d: number of (0) vs (1): (%.5d,%.5d)',ns,n_0,n_1));
n_1_(1+ns) = n_1;
figure;
bar(0:10,n_1_); xlim([-1,11]); xlabel('study index'); ylabel('number'); title('number per study');
set(gcf,'Position',1+[1,1,1024,512]);
print('-depsc','./mds_B_study_index_histogram.eps');

npats_A = length(A_{1}); npats_B = length(B_{1}); 
A_names = cell(npats_A,1); for np=1:npats_A; A_names{np} = sprintf('%s%s%s',A_{1}{np},'&',A_{2}{np}); end;%for np=1:npats_A;
B_names = cell(npats_B,1); for np=1:npats_B; B_names{np} = sprintf('%s%s%s',B_{1}{np},'&',B_{2}{np}); end;%for np=1:npats_B;
assert(length(intersect(A_names,B_names))==0);
npats_mds = npats_A + npats_B;
mds = zeros(npats_mds,20);

rij_cas = zeros(npats_mds,1); rij_con = zeros(npats_mds,1); rij_fam = zeros(npats_mds,1);
for np=1:npats_A; 
rij_cas(np) = strcmp(A_names{np}(1:3),'cas'); 
rij_con(np) = strcmp(A_names{np}(1:3),'con'); 
rij_fam(np) = strcmp(A_names{np}(1:3),'fam'); 
end;%for np=1:npats_A;
for np=1:npats_B; 
rij_cas(npats_A+np) = strcmp(B_names{np}(1:3),'cas'); 
rij_con(npats_A+np) = strcmp(B_names{np}(1:3),'con'); 
rij_fam(npats_A+np) = strcmp(B_names{np}(1:3),'fam'); 
end;%for np=1:npats_B;

rng(0);
for np=1:npats_A;
for ip=1:20; mds(np,ip) = A_{3+ip}(np); end;% for ip=1:20;
end;%for np=1:npats_A;
mds_avg = mean(mds(1:npats_A,:));
mds_std = std(mds(1:npats_A,:));
for np=1:npats_B;
for ip=1:5; mds(npats_A+np,ip) = B_{12+ip}(np); end;% for ip=1:5;
mds(npats_A+np,6:20) = mds_avg(6:20) + mds_std(6:20).*randn(1,15);
end;%for np=1:npats_B;

figure;
v1=1;v2=2;v3=3; %v1=4;v2=5;v3=6; 
subplot(1,2,1);
hold on;
tmp_ij = intersect(find(rij_cas),1:npats_A);
plot3(mds(tmp_ij,v1),mds(tmp_ij,v2),mds(tmp_ij,v3),'r.')
tmp_ij = intersect(find(rij_con),1:npats_A);
plot3(mds(tmp_ij,v1),mds(tmp_ij,v2),mds(tmp_ij,v3),'b.'); 
tmp_ij = intersect(find(rij_fam),1:npats_A);
plot3(mds(tmp_ij,v1),mds(tmp_ij,v2),mds(tmp_ij,v3),'k.'); 
axis vis3d;
hold off;
subplot(1,2,2);
hold on;
tmp_ij = intersect(find(rij_cas),npats_A + [1:npats_B]);
plot3(mds(tmp_ij,v1),mds(tmp_ij,v2),mds(tmp_ij,v3),'r.')
tmp_ij = intersect(find(rij_con),npats_A + [1:npats_B]);
plot3(mds(tmp_ij,v1),mds(tmp_ij,v2),mds(tmp_ij,v3),'b.'); 
tmp_ij = intersect(find(rij_fam),npats_A + [1:npats_B]);
plot3(mds(tmp_ij,v1),mds(tmp_ij,v2),mds(tmp_ij,v3),'k.'); 
axis vis3d;
hold off;

figure;
v1=1;v2=2;v3=3; %v1=4;v2=5;v3=6; 
hold on;
tmp_ij = 1:npats_A;
plot3(mds(tmp_ij,v1),mds(tmp_ij,v2),mds(tmp_ij,v3),'r.')
tmp_ij = npats_A + [1:npats_B];
plot3(mds(tmp_ij,v1),mds(tmp_ij,v2),mds(tmp_ij,v3),'b.'); 
hold off;
axis vis3d;


mds_names = [A_names;B_names];
save('/data/rangan/dir_bcc/dir_PGC_20180304/mds.mat','npats_A','npats_B','npats_mds','mds','rij_cas','rij_con','rij_fam','mds_names');
