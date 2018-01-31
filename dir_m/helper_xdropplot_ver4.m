function helper_xdropplot_ver4(prefix,infix,p_threshold,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,nshuf,nstd)
% plots rank-ordering for trace ; here mds_srt should be the sorted array of mds-components associated with the rows of A_full_n.b16 ;
%{
  %072217; 
  prefix = 'PGC_cl3_maf01'; infix = 'dex_ADZSZDA'; p_threshold = 0.01; rev_flag=0; nmds = 20; mds_used = [1,2]; gamma = 0.004; B_MLT = 34; Ireq = 0; nshuf = 1+64; nstd=0;
  %081517; 
  prefix = 'PGC_cl3_maf01'; infix = 'lak_ADZSZDA'; p_threshold = 0.01; rev_flag=0; nmds = 20; mds_used = [1,2]; gamma = 0.004; B_MLT = 34; Ireq = 0; nshuf = 1+64; nstd=0;
  helper_xdropplot_ver4(prefix,infix,p_threshold,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,nshuf,nstd) ;  
  %081517; checking older runs (not sure if these are completely bug free) ;
  prefix = 'PGC_cl3_maf01'; infix = ''; p_threshold = 0.25; rev_flag=0; nmds = 20; mds_used = [1,2]; gamma = 0.004; B_MLT = 34; Ireq = 0; nshuf = 1+64; nstd=0;
  helper_xdropplot_ver4(prefix,infix,p_threshold,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,nshuf,nstd) ;
  %081517; checking dexcluster runs again ;
  prefix = 'PGC_cl3_maf01'; infix = 'dex_ADZSZDA'; p_threshold = 0.25; rev_flag=0; nmds = 20; mds_used = [1,2]; gamma = 0.004; B_MLT = 34; Ireq = 0; nshuf = 1+64; nstd=0; helper_xdropplot_ver4(prefix,infix,p_threshold,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,nshuf,nstd);
  %081517; double checking dexcluster X runs ;
  prefix = 'PGC_cl3_maf01'; infix = 'dex_ADZSZDA'; p_threshold = 0.10; rev_flag=1; nmds = 20; mds_used = [1,2]; gamma = 0.004; B_MLT = 34; Ireq = 0; nshuf = 1+64; nstd=0; helper_xdropplot_ver4(prefix,infix,p_threshold,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,nshuf,nstd);
  %081517; double checking dexcluster X runs with original gamma value of 0.05;
  prefix = 'PGC_cl3_maf01'; infix = 'dex_ADZSZDA'; p_threshold = 0.10; rev_flag=1; nmds = 20; mds_used = [1,2]; gamma = 0.050; B_MLT = 34; Ireq = 0; nshuf = 1+64; nstd=0; helper_xdropplot_ver4(prefix,infix,p_threshold,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,nshuf,nstd);
  %081517; double checking dexcluster X runs with original gamma value of 0.05 on older PGC_cl3 data set (no maf01);
  prefix = 'PGC_cl3'; infix = 'dex_ADZSZDA'; p_threshold = -1; rev_flag=1; nmds = 20; mds_used = [1,2]; gamma = 0.050; B_MLT = 34; Ireq = 0; nshuf = 1+64; nstd=0; helper_xdropplot_ver4(prefix,infix,p_threshold,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,nshuf,nstd);
  %081917; rechecking PGC_cl3_maf01_lak_ADZSZDA_p10_X_m3_g004_B34_n0_sXXXX.mat from rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/dir_PGC_cl3_maf01/ ;
  prefix = 'PGC_cl3_maf01'; infix = 'lak_ADZSZDA'; p_threshold = 0.10; rev_flag=1; nmds = 20; mds_used = [1,2]; gamma = 0.004; B_MLT = 34; Ireq = 0; nshuf = 1+16; nstd=0; helper_xdropplot_ver4(prefix,infix,p_threshold,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,nshuf,nstd);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % extract rdrop and cdrop for example-2 in manuscript ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  fname_1 = '/data/rangan/dir_bcc/dir_code_022316/dir_PGC_040116/PGC_cl3_D_m3_g05_B33_n0_sXXXX.mat'; % Note that 'D' tag was incorrectly applied during original lakcluster_ver17 runs for manuscript ;
  clear xa_trn xb_trn; load(fname_1,'xa_trn__','xb_trn__'); nf=1;
  xa_trn = xa_trn__{1}{1}{1}{nf}; rdrop_a = xa_trn(find(xa_trn(:,1)>-1),1)+1; cdrop_a = xa_trn(find(xa_trn(:,2)>-1),2)+1;
  xb_trn = xb_trn__{1}{1}{1}{nf}; rdrop_b = xb_trn(find(xb_trn(:,1)>-1),1)+1; cdrop_b = xb_trn(find(xb_trn(:,2)>-1),2)+1;
  ra_1 = rdrop_a; rb_1 = rdrop_b; ca_1 = cdrop_a; cb_1 = cdrop_b;
  fname_2 = '/data/rangan/dir_bcc/dir_code_081717/dir_m/dir_PGC_cl3_maf01/PGC_cl3_maf01_lak_ADZSZDA_p10_X_m3_g004_B34_n0_sXXXX.mat'; % run with lakcluster_ver18 ;
  clear xa_trn xb_trn; load(fname_2,'xa_trn__','xb_trn__'); nf=1;
  xa_trn = xa_trn__{1}{1}{1}{nf}; rdrop_a = xa_trn(find(xa_trn(:,1)>-1),1)+1; cdrop_a = xa_trn(find(xa_trn(:,2)>-1),2)+1;
  xb_trn = xb_trn__{1}{1}{1}{nf}; rdrop_b = xb_trn(find(xb_trn(:,1)>-1),1)+1; cdrop_b = xb_trn(find(xb_trn(:,2)>-1),2)+1;
  ra_2 = rdrop_a; rb_2 = rdrop_b; ca_2 = cdrop_a; cb_2 = cdrop_b;
  % calculate overlap intersection of patients starting from the beginning ;
  n_p = length(rb_2); PS_p = cintersect(rb_1,rb_2)./transpose(1:n_p); 
  iteration_max=128; PS__ = zeros(n_p,iteration_max); for ni=1:iteration_max; PS__(:,ni) = cintersect(rb_1,rb_2(randperm(n_p)))./transpose(1:n_p); end;
  % calculate overlap intersection of patients starting from the end ;
  n_p = length(rb_2); PF_p = cintersect(rb_1(end:-1:1),rb_2(end:-1:1))./transpose(1:n_p); 
  iteration_max=128; PF__ = zeros(n_p,iteration_max); for ni=1:iteration_max; PF__(:,ni) = cintersect(rb_1(end:-1:1),rb_2(randperm(n_p)))./transpose(1:n_p); end;
  % calculate overlap intersection of snps starting from the beginning ;
  n_g = length(cb_2); GS_p = cintersect(cb_1,cb_2)./transpose(1:n_g); 
  iteration_max=32; GS__ = zeros(n_g,iteration_max); for ni=1:iteration_max; GS__(:,ni) = cintersect(cb_1,cb_2(randperm(n_g)))./transpose(1:n_g); end;
  % calculate overlap intersection of snps starting from the end ;
  n_g = length(cb_2); GF_p = cintersect(cb_1(end:-1:1),cb_2(end:-1:1))./transpose(1:n_g); 
  iteration_max=32; GF__ = zeros(n_g,iteration_max); for ni=1:iteration_max; GF__(:,ni) = cintersect(cb_1(end:-1:1),cb_2(randperm(n_g)))./transpose(1:n_g); end;
  subplot(2,3,1); hold on; plot(1:n_p,PS__,'b-','LineWidth',0.5); plot(1:n_p,PS_p,'r-','LineWidth',2); hold off; axis([1,n_p,0,1]); xlabel('patients ejected'); ylabel('overlap'); title('patient intersection from beginning');
  subplot(2,3,2); hold on; plot(1:n_p,PF__,'b-','LineWidth',0.5); plot(1:n_p,PF_p,'r-','LineWidth',2); hold off; axis([1,n_p,0,1]); xlabel('patients retained'); ylabel('overlap'); title('patient intersection from end');
  subplot(2,3,3); hold on; plot(1:n_p,PF__,'b-','LineWidth',0.5); plot(1:n_p,PF_p,'r-','LineWidth',2); line(115*[1,1],[0,1],'LineStyle',':'); hold off; axis([1,115*2,0,1]); xlabel('patients retained'); ylabel('overlap'); title('patient intersection from end');
  subplot(2,3,4); hold on; plot(1:n_g,GS__,'b-','LineWidth',0.5); plot(1:n_g,GS_p,'r-','LineWidth',2); hold off; axis([1,n_g,0,1]); xlabel('alleles ejected'); ylabel('overlap'); title('allele intersection from beginning');
  subplot(2,3,5); hold on; plot(1:n_g,GF__,'b-','LineWidth',0.5); plot(1:n_g,GF_p,'r-','LineWidth',2); hold off; axis([1,n_g,0,1]); xlabel('alleles retained'); ylabel('overlap'); title('allele intersection from end');
  subplot(2,3,6); hold on; plot(1:n_g,GF__,'b-','LineWidth',0.5); plot(1:n_g,GF_p,'r-','LineWidth',2); line(706*[1,1],[0,1],'LineStyle',':'); hold off; axis([1,706*2,0,1]); xlabel('alleles retained'); ylabel('overlap'); title('allele intersection from end');
  set(gcf,'Position',1+[0,0,1024+512,768]);
  print('-depsc','PGC_cl3_maf01_lak_ADZSZDA_p10_X_m3_g05_B32_vs_g004_B34_ranking.eps');

  %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% load mat file ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
p_str = sprintf('p%.2d',floor(100*p_threshold));
if (rev_flag==1); rev_str = 'X'; else rev_str = 'D'; end;
mds_tmp = zeros(1,nmds); mds_tmp(mds_used)=1; mds_code = dot(2.^[0:nmds-1],mds_tmp);
gamma_d = floor(gamma*1000);
if (length(infix)>0);
if (p_threshold>=0); output_string = sprintf('%s_%s_%s_m%d_g%.3d_B%.2d_n%d_sXXXX',infix,p_str,rev_str,mds_code,gamma_d,B_MLT,Ireq); end;
if (p_threshold< 0); output_string = sprintf('%s_%s_m%d_g%.3d_B%.2d_n%d_sXXXX',infix,rev_str,mds_code,gamma_d,B_MLT,Ireq); end;
end;%if (length(infix)>0);
if (length(infix)==0);
if (p_threshold>=0); output_string = sprintf('%s_%s_m%d_g%.3d_B%.2d_n%d_sXXXX',p_str,rev_str,mds_code,gamma_d,B_MLT,Ireq); end;
if (p_threshold< 0); output_string = sprintf('%s_m%d_g%.3d_B%.2d_n%d_sXXXX',rev_str,mds_code,gamma_d,B_MLT,Ireq); end;
end;%if (length(infix)==0);
dir_pre2 = sprintf('/data/rangan/dir_bcc/dir_code_081717/dir_m/dir_%s',prefix);
fname = sprintf('%s/%s_%s.mat',dir_pre2,prefix,output_string),;
load(fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% load metadata ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

dir_pre1 = '/data/rangan/dir_bcc/dir_code_022316/dir_PGC_040116';
load(sprintf('%s/mds.mat',dir_pre1));
fid = fopen(sprintf('%s/%s_fam.ext',dir_pre1,prefix)); A = textscan(fid,'%s %s %s %s %d %d %s','headerlines',0); fclose(fid);
npats_A = length(A{1}); A_names = cell(npats_A,1); for np=1:npats_A; A_names{np} = sprintf('%s%s%s',A{1}{np},'&',A{2}{np}); end;%for np=1:npats_A;
[tmp,iA,iB] = intersect(A_names,mds_names,'stable'); if (length(iA)<length(A_names)); disp('Warning! A_names not a subset of mds_names'); end; mds_srt = mds(iB,:);

fid = fopen(sprintf('%s/PGC_BIP32b-BD1.20160612.pheno',dir_pre1)); BD1 = textscan(fid,'%s %s %d','headerlines',1); fclose(fid);
npats_BD1 = length(BD1{1}); BD1_names = cell(npats_BD1,1); for np=1:npats_BD1; BD1_names{np} = sprintf('%s%s%s',BD1{1}{np},'&',BD1{2}{np}); end;%for np=1:npats_BD1;
BD1_pheno = BD1{3};
BD1_status = zeros(length(A_names),1); [tmp,iA,iB] = intersect(A_names,BD1_names,'stable'); BD1_status(iA) = BD1_pheno(iB);
fid = fopen(sprintf('%s/PGC_BIP32b-BD2.20160612.pheno',dir_pre1)); BD2 = textscan(fid,'%s %s %d','headerlines',1); fclose(fid);
npats_BD2 = length(BD2{1}); BD2_names = cell(npats_BD2,1); for np=1:npats_BD2; BD2_names{np} = sprintf('%s%s%s',BD2{1}{np},'&',BD2{2}{np}); end;%for np=1:npats_BD2;
BD2_pheno = BD2{3};
BD2_status = zeros(length(A_names),1); [tmp,iA,iB] = intersect(A_names,BD2_names,'stable'); BD2_status(iA) = BD2_pheno(iB);

bim_ext_fname = sprintf('%s/%s_bim.ext',dir_pre1,prefix);
fid = fopen(bim_ext_fname); bim_ext = textscan(fid,'%d\t%s\t%d\t%d\t%c\t%c\t%s\t%f\t%f\t%f\n','headerlines',0); fclose(fid);
bim_nsnps = length(bim_ext{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% unpack traces ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
tmp_trn = tr_trn__{1}{1}{1}{1};
[max_iteration,ncols] = size(tmp_trn);
if (max_iteration==1 | ncols~=6); disp(sprintf(' %% Warning! max_iteration %d ncols %d in first trace',max_iteration,ncols)); end;
ni_ = tmp_trn(:,1);
r_rem = (tmp_trn(:,2));
c_rem = (tmp_trn(:,3));
t_rem = (tmp_trn(:,6)); t_ret = find(t_rem>=nstd);

finished_trn__ = zeros(nshuf,1); 
for ns=1:nshuf;
try;
tmp_trn = tr_trn__{1}{1}{1}{ns};
[rtmp,ctmp] = size(tmp_trn);
if (rtmp==max_iteration & ctmp==6); finished_trn__(ns)=1; else finished_trn__(ns)=0; end;
catch; end;%try;
end;%for ns=1:nshuf;

f_ij = find(finished_trn__);
n_fin = length(f_ij); disp(sprintf(' %% found n_fin %d',n_fin));

ZR = zeros(n_fin,max_iteration);
ZC = zeros(n_fin,max_iteration);

for nf=1:n_fin;
ns = f_ij(nf); tmp_trn = tr_trn__{1}{1}{1}{ns};
ZR(nf,:) = transpose(tmp_trn(:,4));
ZC(nf,:) = transpose(tmp_trn(:,5));
end;%for nf=1:n_fin;

PR = zeros(n_fin,max_iteration);
PC = zeros(n_fin,max_iteration);

for ni=1:max_iteration;
[tmp,tmp_ij] = sort(ZR(:,ni),'ascend'); [tmp,PR(:,ni)] = sort(tmp_ij,'ascend');
[tmp,tmp_ij] = sort(ZC(:,ni),'ascend'); [tmp,PC(:,ni)] = sort(tmp_ij,'ascend');
end;%for ni=1:max_iteration;

ZR = transpose(ZR); ZC = transpose(ZC);
PR = transpose(PR); PC = transpose(PC);
%LR = log10(abs(ZR)); LC = log10(abs(ZC));
%LRm = mean(LR,2); LRs = std(LR,[],2); LCm = mean(LC,2); LCs = std(LC,[],2);
%NR = (LR - repmat(LRm,1,n_fin))./repmat(LRs,1,n_fin); NC = (LC - repmat(LCm,1,n_fin))./repmat(LCs,1,n_fin);
ZRm = mean(ZR,2); ZRs = std(ZR,[],2); ZCm = mean(ZC,2); ZCs = std(ZC,[],2);
NR = (ZR - repmat(ZRm,1,n_fin))./repmat(ZRs,1,n_fin); NC = (ZC - repmat(ZCm,1,n_fin))./repmat(ZCs,1,n_fin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
nf=1;
xa_trn = xa_trn__{1}{1}{1}{nf}; rdrop_a = xa_trn(find(xa_trn(:,1)>-1),1)+1; cdrop_a = xa_trn(find(xa_trn(:,2)>-1),2)+1;
xb_trn = xb_trn__{1}{1}{1}{nf}; rdrop_b = xb_trn(find(xb_trn(:,1)>-1),1)+1; cdrop_b = xb_trn(find(xb_trn(:,2)>-1),2)+1;
disp(sprintf(' %% figure 1: out_trace results; %d pats -x- %d snps ; %d iterations',max(rdrop_b),max(cdrop_b),max(t_ret)));
subplot(2,4,1); hold on; 
plot(max(r_rem)-r_rem(t_ret),NR(t_ret,2:end),'b-'); 
plot(max(r_rem)-r_rem(t_ret),NR(t_ret,1),'ro-'); 
xlabel('rows eliminated'); ylabel('(NR)'); xlim([0,max(max(r_rem)-r_rem(t_ret))+1]);
subplot(2,4,2); hold on; 
plot(max(r_rem)-r_rem(t_ret),PR(t_ret,2:end),'b-');  
plot(max(r_rem)-r_rem(t_ret),PR(t_ret,1),'ro-'); 
xlabel('rows eliminated'); ylabel('PR'); xlim([0,max(max(r_rem)-r_rem(t_ret))+1]);ylim([0,n_fin+1]);
subplot(2,4,3); hold on; 
plot(max(t_ret)-t_ret(end:-1:1),NR(t_ret,2:end),'b-'); 
plot(max(t_ret)-t_ret(end:-1:1),NR(t_ret,1),'ro-'); 
xlabel('row iteration number'); ylabel('(NR)'); xlim([0,max((t_ret))+1]);
subplot(2,4,4); hold on; 
plot(max(t_ret)-t_ret(end:-1:1),PR(t_ret,2:end),'b-'); 
plot(max(t_ret)-t_ret(end:-1:1),PR(t_ret,1),'ro-'); 
xlabel('row iteration number'); ylabel('PR'); xlim([0,max((t_ret))+1]);ylim([0,n_fin+1]);
subplot(2,4,5); hold on; 
plot(max(c_rem)-c_rem(t_ret),NC(t_ret,2:end),'b-'); 
plot(max(c_rem)-c_rem(t_ret),NC(t_ret,1),'ro-'); 
xlabel('cols eliminated'); ylabel('(NC)'); xlim([0,max(max(c_rem)-c_rem(t_ret))+1]);
subplot(2,4,6); hold on; 
plot(max(c_rem)-c_rem(t_ret),PC(t_ret,2:end),'b-'); 
plot(max(c_rem)-c_rem(t_ret),PC(t_ret,1),'ro-'); 
xlabel('cols eliminated'); ylabel('PC'); xlim([0,max(max(c_rem)-c_rem(t_ret))+1]);ylim([0,n_fin+1]);
subplot(2,4,7); hold on; 
plot(max(t_ret)-t_ret(end:-1:1),NC(t_ret,2:end),'b-'); 
plot(max(t_ret)-t_ret(end:-1:1),NC(t_ret,1),'ro-'); 
xlabel('col iteration number'); ylabel('(NC)'); xlim([0,max((t_ret))+1]);
subplot(2,4,8); hold on; 
plot(max(t_ret)-t_ret(end:-1:1),PC(t_ret,2:end),'b-'); 
plot(max(t_ret)-t_ret(end:-1:1),PC(t_ret,1),'ro-'); 
xlabel('col iteration number'); ylabel('PC'); xlim([0,max((t_ret))+1]);ylim([0,n_fin+1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
aij = 1:min(100,length(r_rem)-1); dr = diff(r_rem); dc = diff(c_rem); dr2 = transpose(dr(aij))/sum(dr(aij)); dc2 = transpose(dc(aij))/sum(dc(aij));
mij = 30:min(100,length(r_rem)-1);
subplot(2,4,1); title('max-vs-mean-iter NR');
hold on;
plot(max(NR(mij,2:end),[],1),mean(NR(aij,2:end),1),'bo');
plot(max(NR(mij,1    ),[],1),mean(NR(aij,1    ),1),'ro');
hold off;
subplot(2,4,5); title('max-vs-mean-row NR');
hold on;
plot(max(NR(mij,2:end),[],1),dr2*NR(aij,2:end),'bx');
plot(max(NR(mij,1    ),[],1),dr2*NR(aij,1    ),'rx');
hold off;
subplot(2,4,2); title('max-vs-mean-iter NC');
hold on;
plot(max(NC(mij,2:end),[],1),mean(NC(aij,2:end),1),'bo');
plot(max(NC(mij,1    ),[],1),mean(NC(aij,1    ),1),'ro');
hold off;
subplot(2,4,6); title('max-vs-mean-col NC');
hold on;
plot(max(NC(mij,2:end),[],1),dc2*NC(aij,2:end),'bx');
plot(max(NC(mij,1    ),[],1),dc2*NC(aij,1    ),'rx');
hold off;
subplot(2,4,3); title('max-vs-mean-iter PR');
hold on;
plot(max(PR(mij,2:end),[],1),mean(PR(aij,2:end),1),'bo');
plot(max(PR(mij,1    ),[],1),mean(PR(aij,1    ),1),'ro');
hold off;
subplot(2,4,7); title('max-vs-mean-row PR');
hold on;
plot(max(PR(mij,2:end),[],1),dr2*PR(aij,2:end),'bx');
plot(max(PR(mij,1    ),[],1),dr2*PR(aij,1    ),'rx');
hold off;
subplot(2,4,4); title('max-vs-mean-iter PC');
hold on;
plot(max(PC(mij,2:end),[],1),mean(PC(aij,2:end),1),'bo');
plot(max(PC(mij,1    ),[],1),mean(PC(aij,1    ),1),'ro');
hold off;
subplot(2,4,8); title('max-vs-mean-col PC');
hold on;
plot(max(PC(mij,2:end),[],1),dc2*PC(aij,2:end),'bx');
plot(max(PC(mij,1    ),[],1),dc2*PC(aij,1    ),'rx');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

figure;
nf=1;
xa_trn = xa_trn__{1}{1}{1}{nf}; rdrop_a = xa_trn(find(xa_trn(:,1)>-1),1)+1; cdrop_a = xa_trn(find(xa_trn(:,2)>-1),2)+1;
xb_trn = xb_trn__{1}{1}{1}{nf}; rdrop_b = xb_trn(find(xb_trn(:,1)>-1),1)+1; cdrop_b = xb_trn(find(xb_trn(:,2)>-1),2)+1;
disp(sprintf(' %% figure 1: out_trace results; %d pats -x- %d snps ; %d iterations',max(rdrop_b),max(cdrop_b),max(t_ret)));
subplot(2,4,1); plot(max(r_rem)-r_rem(t_ret),ZR(t_ret,1),'ro-',max(r_rem)-r_rem(t_ret),ZR(t_ret,2:end),'b-'); xlabel('rows eliminated'); ylabel('(ZR)'); xlim([0,max(max(r_rem)-r_rem(t_ret))+1]);
subplot(2,4,2); plot(max(r_rem)-r_rem(t_ret),PR(t_ret,1),'ro-',max(r_rem)-r_rem(t_ret),PR(t_ret,2:end),'b-');  xlabel('rows eliminated'); ylabel('PR'); xlim([0,max(max(r_rem)-r_rem(t_ret))+1]);ylim([0,n_fin+1]);
subplot(2,4,3); plot(max(t_ret)-t_ret(end:-1:1),ZR(t_ret,1),'ro-',max(t_ret)-t_ret(end:-1:1),ZR(t_ret,2:end),'b-'); xlabel('row iteration number'); ylabel('(ZR)'); xlim([0,max((t_ret))+1]);
subplot(2,4,4); plot(max(t_ret)-t_ret(end:-1:1),PR(t_ret,1),'ro-',max(t_ret)-t_ret(end:-1:1),PR(t_ret,2:end),'b-'); xlabel('row iteration number'); ylabel('PR'); xlim([0,max((t_ret))+1]);ylim([0,n_fin+1]);
subplot(2,4,5); plot(max(c_rem)-c_rem(t_ret),ZC(t_ret,1),'ro-',max(c_rem)-c_rem(t_ret),ZC(t_ret,2:end),'b-'); xlabel('cols eliminated'); ylabel('(ZC)'); xlim([0,max(max(c_rem)-c_rem(t_ret))+1]);
subplot(2,4,6); plot(max(c_rem)-c_rem(t_ret),PC(t_ret,1),'ro-',max(c_rem)-c_rem(t_ret),PC(t_ret,2:end),'b-'); xlabel('cols eliminated'); ylabel('PC'); xlim([0,max(max(c_rem)-c_rem(t_ret))+1]);ylim([0,n_fin+1]);
subplot(2,4,7); plot(max(t_ret)-t_ret(end:-1:1),ZC(t_ret,1),'ro-',max(t_ret)-t_ret(end:-1:1),ZC(t_ret,2:end),'b-'); xlabel('col iteration number'); ylabel('(ZC)'); xlim([0,max((t_ret))+1]);
subplot(2,4,8); plot(max(t_ret)-t_ret(end:-1:1),PC(t_ret,1),'ro-',max(t_ret)-t_ret(end:-1:1),PC(t_ret,2:end),'b-'); xlabel('col iteration number'); ylabel('PC'); xlim([0,max((t_ret))+1]);ylim([0,n_fin+1]);

figure;
aij = 1:100; dr = diff(r_rem); dc = diff(c_rem); dr2 = transpose(dr(aij))/sum(dr(aij)); dc2 = transpose(dc(aij))/sum(dc(aij));
mij = 30:100;
subplot(2,4,1); title('max-vs-mean ZR');
hold on;
plot(max(ZR(mij,2:end),[],1),mean(ZR(aij,2:end),1),'bo');
plot(max(ZR(mij,1    ),[],1),mean(ZR(aij,1    ),1),'ro');
hold off;
subplot(2,4,5); title('max-vs-mean ZR');
hold on;
plot(max(ZR(mij,2:end),[],1),dr2*ZR(aij,2:end),'bx');
plot(max(ZR(mij,1    ),[],1),dr2*ZR(aij,1    ),'rx');
hold off;
subplot(2,4,2); title('max-vs-mean ZC');
hold on;
plot(max(ZC(mij,2:end),[],1),mean(ZC(aij,2:end),1),'bo');
plot(max(ZC(mij,1    ),[],1),mean(ZC(aij,1    ),1),'ro');
hold off;
subplot(2,4,6); title('max-vs-mean ZC');
hold on;
plot(max(ZC(mij,2:end),[],1),dc2*ZC(aij,2:end),'bx');
plot(max(ZC(mij,1    ),[],1),dc2*ZC(aij,1    ),'rx');
hold off;
subplot(2,4,3); title('max-vs-mean PR');
hold on;
plot(max(PR(mij,2:end),[],1),mean(PR(aij,2:end),1),'bo');
plot(max(PR(mij,1    ),[],1),mean(PR(aij,1    ),1),'ro');
hold off;
subplot(2,4,7); title('max-vs-mean PR');
hold on;
plot(max(PR(mij,2:end),[],1),dr2*PR(aij,2:end),'bx');
plot(max(PR(mij,1    ),[],1),dr2*PR(aij,1    ),'rx');
hold off;
subplot(2,4,4); title('max-vs-mean PC');
hold on;
plot(max(PC(mij,2:end),[],1),mean(PC(aij,2:end),1),'bo');
plot(max(PC(mij,1    ),[],1),mean(PC(aij,1    ),1),'ro');
hold off;
subplot(2,4,8); title('max-vs-mean PC');
hold on;
plot(max(PC(mij,2:end),[],1),dc2*PC(aij,2:end),'bx');
plot(max(PC(mij,1    ),[],1),dc2*PC(aij,1    ),'rx');
hold off;

 %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nf=1;
xa_trn = xa_trn__{1}{1}{1}{nf}; rdrop_a = xa_trn(find(xa_trn(:,1)>-1),1)+1; cdrop_a = xa_trn(find(xa_trn(:,2)>-1),2)+1;
xb_trn = xb_trn__{1}{1}{1}{nf}; rdrop_b = xb_trn(find(xb_trn(:,1)>-1),1)+1; cdrop_b = xb_trn(find(xb_trn(:,2)>-1),2)+1;
rdrop_a = xa_trn(find(xa_trn(:,1)>-1),1)+1; cdrop_a = xa_trn(find(xa_trn(:,2)>-1),2)+1;
rdrop_b = xb_trn(find(xb_trn(:,1)>-1),1)+1; cdrop_b = xb_trn(find(xb_trn(:,2)>-1),2)+1;
figure;cla;
subplot(2,2,1); plot(rdrop_a,'.'); ll=size(mds_srt,1); xlim([1,length(rdrop_a)]);ylim([1,ll]);
subplot(2,2,3); ll=size(mds_srt,1); imagesc(log2(1+adi_hist2d(1:length(rdrop_a),rdrop_a,[1,length(rdrop_a)],256,[1,ll],256))); 
subplot(2,2,2);plot(cdrop_a,'.'); ll=bim_nsnps; xlim([1,length(cdrop_a)]);ylim([1,ll]);
subplot(2,2,4); ll=bim_nsnps; imagesc(log2(1+adi_hist2d(1:length(cdrop_a),cdrop_a,[1,length(cdrop_a)],256,[1,ll],256)));

nf=1;
xa_trn = xa_trn__{1}{1}{1}{nf}; rdrop_a = xa_trn(find(xa_trn(:,1)>-1),1)+1; cdrop_a = xa_trn(find(xa_trn(:,2)>-1),2)+1;
xb_trn = xb_trn__{1}{1}{1}{nf}; rdrop_b = xb_trn(find(xb_trn(:,1)>-1),1)+1; cdrop_b = xb_trn(find(xb_trn(:,2)>-1),2)+1;
rdrop_a = xa_trn(find(xa_trn(:,1)>-1),1)+1; cdrop_a = xa_trn(find(xa_trn(:,2)>-1),2)+1;
rdrop_b = xb_trn(find(xb_trn(:,1)>-1),1)+1; cdrop_b = xb_trn(find(xb_trn(:,2)>-1),2)+1;
figure;cla;
xl = 0.1*[-0.5,1]; yl = 0.065*[-1,1];
v1=1;v2=2;cra = colormap('jet'); ncra=size(cra,1);
for np=1:8;
subplot(2,4,np); rcut = max(1,min(length(rdrop_a),floor(length(rdrop_a)/2.^(np-1)))); rij = rdrop_a(max(1,length(rdrop_a)-rcut):end); 
imagesc(log2(1+adi_hist2d(mds_srt(rij,v1),mds_srt(rij,v2),xl,256,yl,256))); 
title(sprintf('np %d --> npats %d',np,length(rij)));
end;%for np=1:8;
%figure;cla;
%xl = 0.1*[-0.5,1]; yl = 0.065*[-1,1];
%v1=1;v2=2;cra = colormap('jet'); ncra=size(cra,1);
%subplot(1,1,1); rcut = 115; rij = rdrop_a(max(1,end-rcut):end); 
%imagesc(log2(1+adi_hist2d(mds_srt(rij,v1),mds_srt(rij,v2),xl,256,yl,256))); 
%axis off;
%title(sprintf('npats %d',length(rij)));

%{
nf=1;
xa_trn = xa_trn__{1}{1}{1}{nf}; rdrop_a = xa_trn(find(xa_trn(:,1)>-1),1)+1; cdrop_a = xa_trn(find(xa_trn(:,2)>-1),2)+1;
xb_trn = xb_trn__{1}{1}{1}{nf}; rdrop_b = xb_trn(find(xb_trn(:,1)>-1),1)+1; cdrop_b = xb_trn(find(xb_trn(:,2)>-1),2)+1;
figure;cla;
v1=1;v2=2;v3=3;cra = colormap('jet'); ncra=size(cra,1);
hold on;
for nr=1:length(rdrop_a);
nb = max(1,min(ncra,1+floor(ncra*nr/length(rdrop_a))));
plot3(mds_srt(rdrop_a(nr),v1),mds_srt(rdrop_a(nr),v2),mds_srt(rdrop_a(nr),v3),'.','Color',cra(nb,:));
end;%for nr=1:length(rdrop_a);
hold off;
axis vis3d;
 %}
