function bed_to_b16_ver2(cl_num);
%{
  scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_022316/dir_PGC_040116/intersectall.m ./ ;
  scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_022316/dir_PGC_040116/bed_to_b16_ver0.m ./ ;
  scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_022316/tutorial_binary_uncompress.m ./ ;
  scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_022316/tutorial_binary_compress.m ./ ;
  scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_022316/dir_PGC_040116/dosage_gen_ver3.m ./ ;

  !scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_022316/bed_to_b16_ver2.m ./;
   %}

%clear;
case_tag=2; ctrl_tag=1; nkhr=22; 
bitj=16; bit8=8; bit4=4; snp_and_tag=3; snp__or_tag=2; snp_nnd_tag=0; snp_mss_tag=1; snp_not_tag=5;

ent_cutoff = 0.03; frq_cutoff = 0.01; mss_cutoff = 0.02;

if nargin<1; cl_num = input('cluster number? [1,2,3]'); end;

if cl_num==1;
nstd = 11; ns = 1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_may1_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'ms.bip_butr_eur-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_top8_eur_sr-qc2'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'ms.bip_uktr_eur-qc'; ns = ns+1;
dir_base{ns} = 'v1_boma'; fpre{ns} = 'bip_bmg2_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_ume4_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_ucl2_eur_sr-qc2'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_usc2_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1_boma'; fpre{ns} = 'bip_rom3_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_swei_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1_boma'; fpre{ns} = 'bip_hal2_eur_sr-qc'; ns = ns+1;
end;% cl1;
if cl_num==2;
nstd = 7; ns = 1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_gsk1_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_mich_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_bonn_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1_boma'; fpre{ns} = 'bip_bmau_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1_boma'; fpre{ns} = 'bip_bmg3_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_fran_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1_boma'; fpre{ns} = 'bip_bmpo_eur_sr-qc'; ns = ns+1;
end;%cl2;
if cl_num==3;
nstd = 10; ns = 1;
dir_base{ns} = 'v1_gain'; fpre{ns} = 'bip_gain_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_dub1_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_top7_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_swa2_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_fat2_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_wtcc_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_edi1_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_uclo_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_stp1_eur_sr-qc'; ns = ns+1;
dir_base{ns} = 'v1'; fpre{ns} = 'bip_st2c_eur_sr-qc'; ns = ns+1;
end;%cl3;

tmp = sprintf('PGC_cl%d_maf%.2d',cl_num,floor(100*frq_cutoff));
snp_pre = sprintf('%s_',tmp);
dir_out = sprintf('%s/dir_%s',pwd,tmp); if ~exist(dir_out,'dir'); mkdir(dir_out); end;
disp(sprintf('%% snp_pre %s\n%% dir_out %s\n%% tmp_char %s/%s\n',snp_pre,dir_out,dir_out,snp_pre));

% fix later ;
% nsnps_crop = 1000 + [1:nstd];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Loading mds.mat ;'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('mds.mat');

for ns=1:nstd;
disp(sprintf(' %% ns %d: %s/%s',ns,dir_base{ns},fpre{ns}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Reading fam file (patients) ;'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnfam{ns} = sprintf('%s/%s.fam',dir_base{ns},fpre{ns});
%[status,npats_tmp] = system(sprintf('wc -l %s | cut -f 1 -d '' ''',fnfam{ns})); npats(ns) = str2num(npats_tmp); 
%tmp_rem = mod(npats(ns),bit4); tmp_quo = (npats(ns)-tmp_rem)/bit4;
%npats_rup(ns) = tmp_quo + (tmp_rem>0);
fpfam{ns} = fopen(fnfam{ns}); cafam{ns} = textscan(fpfam{ns},'%s %s %s %s %d %d'); fclose(fpfam{ns});
npats(ns) = length(cafam{ns}{1}); 
tmp_rem = mod(npats(ns),bit4); tmp_quo = (npats(ns)-tmp_rem)/bit4;
npats_rup(ns) = tmp_quo + (tmp_rem>0);
rij_case_{ns} = find(cafam{ns}{6}==case_tag);rij_ctrl_{ns} = find(cafam{ns}{6}==ctrl_tag);
disp(sprintf('%% npats %d (%d cases, %d ctrls) --> %d',npats(ns),length(rij_case_{ns}),length(rij_ctrl_{ns}),npats_rup(ns)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Reading bim file (snp names) ;'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnbim{ns} = sprintf('%s/%s.bim',dir_base{ns},fpre{ns});
%[status,nsnps_tmp] = system(sprintf('wc -l %s | cut -f 1 -d '' ''',fnbim{ns})); nsnps(ns) = str2num(nsnps_tmp);
fpbim{ns} = fopen(fnbim{ns}); cabim{ns} = textscan(fpbim{ns},'%d %s %d %d %c %c'); fclose(fpbim{ns});
nsnps(ns) = length(cabim{ns}{1});
% fix later ;
% nsnps(ns) = nsnps_crop(ns); for nl=1:6; cabim{ns}{nl} = cabim{ns}{nl}(1:nsnps(ns),:); end;
clear tmp; for nk=1:nkhr; tmp(nk) = length(find(cabim{ns}{1}==nk)); end;%for nk=1:nkhr;
khrdist_{ns} = tmp;
disp(sprintf('%% nsnps %d --> [%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d]',nsnps(ns),khrdist_{ns}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Reading bed file (binarized snp data) ;'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnbed{ns} = sprintf('%s/%s.bed',dir_base{ns},fpre{ns});
fpbed{ns} = fopen(fnbed{ns}); 
key1 = fread(fpbed{ns},1,'uint8'); key2 = fread(fpbed{ns},1,'uint8'); key3 = fread(fpbed{ns},1,'uint8'); assert(key1==108); assert(key2==27); assert(key3==1);
ver_flag=0;
if ver_flag==0;
disp(sprintf(' %% version 0: simultaneous read'));tic;
n_and{ns} = zeros(1,nsnps(ns));n__or{ns} = zeros(1,nsnps(ns));n_nnd{ns} = zeros(1,nsnps(ns));n_mss{ns} = zeros(1,nsnps(ns));
m_and{ns} = zeros(1,4*npats_rup(ns));m__or{ns} = zeros(1,4*npats_rup(ns));m_nnd{ns} = zeros(1,4*npats_rup(ns));m_mss{ns} = zeros(1,4*npats_rup(ns));
uchar_array = transpose(fread(fpbed{ns},[npats_rup(ns),nsnps(ns)],'uint8=>uint8'));
tmp_t = toc; disp(sprintf(' %% simultaneous read %0.2f',tmp_t)); tic;
disp(sprintf(' %% assuming more than 8G of ram... '));
uchar_00000011 = mod(uchar_array,  4);
uchar_00001111 = mod(uchar_array, 16); uchar_00001100 = idivide(uchar_00001111 - uchar_00000011, 4) ; 
uchar_00111111 = mod(uchar_array, 64); uchar_00110000 = idivide(uchar_00111111 - uchar_00001111,16) ; clear uchar_00001111;
                                       uchar_11000000 = idivide(uchar_array    - uchar_00111111,64) ; clear uchar_00111111;
clear uchar_array;
tmp_t = toc; disp(sprintf(' %% mod %0.2f',tmp_t)); tic;
stride1 = 1:4:4*npats_rup(ns); tmp_ij = find(stride1>npats(ns)); uchar_00000011(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride2 = 2:4:4*npats_rup(ns); tmp_ij = find(stride2>npats(ns)); uchar_00001100(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride3 = 3:4:4*npats_rup(ns); tmp_ij = find(stride3>npats(ns)); uchar_00110000(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride4 = 4:4:4*npats_rup(ns); tmp_ij = find(stride4>npats(ns)); uchar_11000000(:,tmp_ij)=snp_not_tag; clear tmp_ij;
clear uchar_tmp;
uchar_tmp = uchar_00000011==snp_and_tag; sum1_and_00000011 = sum(uchar_tmp,1); sum2_and_00000011 = sum(uchar_tmp,2);
uchar_tmp = uchar_00000011==snp__or_tag; sum1__or_00000011 = sum(uchar_tmp,1); sum2__or_00000011 = sum(uchar_tmp,2);
uchar_tmp = uchar_00000011==snp_nnd_tag; sum1_nnd_00000011 = sum(uchar_tmp,1); sum2_nnd_00000011 = sum(uchar_tmp,2);
uchar_tmp = uchar_00000011==snp_mss_tag; sum1_mss_00000011 = sum(uchar_tmp,1); sum2_mss_00000011 = sum(uchar_tmp,2);
clear uchar_tmp;
uchar_tmp = uchar_00001100==snp_and_tag; sum1_and_00001100 = sum(uchar_tmp,1); sum2_and_00001100 = sum(uchar_tmp,2);
uchar_tmp = uchar_00001100==snp__or_tag; sum1__or_00001100 = sum(uchar_tmp,1); sum2__or_00001100 = sum(uchar_tmp,2);
uchar_tmp = uchar_00001100==snp_nnd_tag; sum1_nnd_00001100 = sum(uchar_tmp,1); sum2_nnd_00001100 = sum(uchar_tmp,2);
uchar_tmp = uchar_00001100==snp_mss_tag; sum1_mss_00001100 = sum(uchar_tmp,1); sum2_mss_00001100 = sum(uchar_tmp,2);
clear uchar_tmp;
uchar_tmp = uchar_00110000==snp_and_tag; sum1_and_00110000 = sum(uchar_tmp,1); sum2_and_00110000 = sum(uchar_tmp,2);
uchar_tmp = uchar_00110000==snp__or_tag; sum1__or_00110000 = sum(uchar_tmp,1); sum2__or_00110000 = sum(uchar_tmp,2);
uchar_tmp = uchar_00110000==snp_nnd_tag; sum1_nnd_00110000 = sum(uchar_tmp,1); sum2_nnd_00110000 = sum(uchar_tmp,2);
uchar_tmp = uchar_00110000==snp_mss_tag; sum1_mss_00110000 = sum(uchar_tmp,1); sum2_mss_00110000 = sum(uchar_tmp,2);
clear uchar_tmp;
uchar_tmp = uchar_11000000==snp_and_tag; sum1_and_11000000 = sum(uchar_tmp,1); sum2_and_11000000 = sum(uchar_tmp,2);
uchar_tmp = uchar_11000000==snp__or_tag; sum1__or_11000000 = sum(uchar_tmp,1); sum2__or_11000000 = sum(uchar_tmp,2);
uchar_tmp = uchar_11000000==snp_nnd_tag; sum1_nnd_11000000 = sum(uchar_tmp,1); sum2_nnd_11000000 = sum(uchar_tmp,2);
uchar_tmp = uchar_11000000==snp_mss_tag; sum1_mss_11000000 = sum(uchar_tmp,1); sum2_mss_11000000 = sum(uchar_tmp,2);
clear uchar_tmp;
clear uchar_00000011;clear uchar_00001100;clear uchar_00110000;clear uchar_11000000;
tmp_t = toc; disp(sprintf(' %% sums %0.2f',tmp_t)); tic;
n_and{ns} = sum2_and_00000011 + sum2_and_00001100 + sum2_and_00110000 + sum2_and_11000000;
n__or{ns} = sum2__or_00000011 + sum2__or_00001100 + sum2__or_00110000 + sum2__or_11000000;
n_nnd{ns} = sum2_nnd_00000011 + sum2_nnd_00001100 + sum2_nnd_00110000 + sum2_nnd_11000000;
n_mss{ns} = sum2_mss_00000011 + sum2_mss_00001100 + sum2_mss_00110000 + sum2_mss_11000000;
tmp_t = toc; disp(sprintf(' %% n_xxx %0.2f',tmp_t)); tic;
m_and{ns}(stride1) = m_and{ns}(stride1) + sum1_and_00000011;
m_and{ns}(stride2) = m_and{ns}(stride2) + sum1_and_00001100;
m_and{ns}(stride3) = m_and{ns}(stride3) + sum1_and_00110000;
m_and{ns}(stride4) = m_and{ns}(stride4) + sum1_and_11000000;
m__or{ns}(stride1) = m__or{ns}(stride1) + sum1__or_00000011;
m__or{ns}(stride2) = m__or{ns}(stride2) + sum1__or_00001100;
m__or{ns}(stride3) = m__or{ns}(stride3) + sum1__or_00110000;
m__or{ns}(stride4) = m__or{ns}(stride4) + sum1__or_11000000;
m_nnd{ns}(stride1) = m_nnd{ns}(stride1) + sum1_nnd_00000011;
m_nnd{ns}(stride2) = m_nnd{ns}(stride2) + sum1_nnd_00001100;
m_nnd{ns}(stride3) = m_nnd{ns}(stride3) + sum1_nnd_00110000;
m_nnd{ns}(stride4) = m_nnd{ns}(stride4) + sum1_nnd_11000000;
m_mss{ns}(stride1) = m_mss{ns}(stride1) + sum1_mss_00000011;
m_mss{ns}(stride2) = m_mss{ns}(stride2) + sum1_mss_00001100;
m_mss{ns}(stride3) = m_mss{ns}(stride3) + sum1_mss_00110000;
m_mss{ns}(stride4) = m_mss{ns}(stride4) + sum1_mss_11000000;
tmp_t = toc; disp(sprintf(' %% m_xxx %0.2f',tmp_t)); tic;
clear sum1_and_00000011;clear sum1_and_00001100;clear sum1_and_00110000;clear sum1_and_11000000;
clear sum1__or_00000011;clear sum1__or_00001100;clear sum1__or_00110000;clear sum1__or_11000000;
clear sum1_nnd_00000011;clear sum1_nnd_00001100;clear sum1_nnd_00110000;clear sum1_nnd_11000000;
clear sum1_mss_00000011;clear sum1_mss_00001100;clear sum1_mss_00110000;clear sum1_mss_11000000;
clear sum2_and_00000011;clear sum2_and_00001100;clear sum2_and_00110000;clear sum2_and_11000000;
clear sum2__or_00000011;clear sum2__or_00001100;clear sum2__or_00110000;clear sum2__or_11000000;
clear sum2_nnd_00000011;clear sum2_nnd_00001100;clear sum2_nnd_00110000;clear sum2_nnd_11000000;
clear sum2_mss_00000011;clear sum2_mss_00001100;clear sum2_mss_00110000;clear sum2_mss_11000000;
clear stride1 stride2 stride3 stride4;
end;%if ver_flag==0;
if ver_flag==1;
disp(sprintf(' %% version 1: sequential read'));
data_snp = zeros(1,4*npats_rup(ns),'uint8');
n_and{ns} = zeros(1,nsnps(ns));n__or{ns} = zeros(1,nsnps(ns));n_nnd{ns} = zeros(1,nsnps(ns));n_mss{ns} = zeros(1,nsnps(ns));
m_and{ns} = zeros(1,npats(ns));m__or{ns} = zeros(1,npats(ns));m_nnd{ns} = zeros(1,npats(ns));m_mss{ns} = zeros(1,npats(ns));
tmp_t_tot=0; tic;
for nl=1:nsnps(ns);
if (mod(nl,1000)==0); 
tmp_t = toc; tmp_t_tot = tmp_t_tot + tmp_t; tmp_t_full = tmp_t_tot*nsnps(ns)/nl; tmp_t_rem = tmp_t_full - tmp_t_tot;
disp(sprintf(' %% nl %d/%d, %0.2f elapsed (%0.2f estimated, %0.2f remaining)',nl,nsnps(ns),tmp_t,tmp_t_full,tmp_t_rem)); tic; 
end;%if (mod(nl,1000)==0); 
uchar_array = transpose(uint8(fread(fpbed{ns},npats_rup(ns),'uint8')));
uchar_00000011 = mod(uchar_array,  4);
uchar_00001111 = mod(uchar_array, 16); uchar_00001100 = idivide(uchar_00001111 - uchar_00000011, 4) ;
uchar_00111111 = mod(uchar_array, 64); uchar_00110000 = idivide(uchar_00111111 - uchar_00001111,16) ;
uchar_11111111 = mod(uchar_array,256); uchar_11000000 = idivide(uchar_11111111 - uchar_00111111,64) ;
data_snp = zeros(1,4*npats_rup(ns));
data_snp(1:4:4*npats_rup(ns)) = uchar_00000011;
data_snp(2:4:4*npats_rup(ns)) = uchar_00001100;
data_snp(3:4:4*npats_rup(ns)) = uchar_00110000;
data_snp(4:4:4*npats_rup(ns)) = uchar_11000000;
data_snp_tmp = data_snp(1:npats(ns));
n_and{ns}(nl) = sum(data_snp_tmp==snp_and_tag);
n__or{ns}(nl) = sum(data_snp_tmp==snp__or_tag);
n_nnd{ns}(nl) = sum(data_snp_tmp==snp_nnd_tag);
n_mss{ns}(nl) = sum(data_snp_tmp==snp_mss_tag);
m_and{ns} = m_and{ns} + (data_snp_tmp==snp_and_tag);
m__or{ns} = m__or{ns} + (data_snp_tmp==snp__or_tag);
m_nnd{ns} = m_nnd{ns} + (data_snp_tmp==snp_nnd_tag);
m_mss{ns} = m_mss{ns} + (data_snp_tmp==snp_mss_tag);
end;%for nl=1:nsnps(ns);
end;%if ver_flag==1;
fclose(fpbed{ns});
end;%for ns=1:nstd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Calculating entropy (i.e., KL-divergence) for each snp within study ;'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ns=1:nstd;
pp{ns} = n_and{ns}/npats(ns) + 0.5*n__or{ns}/npats(ns);
qq{ns} = n_nnd{ns}/npats(ns) + 0.5*n__or{ns}/npats(ns);
ii{ns} = n_nnd{ns}/npats(ns).*log((n_nnd{ns}/npats(ns))./(qq{ns}.^2)) ...
       + n__or{ns}/npats(ns).*log((n__or{ns}/npats(ns))./(2*pp{ns}.*qq{ns})) ...
       + n_and{ns}/npats(ns).*log((n_and{ns}/npats(ns))./(pp{ns}.^2));
disp(sprintf(' %% ns %d: [05 25 50 75 95] percentile: ',ns));
tmp = [prctile(pp{ns}, 5) , prctile(pp{ns},25) , prctile(pp{ns},50) , prctile(pp{ns},75) , prctile(pp{ns},95) ] ; 
disp(sprintf(' %% pp [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp));
tmp = [prctile(qq{ns}, 5) , prctile(qq{ns},25) , prctile(qq{ns},50) , prctile(qq{ns},75) , prctile(qq{ns},95) ] ; 
disp(sprintf(' %% qq [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp));
tmp = [prctile(ii{ns}, 5) , prctile(ii{ns},25) , prctile(ii{ns},50) , prctile(ii{ns},75) , prctile(ii{ns},95) ] ; 
disp(sprintf(' %% ii [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp));
end;%for ns=1:nstd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Calculating intersection of snps ;'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snp__id_ = cell(nstd,1); snp_pos_ = cell(nstd,1); for ns=1:nstd; snp__id_{ns} = cabim{ns}{2}; snp_pos_{ns} = cabim{ns}{4}; end;
[snp_cap] = intersectall(snp__id_); snp_cap_length = length(snp_cap);
snp_cap_ = cell(nstd,1);
for ns=1:nstd; [tmpa,tmpb,snp_cap_{ns}] = intersect(snp_cap,snp__id_{ns},'stable'); end;%for ns=1:nstd;
disp_flag=0;
if disp_flag;
figure;cla; 
for ns1=1:nstd;
for ns2=ns1+1:nstd;
subplot(1,2,1); hold on;
plot(1:snp_cap_length,snp_pos_{ns1}(snp_cap_{ns1}),'ro',1:snp_cap_length,snp_pos_{ns2}(snp_cap_{ns2}),'b.'); xlabel('snp index number'); ylabel('snp pos'); title('snp ij vs pos');
subplot(1,2,2); hold on;
plot(snp_pos_{ns1}(snp_cap_{ns1}),snp_pos_{ns2}(snp_cap_{ns2}),'.'); xlabel('snp pos'); ylabel('snp pos'); title('pos vs pos');
end;%for ns2=ns1+1:nstd;
end;%for ns1=1:nstd;
end;%if disp_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Calculating pp, qq and entropy for each snp across studies ;'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_and_tot = zeros(snp_cap_length,1); n__or_tot = zeros(snp_cap_length,1); n_nnd_tot = zeros(snp_cap_length,1); n_mss_tot = zeros(snp_cap_length,1); ii_max = zeros(snp_cap_length,1);
for ns=1:nstd; 
n_and_tot = n_and_tot + n_and{ns}(snp_cap_{ns}); n__or_tot = n__or_tot + n__or{ns}(snp_cap_{ns}); n_nnd_tot = n_nnd_tot + n_nnd{ns}(snp_cap_{ns}); n_mss_tot = n_mss_tot + n_mss{ns}(snp_cap_{ns}); 
ii_max = max(ii_max,ii{ns}(snp_cap_{ns}));
end;%for ns=1:nstd;
frq_and_tot = n_and_tot / sum(npats); frq__or_tot = n__or_tot / sum(npats); frq_nnd_tot = n_nnd_tot / sum(npats); frq_mss_tot = n_mss_tot / sum(npats);
fr2_and_tot = frq_and_tot; fr2__or_tot = frq__or_tot; fr2_nnd_tot = frq_nnd_tot; fr2_mode = zeros(snp_cap_length,1);
tmp_ij = find(frq_and_tot >frq__or_tot & frq_and_tot >frq_nnd_tot); fr2_mode(tmp_ij) = snp_and_tag; fr2_and_tot(tmp_ij) = fr2_and_tot(tmp_ij) + frq_mss_tot(tmp_ij);
tmp_ij = find(frq__or_tot>=frq_and_tot & frq__or_tot>=frq_nnd_tot); fr2_mode(tmp_ij) = snp__or_tag; fr2__or_tot(tmp_ij) = fr2__or_tot(tmp_ij) + frq_mss_tot(tmp_ij);
tmp_ij = find(frq_nnd_tot>=frq_and_tot & frq_nnd_tot >frq__or_tot); fr2_mode(tmp_ij) = snp_nnd_tag; fr2_nnd_tot(tmp_ij) = fr2_nnd_tot(tmp_ij) + frq_mss_tot(tmp_ij);
clear tmp_ij;
pp_tot = fr2_and_tot + 0.5*fr2__or_tot;
qq_tot = fr2_nnd_tot + 0.5*fr2__or_tot;
ii_tot = fr2_nnd_tot.*log(fr2_nnd_tot./(qq_tot.^2)) ...
  + fr2__or_tot.*log(fr2__or_tot./(2*pp_tot.*qq_tot)) ...
  + fr2_and_tot.*log(fr2_and_tot./(pp_tot.^2));
ii_tot(find(fr2_nnd_tot==0))=0;
disp(sprintf(' %% all studies: [05 25 50 75 95] percentile: ',ns));
tmp = [prctile(pp_tot, 5) , prctile(pp_tot,25) , prctile(pp_tot,50) , prctile(pp_tot,75) , prctile(pp_tot,95) ] ; 
disp(sprintf(' %% pp [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp));
tmp = [prctile(qq_tot, 5) , prctile(qq_tot,25) , prctile(qq_tot,50) , prctile(qq_tot,75) , prctile(qq_tot,95) ] ; 
disp(sprintf(' %% qq [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp));
tmp = [prctile(ii_tot, 5) , prctile(ii_tot,25) , prctile(ii_tot,50) , prctile(ii_tot,75) , prctile(ii_tot,95) ] ; 
disp(sprintf(' %% ii [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Calculating subset of col-intersection which satisfies cutoffs ; '));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snp_cap_sub = find(ii_tot < ent_cutoff & ii_max < ent_cutoff & min(pp_tot,qq_tot) > frq_cutoff & frq_mss_tot < mss_cutoff); snp_cap_sub_length = length(snp_cap_sub);
disp(sprintf(' %% retaining %d/%d snps',length(snp_cap_sub),length(snp_cap)));
snp_cap_sub_ = cell(nstd,1);
for ns=1:nstd; snp_cap_sub_{ns} = snp_cap_{ns}(snp_cap_sub); end;%for ns=1:nstd;
disp_flag=0;
if disp_flag;
figure;cla; 
for ns1=1:nstd;
for ns2=ns1+1:nstd;
subplot(1,2,1); hold on;
plot(1:snp_cap_sub_length,snp_pos_{ns1}(snp_cap_sub_{ns1}),'ro',1:snp_cap_sub_length,snp_pos_{ns2}(snp_cap_sub_{ns2}),'b.'); xlabel('snp index number'); ylabel('snp pos'); title('snp cij vs pos');
subplot(1,2,2); hold on;
plot(snp_pos_{ns1}(snp_cap_sub_{ns1}),snp_pos_{ns2}(snp_cap_sub_{ns2}),'.'); xlabel('snp pos'); ylabel('snp pos'); title('pos vs pos');
end;%for ns2=ns1+1:nstd;
end;%for ns1=1:nstd;
end;%if disp_flag;

brk=3; col_fr2 = nstd+1; col_srt = nstd+2; col_trs = nstd+3;
cij_index_array = zeros(2*snp_cap_sub_length,nstd+brk);
cij_and = transpose(1:snp_cap_sub_length);
cij__or = snp_cap_sub_length + transpose(1:snp_cap_sub_length);
for ns=1:nstd; cij_index_array(cij_and,ns) = snp_cap_sub_{ns}; cij_index_array(cij__or,ns) = snp_cap_sub_{ns}; end;% for ns=1:nstd;
cij_index_array(cij_and,col_fr2) = fr2_and_tot(snp_cap_sub);
cij_index_array(cij__or,col_fr2) = fr2_and_tot(snp_cap_sub) + fr2__or_tot(snp_cap_sub);
[tmp_v,tmp_r] = sort(cij_index_array(:,col_fr2),'ascend');
cij_index_array(:,col_srt) = tmp_r;
[tmp_v,tmp_s] = sort(tmp_r,'ascend');
cij_index_array(:,col_trs) = tmp_s;
cij_and_rpl = cij_index_array(cij_and,col_trs);
cij__or_rpl = cij_index_array(cij__or,col_trs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Calculating subset of rows which satisfies cutoffs ; '));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ns=1:nstd; 
rij_cap_{ns} = find(m_mss{ns}(1,1:npats(ns))/nsnps(ns) < mss_cutoff); 
disp(sprintf(' %% ns %d: retaining %d/%d patients',ns,length(rij_cap_{ns}),npats(ns)));
end;%for ns=1:nstd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Grabbing mds components ; '));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrows_T_full=0; ncols_Tn = 1 + size(mds,2);
for ns=1:nstd;
tmp_fam = cell(length(rij_cap_{ns}),1);
for np=1:length(rij_cap_{ns}); 
tmp_fam{np} = sprintf('%s%s%s',cafam{ns}{1}{rij_cap_{ns}(np)},'&',cafam{ns}{2}{rij_cap_{ns}(np)}); 
end;%for np=1:length(rij_cap_{ns});
[tmp_v,tmp_r,tmp_m] = intersect(tmp_fam,mds_names,'stable');
if length(tmp_r)<length(rij_cap_{ns}); 
disp(sprintf(' %% Warning! ns %d, tmp_r length %d/%d',ns,length(tmp_r),length(rij_cap_{ns}))); 
rij_cap_{ns} = rij_cap_{ns}(tmp_r);
 else; disp(sprintf(' %% ns %d, tmp_r %d == length(rij_cap_{ns}) %d',ns,length(tmp_r),length(rij_cap_{ns})));
end;%if length(tmp_r)<length(rij_cap_{ns}); 
T_{ns} = mds(tmp_m,:);
nrows_T_full = nrows_T_full + length(tmp_r);
end;%for ns=1:nstd;
T_full = zeros(nrows_T_full,size(mds,2));
nrows_T_sum=0;
for ns=1:nstd;
T_full(nrows_T_sum + (1:size(T_{ns},1)),1:size(T_{ns},2)) = T_{ns};
nrows_T_sum = nrows_T_sum + size(T_{ns},1);
end;%for ns=1:nstd;
T_med = median(T_full);
T_full = T_full - repmat(T_med,nrows_T_sum,1);
for ns=1:nstd;
T_{ns} = T_{ns} - repmat(T_med,length(rij_cap_{ns}),1);
end;%for ns=1:nstd;
mc_T = ones(ncols_Tn,1);
tmpchar = sprintf('%s/%smc_T.b16',dir_out,snp_pre); tutorial_binary_compress(bitj,mc_T(:)>0,tmpchar); 
tmpchar = sprintf('%s/%sT_full_n.b16',dir_out,snp_pre);tutorial_binary_compress(bitj,[ones(nrows_T_sum,1) , (T_full>0)],tmpchar);
tmpchar = sprintf('%s/%sT_full_t.b16',dir_out,snp_pre);tutorial_binary_compress(bitj,transpose([ones(nrows_T_sum,1) , (T_full>0)]),tmpchar);
for ns=1:nstd;
tmpchar = sprintf('%s/%sT_%.2d_n.b16',dir_out,snp_pre,ns);tutorial_binary_compress(bitj,[ones(length(rij_cap_{ns}),1) , (T_{ns}>0)],tmpchar);
tmpchar = sprintf('%s/%sT_%.2d_t.b16',dir_out,snp_pre,ns);tutorial_binary_compress(bitj,transpose([ones(length(rij_cap_{ns}),1) , (T_{ns}>0)]),tmpchar);
end;%for ns=1:nstd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Initializing An_ and At_ ;'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncols_An = 2*snp_cap_sub_length; ncols_An_extend = mod(bitj - mod(ncols_An,bitj),bitj); lcols_An = (ncols_An + ncols_An_extend)/bit8;
nrows_An_ = zeros(1,nstd);
for ns=1:nstd;
nrows_An_(ns) = length(rij_cap_{ns}); nrows_An_extend_(ns) = mod(bitj - mod(nrows_An_(ns),bitj),bitj); lrows_An_(ns) = (nrows_An_(ns) + nrows_An_extend_(ns))/bit8;
end;%for ns=1:nstd;
nrows_An_full = sum(nrows_An_); nrows_An_full_extend = mod(bitj - mod(nrows_An_full,bitj),bitj); lrows_An_full = (nrows_An_full + nrows_An_full_extend)/bit8;
nrows_An_csum = cumsum([0,nrows_An_]);
clear An_ At_ ;
An_full = zeros(lrows_An_full,ncols_An,'uint8');
At_full = zeros(lcols_An,nrows_An_full,'uint8');
for ns=1:nstd;
An_{ns} = zeros(lrows_An_(ns),ncols_An,'uint8');
At_{ns} = zeros(lcols_An,nrows_An_(ns),'uint8');
end;%for ns=1:nstd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Writing extended bim and fam files (i.e., bim.ext and fam.ext) ; '));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bim_flag=1;
if bim_flag;
tmp_bim = cell(6 + 4);
tmp_bim{1} = zeros(ncols_An,1);
tmp_bim{2} = cell(ncols_An,1);
tmp_bim{3} = zeros(ncols_An,1);
tmp_bim{4} = zeros(ncols_An,1);
tmp_bim{5} = zeros(ncols_An,1);
tmp_bim{6} = zeros(ncols_An,1);
tmp_bim{7} = cell(ncols_An,1);
tmp_bim{8} = zeros(ncols_An,1);
tmp_bim{9} = zeros(ncols_An,1);
tmp_bim{10} = zeros(ncols_An,1);
for nl3=1:snp_cap_sub_length; 
nl2 = snp_cap_sub(nl3);
ns=1; nla = snp_cap_sub_{ns}(nl3); nlb = snp_cap_{ns}(nl2); assert(nla==nlb); nl1=nla;
tmp_bim{1}(cij_and_rpl(nl3)) = cabim{ns}{1}(nl1);
tmp_bim{2}(cij_and_rpl(nl3)) = cabim{ns}{2}(nl1);
tmp_bim{3}(cij_and_rpl(nl3)) = cabim{ns}{3}(nl1);
tmp_bim{4}(cij_and_rpl(nl3)) = cabim{ns}{4}(nl1);
tmp_bim{5}(cij_and_rpl(nl3)) = cabim{ns}{5}(nl1);
tmp_bim{6}(cij_and_rpl(nl3)) = cabim{ns}{6}(nl1);
tmp_bim{7}(cij_and_rpl(nl3)) = {'and'};
tmp_bim{8}(cij_and_rpl(nl3)) = ii_tot(nl2);
tmp_bim{9}(cij_and_rpl(nl3)) = fr2_and_tot(nl2);
tmp_bim{10}(cij_and_rpl(nl3)) = frq_mss_tot(nl2);
tmp_bim{1}(cij__or_rpl(nl3)) = cabim{ns}{1}(nl1);
tmp_bim{2}(cij__or_rpl(nl3)) = cabim{ns}{2}(nl1);
tmp_bim{3}(cij__or_rpl(nl3)) = cabim{ns}{3}(nl1);
tmp_bim{4}(cij__or_rpl(nl3)) = cabim{ns}{4}(nl1);
tmp_bim{5}(cij__or_rpl(nl3)) = cabim{ns}{5}(nl1);
tmp_bim{6}(cij__or_rpl(nl3)) = cabim{ns}{6}(nl1);
tmp_bim{7}(cij__or_rpl(nl3)) = {'or'};
tmp_bim{8}(cij__or_rpl(nl3)) = ii_tot(nl2);
tmp_bim{9}(cij__or_rpl(nl3)) = fr2_and_tot(nl2) + fr2__or_tot(nl2);
tmp_bim{10}(cij__or_rpl(nl3)) = frq_mss_tot(nl2);
end;%for nl3=1:snp_cap_sub_length;
tmpchar = sprintf('%s/%sbim.ext',dir_out,snp_pre);
fid = fopen(tmpchar,'w');
for nc=1:ncols_An;
%fprintf(fid,'%d\t%s\t%d\t%d\t%c\t%c\t%s\t%f\t%f\t%f\n',tmp_bim{1}(nc),tmp_bim{2}(nc),tmp_bim{3}(nc),tmp_bim{4}(nc),tmp_bim{5}(nc),tmp_bim{6}(nc),tmp_bim{7}(nc),tmp_bim{8}(nc),tmp_bim{9}(nc),tmp_bim{10}(nc));
% strings must be input differently ;
% fields are: ; 
%{
  Chromosome
  Marker ID
  Genetic distance
  Physical position
  Allele 1
  Allele 2
  bit type (and vs or) 
  ii_tot (entropy)
  fr2_xxx_tot (sparsity of column)
  frq_mss_tot (missingness of snp)
    %}
fprintf(fid,'%d\t%s\t%d\t%d\t%c\t%c\t%s\t%f\t%f\t%f\n',tmp_bim{1}(nc),tmp_bim{2}{nc},tmp_bim{3}(nc),tmp_bim{4}(nc),tmp_bim{5}(nc),tmp_bim{6}(nc),tmp_bim{7}{nc},tmp_bim{8}(nc),tmp_bim{9}(nc),tmp_bim{10}(nc));
end;%for nc=1:ncols_An;
fclose(fid);
end;% if bim_flag;
fam_flag=1;
if fam_flag;
tmp_fam = cell(6 + 1);
tmp_fam{1} = cell(nrows_An_full,1);
tmp_fam{2} = cell(nrows_An_full,1);
tmp_fam{3} = cell(nrows_An_full,1);
tmp_fam{4} = cell(nrows_An_full,1);
tmp_fam{5} = zeros(nrows_An_full,1);
tmp_fam{6} = zeros(nrows_An_full,1);
tmp_fam{7} = cell(nrows_An_full,1);
nr=1;
for ns=1:nstd;
for np2=1:nrows_An_(ns);
np1 = rij_cap_{ns}(np2);
tmp_fam{1}(nr) = cafam{ns}{1}(np1);
tmp_fam{2}(nr) = cafam{ns}{2}(np1);
tmp_fam{3}(nr) = cafam{ns}{3}(np1);
tmp_fam{4}(nr) = cafam{ns}{4}(np1);
tmp_fam{5}(nr) = cafam{ns}{5}(np1);
tmp_fam{6}(nr) = cafam{ns}{6}(np1);
tmp_fam{7}(nr) = {sprintf('%s/%s',dir_base{ns},fpre{ns})};
%disp(sprintf('%s %s %s %s %d %d %s',cafam{ns}{1}{np1},cafam{ns}{2}{np1},cafam{ns}{3}{np1},cafam{ns}{4}{np1},cafam{ns}{5}(np1),cafam{ns}{6}(np1),sprintf('%s/%s',dir_base{ns},fpre{ns})));
nr = nr+1;
end;%for np2=1:nrows_An_(ns);
end;%for ns=1:nstd;
tmpchar = sprintf('%s/%sfam.ext',dir_out,snp_pre);
fid = fopen(tmpchar,'w');
for nr=1:nrows_An_full;
fprintf(fid,'%s\t%s\t%s\t%s\t%d\t%d\t%s\n',tmp_fam{1}{nr},tmp_fam{2}{nr},tmp_fam{3}{nr},tmp_fam{4}{nr},tmp_fam{5}(nr),tmp_fam{6}(nr),tmp_fam{7}{nr});
end;%for nr=1:nrows_An_full;
fclose(fid);
end;%if fam_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Writing mc and mr masks ; '));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mc_A = ones(ncols_An,1);
tmpchar = sprintf('%s/%smc_A.b16',dir_out,snp_pre); tutorial_binary_compress(bitj,mc_A(:)>0,tmpchar); 
mr_A_full = zeros(nrows_An_full,1); mr_Z_full = zeros(nrows_An_full,1); 
for ns=1:nstd;
mr_A_{ns} = zeros(nrows_An_(ns),1); mr_Z_{ns} = zeros(nrows_An_(ns),1);
tmp_ij = find(cafam{ns}{6}(rij_cap_{ns})==case_tag);
mr_A_{ns}(tmp_ij) = 1; mr_A_full(nrows_An_csum(ns) + tmp_ij) = 1;
tmp_ij = find(cafam{ns}{6}(rij_cap_{ns})==ctrl_tag);
mr_Z_{ns}(tmp_ij) = 1; mr_Z_full(nrows_An_csum(ns) + tmp_ij) = 1;
end;%for ns=1:nstd;
tmpchar = sprintf('%s/%smr_A_full.b16',dir_out,snp_pre);tutorial_binary_compress(bitj,mr_A_full(:)>0,tmpchar); pause(1);
tmpchar = sprintf('%s/%smr_Z_full.b16',dir_out,snp_pre);tutorial_binary_compress(bitj,mr_Z_full(:)>0,tmpchar); pause(1);
for ns=1:nstd;
tmpchar = sprintf('%s/%smr_A_%.2d.b16',dir_out,snp_pre,ns);tutorial_binary_compress(bitj,mr_A_{ns}(:)>0,tmpchar); pause(1);
tmpchar = sprintf('%s/%smr_Z_%.2d.b16',dir_out,snp_pre,ns);tutorial_binary_compress(bitj,mr_Z_{ns}(:)>0,tmpchar); pause(1);
end;%for ns=1:nstd;

%return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('Writing An_ and At_ ; '));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ns=1:nstd;
disp(sprintf(' %% ns %d: %s/%s',ns,dir_base{ns},fpre{ns}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reading bed file (binarized snp data) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpbed{ns} = fopen(fnbed{ns}); 
key1 = fread(fpbed{ns},1,'uint8'); key2 = fread(fpbed{ns},1,'uint8'); key3 = fread(fpbed{ns},1,'uint8'); assert(key1==108); assert(key2==27); assert(key3==1);
ver_flag=0;
if ver_flag==0;
disp(sprintf(' %% version 0: simultaneous read'));tic;
uchar_array = transpose(fread(fpbed{ns},[npats_rup(ns),nsnps(ns)],'uint8=>uint8')); fclose(fpbed{ns});
tmp_t = toc; disp(sprintf(' %% simultaneous read %0.2f',tmp_t)); tic;
disp(sprintf(' %% assuming more than 8G of ram... '));
uchar_00000011 = mod(uchar_array,  4);
uchar_00001111 = mod(uchar_array, 16); uchar_00001100 = idivide(uchar_00001111 - uchar_00000011, 4) ; 
uchar_00111111 = mod(uchar_array, 64); uchar_00110000 = idivide(uchar_00111111 - uchar_00001111,16) ; clear uchar_00001111;
                                       uchar_11000000 = idivide(uchar_array    - uchar_00111111,64) ; clear uchar_00111111;
clear uchar_array;
stride1 = 1:4:4*npats_rup(ns); tmp_ij = find(stride1>npats(ns)); uchar_00000011(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride2 = 2:4:4*npats_rup(ns); tmp_ij = find(stride2>npats(ns)); uchar_00001100(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride3 = 3:4:4*npats_rup(ns); tmp_ij = find(stride3>npats(ns)); uchar_00110000(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride4 = 4:4:4*npats_rup(ns); tmp_ij = find(stride4>npats(ns)); uchar_11000000(:,tmp_ij)=snp_not_tag; clear tmp_ij;
tmp_t = toc; disp(sprintf(' %% mod %0.2f',tmp_t)); tic;
clear uchar_tmp;
%disp(sprintf(' %% ns %d: filling missing data with mode',ns));
%for nl2=1:snp_cap_length;
%if (mod(nl2,1000)==0); disp(sprintf(' %% nl2 %d/%d',nl2,snp_cap_length)); end;
%nl1 = snp_cap_{ns}(nl2);
%tmp_ij = find(uchar_00000011(nl1,:)==snp_mss_tag); uchar_00000011(nl1,tmp_ij) = fr2_mode(nl2); clear tmp_ij;
%tmp_ij = find(uchar_00001100(nl1,:)==snp_mss_tag); uchar_00001100(nl1,tmp_ij) = fr2_mode(nl2); clear tmp_ij;
%tmp_ij = find(uchar_00110000(nl1,:)==snp_mss_tag); uchar_00110000(nl1,tmp_ij) = fr2_mode(nl2); clear tmp_ij;
%tmp_ij = find(uchar_11000000(nl1,:)==snp_mss_tag); uchar_11000000(nl1,tmp_ij) = fr2_mode(nl2); clear tmp_ij;
%end;%for nl2=1:snp_cap_length;
br_u = zeros(1,bit8,'uint8'); br_u = cast(2.^(bit8-1:-1:0),'uint8');
br_d = zeros(1,bit8,'double'); br_d = cast(2.^(bit8-1:-1:0),'double');
for nl3=1:snp_cap_sub_length; 
nl2 = snp_cap_sub(nl3);
nla = snp_cap_sub_{ns}(nl3); nlb = snp_cap_{ns}(nl2); assert(nla==nlb); nl1=nla;
if (mod(nl3,1000)==0); disp(sprintf(' %% nl3 %d/%d nl2 %d/%d nl1 %d/%d',nl3,snp_cap_sub_length,nl2,snp_cap_length,nl1,nsnps(ns))); end;
data_snp = zeros(1,4*npats_rup(ns));
data_snp(stride1) = uchar_00000011(nl1,:);
data_snp(stride2) = uchar_00001100(nl1,:);
data_snp(stride3) = uchar_00110000(nl1,:);
data_snp(stride4) = uchar_11000000(nl1,:);
data_snp = data_snp(1,1:npats(ns));
tmp_ij = find(data_snp==snp_mss_tag); data_snp(1,tmp_ij) = fr2_mode(nl2); clear tmp_ij;
data_snp_and = (data_snp==snp_and_tag) ; 
data_snp_and_sub = data_snp_and(1,rij_cap_{ns}); data_snp_and_bin_u = cast(data_snp_and_sub(1,:),'uint8');
data_snp_and_sub_full = zeros(1,nrows_An_full); data_snp_and_sub_full(1,nrows_An_csum(ns) + (1:nrows_An_(ns))) = data_snp_and_sub;
data_snp_and_bin_d = cast(data_snp_and_sub(1,:),'double'); 
data_snp_and_full_bin_d = cast(data_snp_and_sub_full(1,:),'double'); 
An_full(1:lrows_An_full,1+(cij_and_rpl(nl3)-1)) = An_full(1:lrows_An_full,1+(cij_and_rpl(nl3)-1)) + cast(transpose(br_d*reshape([data_snp_and_full_bin_d,zeros(1,nrows_An_full_extend)],bit8,lrows_An_full)),'uint8');
At_full(1+floor((cij_and_rpl(nl3)-1)/bit8),nrows_An_csum(ns) + (1:nrows_An_(ns))) = At_full(1+floor((cij_and_rpl(nl3)-1)/bit8),nrows_An_csum(ns) + (1:nrows_An_(ns))) + br_u(1+mod((cij_and_rpl(nl3)-1),bit8))*data_snp_and_bin_u;
An_{ns}(1:lrows_An_(ns),1+(cij_and_rpl(nl3)-1)) = An_{ns}(1:lrows_An_(ns),1+(cij_and_rpl(nl3)-1)) + cast(transpose(br_d*reshape([data_snp_and_bin_d,zeros(1,nrows_An_extend_(ns))],bit8,lrows_An_(ns))),'uint8');
At_{ns}(1+floor((cij_and_rpl(nl3)-1)/bit8),(1:nrows_An_(ns))) = At_{ns}(1+floor((cij_and_rpl(nl3)-1)/bit8),(1:nrows_An_(ns))) + br_u(1+mod((cij_and_rpl(nl3)-1),bit8))*data_snp_and_bin_u;
data_snp__or = (data_snp==snp_and_tag) | (data_snp==snp__or_tag) ;
data_snp__or_sub = data_snp__or(1,rij_cap_{ns}); data_snp__or_bin_u = cast(data_snp__or_sub(1,:),'uint8');
data_snp__or_sub_full = zeros(1,nrows_An_full); data_snp__or_sub_full(1,nrows_An_csum(ns) + (1:nrows_An_(ns))) = data_snp__or_sub;
data_snp__or_bin_d = cast(data_snp__or_sub(1,:),'double'); 
data_snp__or_full_bin_d = cast(data_snp__or_sub_full(1,:),'double'); 
An_full(1:lrows_An_full,1+(cij__or_rpl(nl3)-1)) = An_full(1:lrows_An_full,1+(cij__or_rpl(nl3)-1)) + cast(transpose(br_d*reshape([data_snp__or_full_bin_d,zeros(1,nrows_An_full_extend)],bit8,lrows_An_full)),'uint8');
At_full(1+floor((cij__or_rpl(nl3)-1)/bit8),nrows_An_csum(ns) + (1:nrows_An_(ns))) = At_full(1+floor((cij__or_rpl(nl3)-1)/bit8),nrows_An_csum(ns) + (1:nrows_An_(ns))) + br_u(1+mod((cij__or_rpl(nl3)-1),bit8))*data_snp__or_bin_u;
An_{ns}(1:lrows_An_(ns),1+(cij__or_rpl(nl3)-1)) = An_{ns}(1:lrows_An_(ns),1+(cij__or_rpl(nl3)-1)) + cast(transpose(br_d*reshape([data_snp__or_bin_d,zeros(1,nrows_An_extend_(ns))],bit8,lrows_An_(ns))),'uint8');
At_{ns}(1+floor((cij__or_rpl(nl3)-1)/bit8),(1:nrows_An_(ns))) = At_{ns}(1+floor((cij__or_rpl(nl3)-1)/bit8),(1:nrows_An_(ns))) + br_u(1+mod((cij__or_rpl(nl3)-1),bit8))*data_snp__or_bin_u;
end;%for nl3=1:snp_cap_sub_length;
clear uchar_00000011;clear uchar_00001100;clear uchar_00110000;clear uchar_11000000;
clear stride1 stride2 stride3 stride4;
end;%if ver_flag==0;
end;%for ns=1:nstd;

tmpAnchar = sprintf('%s/%sA_full_n.b16',dir_out,snp_pre);
fid = fopen(tmpAnchar,'w');
fwrite(fid,bitj,'int');
fwrite(fid,nrows_An_full,'int');
fwrite(fid,ncols_An,'int');
fwrite(fid,An_full(:),'uint8');
fclose(fid);
tmpAtchar = sprintf('%s/%sA_full_t.b16',dir_out,snp_pre);
fid = fopen(tmpAtchar,'w');
fwrite(fid,bitj,'int');
fwrite(fid,ncols_An,'int');
fwrite(fid,nrows_An_full,'int');
fwrite(fid,At_full(:),'uint8');
fclose(fid);
for ns=1:nstd;
tmpAnchar = sprintf('%s/%sA_%.2d_n.b16',dir_out,snp_pre,ns);
fid = fopen(tmpAnchar,'w');
fwrite(fid,bitj,'int');
fwrite(fid,nrows_An_(ns),'int');
fwrite(fid,ncols_An,'int');
fwrite(fid,An_{ns}(:),'uint8');
fclose(fid);
tmpAtchar = sprintf('%s/%sA_%.2d_t.b16',dir_out,snp_pre,ns);
fid = fopen(tmpAtchar,'w');
fwrite(fid,bitj,'int');
fwrite(fid,ncols_An,'int');
fwrite(fid,nrows_An_(ns),'int');
fwrite(fid,At_{ns}(:),'uint8');
fclose(fid);
end;%for ns=1:nstd;

check_flag=0;
if check_flag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('checking .b16 files for consistency ;'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpAnchar = sprintf('%s/%sA_full_n.b16',dir_out,snp_pre);
Bn_full = tutorial_binary_uncompress(tmpAnchar,1:nrows_An_full,1:ncols_An);
tmpAtchar = sprintf('%s/%sA_full_t.b16',dir_out,snp_pre);
Bt_full = tutorial_binary_uncompress(tmpAtchar,1:ncols_An,1:nrows_An_full);
disp(sprintf(' %% ns all: error %f',norm(Bn_full-transpose(Bt_full))));
subplot(1,2,1);imagesc(Bn_full,[-1,1]); subplot(1,2,2);imagesc(transpose(Bt_full),[-1,1]);
print('-djpeg',sprintf('%s/%sA_full_x.jpg',dir_out,snp_pre));
for ns=1:nstd;
tmpAnchar = sprintf('%s/%sA_%.2d_n.b16',dir_out,snp_pre,ns);
Bn_{ns} = tutorial_binary_uncompress(tmpAnchar,1:nrows_An_(ns),1:ncols_An);
tmpAtchar = sprintf('%s/%sA_%.2d_t.b16',dir_out,snp_pre,ns);
Bt_{ns} = tutorial_binary_uncompress(tmpAtchar,1:ncols_An,1:nrows_An_(ns));
disp(sprintf(' %% ns %d: error %f',ns,norm(Bn_{ns}-transpose(Bt_{ns}))));
subplot(1,2,1);imagesc(Bn_{ns},[-1,1]); subplot(1,2,2);imagesc(transpose(Bt_{ns}),[-1,1]);
print('-djpeg',sprintf('%s/%sA_%.2d_x.jpg',dir_out,snp_pre,ns));
end;%for ns=1:nstd;
end;%if check_flag;




