function test_loader_dexcluster_nonbinary_dr_0(E_array_,gamma,n_shuffle,prefix);
if nargin<2; gamma = 0.01; end;
if nargin<3; n_shuffle = 0; end;
if nargin<4; prefix = 'test_dexcluster'; end;

tmp_n_u = size(E_array_,1);
tmp_n_g = size(E_array_,2);
dir_code = '/data/rangan/dir_bcc/dir_lakcluster_c_dev'; dir_trunk = pwd;
rev_flag = 0; A_n_rind_ = {[1:tmp_n_u]}; A_n_cind = 1:tmp_n_g; Z_n_rind_ = {[]}; T_n_ = {ones(tmp_n_u,1)}; T_n_cind = 1;
GLOBAL_TEST_sparse = 0; B_MLT = 34; Ireq = 0;  verbose_flag = 0; flag_force_create = 0; flag_force_fig = 0;
%%%%%%%%;
shuffle_num=0;
test_string = sprintf('%s_%s',prefix,lakcluster_uADZSZDA_xfix_gen_ver1(rev_flag,A_n_rind_,Z_n_rind_,T_n_cind,GLOBAL_TEST_sparse,gamma,B_MLT,Ireq,shuffle_num));
disp(sprintf(' test_string: %s',test_string));
dir__in = sprintf('%s/dir_%s',dir_trunk,prefix);
disp(sprintf(' dir__in: %s',dir__in));
dir_out = sprintf('%s/dir_%s',dir__in,test_string); 
disp(sprintf(' dir_out: %s',dir_out));
dir_trace = sprintf('%s/dir_trace',dir_out); 
if (~exist(dir_trace,'dir')); disp(sprintf(' %% mkdir %s',dir_trace)); mkdir(dir_trace); end;
dir_jpg = sprintf('%s/dir_jpg',dir_out); 
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

str_trace = sprintf('%s/out_trace_dexnb_s%.4d.txt',dir_trace,shuffle_num);
str_xdrop = sprintf('%s/out_xdrop_a_dexnb.txt',dir_out);
if ( exist(str_trace,'file') & ~flag_force_create);
disp(sprintf(' %% %s found, not creating',str_trace));
out_trace_E_0_ = textread(str_trace);
end;%if ( exist(str_trace,'file') & ~flag_force_create);
if (~exist(str_trace,'file') |  flag_force_create);
disp(sprintf(' %% %s not found, creating',str_trace));
[out_xdrop_E_0_,out_trace_E_0_] = dexcluster_nonbinary_AAAA_ver0(E_array_,gamma);
fp=fopen(str_trace,'w');
fprintf(fp,'%d %d %d %0.16f %0.16f %d\n',transpose(out_trace_E_0_));
fclose(fp);
fp=fopen(str_xdrop,'w');
fprintf(fp,'%d %d\n',transpose(out_xdrop_E_0_));
fclose(fp);
end;%if (~exist(str_trace,'file') |  flag_force_create);
%%%%%%%%;
for nshuffle=1:n_shuffle;
shuffle_num = nshuffle;
test_string_s1 = sprintf('%s_%s',prefix,lakcluster_uADZSZDA_xfix_gen_ver1(rev_flag,A_n_rind_,Z_n_rind_,T_n_cind,GLOBAL_TEST_sparse,gamma,B_MLT,Ireq,shuffle_num));
dir_out_s1 = sprintf('%s/dir_%s',dir__in,test_string_s1); 
str_trace = sprintf('%s/out_trace_dexnb_s%.4d.txt',dir_trace,shuffle_num); 
if ( exist(str_trace,'file') & ~flag_force_create);
disp(sprintf(' %% %s found, not creating',str_trace));
end;%if ( exist(str_trace,'file') & ~flag_force_create);
if (~exist(str_trace,'file') |  flag_force_create);
disp(sprintf(' %% %s not found, creating',str_trace));
rng(nshuffle); [tmp_Q_,~] = qr(randn(tmp_n_u));
L_array_ = tmp_Q_*E_array_;
[out_xdrop_E_0_,out_trace_E_0_] = dexcluster_nonbinary_AAAA_ver0(L_array_,gamma);
clear tmp_Q_ L_array_;
fp=fopen(str_trace,'w');
fprintf(fp,'%d %d %d %0.16f %0.16f %d\n',transpose(out_trace_E_0_));
fclose(fp);
end;%if (~exist(str_trace,'file') |  flag_force_create);
end;%for nshuffle=1:n_shuffle;
%%%%%%%%;

%%%%%%%%;
str_xdrop = sprintf('%s/out_xdrop_a_dexnb.txt',dir_out);
out_xdrop_E_0_ = textread(str_xdrop);
rdrop_E_0_ = 1+out_xdrop_E_0_(find(out_xdrop_E_0_(:,1)>-1),1); rkeep_E_0_ = rdrop_E_0_(end:-1:1);
cdrop_E_0_ = 1+out_xdrop_E_0_(find(out_xdrop_E_0_(:,2)>-1),2); ckeep_E_0_ = cdrop_E_0_(end:-1:1);
%%%%%%%%;
trace_ = cell(1+n_shuffle,1);
for nshuffle=0:n_shuffle;
trace_{1+nshuffle} = textread(sprintf('%s/out_trace_dexnb_s%0.4d.txt',dir_trace,nshuffle));
end;%for nshuffle=1:n_shuffle;
QR_ = zeros(size(trace_{1},1),1+n_shuffle);
QC_ = zeros(size(trace_{1},1),1+n_shuffle);
for nshuffle=0:n_shuffle;
QR_(:,1+nshuffle) = trace_{1+nshuffle}(:,4);
QC_(:,1+nshuffle) = trace_{1+nshuffle}(:,5);
end;%for nshuffle=1:n_shuffle;
n_iter_ = trace_{1}(:,1);
r_rem_ = trace_{1}(:,2);
c_rem_ = trace_{1}(:,3);
%%%%%%%%;
QR_avg_ = mean(QR_(:,2:end),2); QR_std_ = std(QR_(:,2:end),[],2);
ZR_ = (QR_ - repmat(QR_avg_,1,1+n_shuffle))./repmat(QR_std_,1,1+n_shuffle);
QC_avg_ = mean(QC_(:,2:end),2); QC_std_ = std(QC_(:,2:end),[],2);
ZC_ = (QC_ - repmat(QC_avg_,1,1+n_shuffle))./repmat(QC_std_,1,1+n_shuffle);
%%%%%%%%;
[~,ni_cut] = max(ZR_(:,1));
r_rem_cut = r_rem_(ni_cut);
c_rem_cut = c_rem_(ni_cut);
%%%%%%%%;

%%%%%%%%;
flag_plot=1;
if flag_plot;
niteration_ = round(1:10:n_iter_(end));
for nl=1:length(niteration_);
niteration = niteration_(nl);
%%%%%%%%;
tmp_r_rem_cut = r_rem_(niteration);
tmp_c_rem_cut = c_rem_(niteration);
fname_base = sprintf('%s/%s_PCA_dexnb_ni%.4d_r%d_c%d_scatterplot',dir_jpg,prefix,niteration,tmp_r_rem_cut,tmp_c_rem_cut);
if flag_plot & (~exist(sprintf('%s.jpg',fname_base),'file') | flag_force_fig);
figure(1);
%%%%%%%%;
[tmp_U_,~,tmp_V_] = svds( E_array_(rkeep_E_0_(1:tmp_r_rem_cut),ckeep_E_0_(1:tmp_c_rem_cut)) ,3);
tmp_E1_pos_ = E_array_(:,ckeep_E_0_(1:tmp_c_rem_cut))*tmp_V_; %<-- original. ;
h1bin_ = linspace(min(tmp_E1_pos_(:,1)),max(tmp_E1_pos_(:,1)),128);
hA1_ = hist(tmp_E1_pos_(rkeep_E_0_(1:end),1),h1bin_);
hD1_ = hist(tmp_E1_pos_(rkeep_E_0_(1:tmp_r_rem_cut),1),h1bin_);
hX1_ = hist(tmp_E1_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),1),h1bin_);
[tmp_U_,~,tmp_V_] = svds( 1.0*(E_array_(rkeep_E_0_(1:tmp_r_rem_cut),ckeep_E_0_(1:tmp_c_rem_cut))>0) ,3);
tmp_E0_pos_ = (1.0*(E_array_(:,ckeep_E_0_(1:tmp_c_rem_cut))>0))*tmp_V_; %<-- binarized. ;
h0bin_ = linspace(min(tmp_E0_pos_(:,1)),max(tmp_E0_pos_(:,1)),128);
hA0_ = hist(tmp_E0_pos_(rkeep_E_0_(1:end),1),h0bin_);
hD0_ = hist(tmp_E0_pos_(rkeep_E_0_(1:tmp_r_rem_cut),1),h0bin_);
hX0_ = hist(tmp_E0_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),1),h0bin_);
%%%%%%%%;
Msize = 5;
clf;
subplot(2,2,1);hold on;
plot(tmp_E1_pos_(rkeep_E_0_(1:tmp_r_rem_cut),1),tmp_E1_pos_(rkeep_E_0_(1:tmp_r_rem_cut),2),'ro','MarkerSize',Msize);
%plot(tmp_E1_pos_(rkeep_E_0_(1:tmp_r_rem_cut),1),tmp_E1_pos_(rkeep_E_0_(1:tmp_r_rem_cut),2),'rx','MarkerSize',Msize);
plot(tmp_E1_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),1),tmp_E1_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),2),'bo','MarkerSize',Msize);
%plot(tmp_E1_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),1),tmp_E1_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),2),'bx','MarkerSize',Msize);
hold off;
xlabel('pc1'); ylabel('pc2'); title('E original');
subplot(2,2,2);hold on;
plot(tmp_E0_pos_(rkeep_E_0_(1:tmp_r_rem_cut),1),tmp_E0_pos_(rkeep_E_0_(1:tmp_r_rem_cut),2),'rx','MarkerSize',Msize);
%plot(tmp_E0_pos_(rkeep_E_0_(1:tmp_r_rem_cut),1),tmp_E0_pos_(rkeep_E_0_(1:tmp_r_rem_cut),2),'ro','MarkerSize',Msize);
plot(tmp_E0_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),1),tmp_E0_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),2),'bo','MarkerSize',Msize);
%plot(tmp_E0_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),1),tmp_E0_pos_(rkeep_E_0_(tmp_r_rem_cut+1:end),2),'bx','MarkerSize',Msize);
hold off;
xlabel('pc1'); ylabel('pc2'); title('E binarized');
subplot(2,2,3);hold on;
stairs(h1bin_,hA1_,'k-','LineWidth',1);
stairs(h1bin_,hD1_,'r-','LineWidth',2);
stairs(h1bin_,hX1_,'b-','LineWidth',1);
hold off;
xlabel('pc1'); ylabel('#'); title('E original histogram');
subplot(2,2,4);hold on;
stairs(h0bin_,hA0_,'k-','LineWidth',1);
stairs(h0bin_,hD0_,'r-','LineWidth',2);
stairs(h0bin_,hX0_,'b-','LineWidth',1);
hold off;
xlabel('pc1'); ylabel('#'); title('E binarized histogram');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot & (~exist(sprintf('%s.jpg',fname_base),'file') | flag_force_fig);
%%%%%%%%;
end;%for nl=1:length(niteration_);
end;%if flag_plot;

flag_plot=1;
fname_base = sprintf('%s/%s_QR_QC_ZR_ZC_dexnb',dir_jpg,prefix);
if flag_plot & (~exist(sprintf('%s.jpg',fname_base),'file') | flag_force_fig);
figure;clf;
subplot(2,4,1);hold on; plot(n_iter_,QR_(:,1),'r-','LineWidth',2); plot(n_iter_,QR_(:,2:end),'b-','LineWidth',0.5); hold off; xlabel('iteration'); ylabel('QR'); xlim([n_iter_(1),n_iter_(end)]); %ylim([-0.125,1.125]);
subplot(2,4,2);hold on; plot(r_rem_(1)-r_rem_,QR_(:,1),'r-','LineWidth',2); plot(r_rem_(1)-r_rem_,QR_(:,2:end),'b-','LineWidth',0.5); hold off; xlabel('rows removed'); ylabel('QR'); xlim([0,tmp_n_u]); %ylim([-0.125,1.125]);
subplot(2,4,3);hold on; plot(n_iter_,QC_(:,1),'r-','LineWidth',2); plot(n_iter_,QC_(:,2:end),'b-','LineWidth',0.5); hold off; xlabel('iteration'); ylabel('QC'); xlim([n_iter_(1),n_iter_(end)]); %ylim([-0.125,1.125]);
subplot(2,4,4);hold on; plot(c_rem_(1)-c_rem_,QC_(:,1),'r-','LineWidth',2); plot(c_rem_(1)-c_rem_,QC_(:,2:end),'b-','LineWidth',0.5); hold off; xlabel('cols removed'); ylabel('QC'); xlim([0,tmp_n_g]); %ylim([-0.125,1.125]);
subplot(2,4,5);hold on; plot(n_iter_,ZR_(:,1),'r-','LineWidth',2); plot(n_iter_,ZR_(:,2:end),'b-','LineWidth',0.5); hold off; xlabel('iteration'); ylabel('ZR'); xlim([n_iter_(1),n_iter_(end)]); %ylim([-0.125,1.125]);
subplot(2,4,6);hold on; plot(r_rem_(1)-r_rem_,ZR_(:,1),'r-','LineWidth',2); plot(r_rem_(1)-r_rem_,ZR_(:,2:end),'b-','LineWidth',0.5); hold off; xlabel('rows removed'); ylabel('ZR'); xlim([0,tmp_n_u]); %ylim([-0.125,1.125]);
subplot(2,4,7);hold on; plot(n_iter_,ZC_(:,1),'r-','LineWidth',2); plot(n_iter_,ZC_(:,2:end),'b-','LineWidth',0.5); hold off; xlabel('iteration'); ylabel('ZC'); xlim([n_iter_(1),n_iter_(end)]); %ylim([-0.125,1.125]);
subplot(2,4,8);hold on; plot(c_rem_(1)-c_rem_,ZC_(:,1),'r-','LineWidth',2); plot(c_rem_(1)-c_rem_,ZC_(:,2:end),'b-','LineWidth',0.5); hold off; xlabel('cols removed'); ylabel('ZC'); xlim([0,tmp_n_g]); %ylim([-0.125,1.125]);
figbig;
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot & (~exist(sprintf('%s.jpg',fname_base),'file') | flag_force_fig);

flag_plot=1;
fname_base = sprintf('%s/%s_A_n_sorted_dexnb',dir_jpg,prefix);
if flag_plot & (~exist(sprintf('%s.jpg',fname_base),'file') | flag_force_fig);
figure;clf;
subplot(1,2,1);colormap(colormap_beach()); imagesc(E_array_(rkeep_E_0_,ckeep_E_0_)); colorbar; title('dexcluster'); 
subplot(1,2,2);colormap(colormap_beach()); imagesc(E_array_(rkeep_E_0_,ckeep_E_0_)>0,[0,1]); colorbar; title('dexcluster'); 
figbig;
disp(sprintf(' %% writing %s.jpg',fname_base));
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
end;%if flag_plot & (~exist(sprintf('%s.jpg',fname_base),'file') | flag_force_fig);

