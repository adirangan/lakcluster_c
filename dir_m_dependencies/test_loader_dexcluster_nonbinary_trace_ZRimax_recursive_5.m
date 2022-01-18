function ...
[ ...
 output_label_ ...
,lpFmax_label_ ...
,lpnext_label_ ...
] = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_5( ...
 dir_trunk ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,p_prev ...
,flag_force_create ...
);

%%%%%%%%;
% Two p-value adjustments. ;
% First: instead of using p_set directly, we use p_use, which is defined so that: ;
% p_use + 2*p_use^2 + 4*p_use^3 + 8*p_uye^4 + ... = p_set. ;
% i.e., p_use = p_set/(1+2*p_set). ;
% Second: if we define q=1-p to be the probability of rejecting the null-hypothesis, ;
% then we only proceed if p_next := 1 - q_prev*q_branch <= p_use, ;
% where q_prev = cumulative rejection probability so far, ;
% and q_branch = conditional rejection for current branch. ;
%%%%%%%%;
% Fitting sample distribution with gumbel (p_gbranch). ;
%%%%%%%%;

na=0;
if (nargin<1+na); dir_trunk = []; end; na=na+1;
if (nargin<1+na); dir_out = []; end; na=na+1;
if (nargin<1+na); prefix_base = []; end; na=na+1;
if (nargin<1+na); E_array_base_ = []; end; na=na+1;
if (nargin<1+na); E_array_r_ij_ = []; end; na=na+1;
if (nargin<1+na); E_array_c_ij_ = []; end; na=na+1;
if (nargin<1+na); gamma = []; end; na=na+1;
if (nargin<1+na); n_shuffle = []; end; na=na+1;
if (nargin<1+na); p_set = []; end; na=na+1;
if (nargin<1+na); n_member_lob = []; end; na=na+1;
if (nargin<1+na); p_prev = []; end; na=na+1;
if (nargin<1+na); flag_force_create = []; end; na=na+1;

%%%%%%%%;
verbose=2;
if (isempty(dir_trunk)); dir_trunk = pwd(); end;
if (isempty(dir_out)); dir_out = []; end; %<-- will generate below. ;
if (isempty(prefix_base)); prefix_base = 'test'; end;
if (isempty(E_array_base_)); rng(0); E_array_base_ = randn(128,128); E_array_base_(1:24,1:24) = E_array_base_(1:24,1:24) + 2; end;
if (isempty(E_array_r_ij_)); E_array_r_ij_ = 1:size(E_array_base_,1); end;
if (isempty(E_array_c_ij_)); E_array_c_ij_ = 1:size(E_array_base_,2); end;
if (isempty(gamma)); gamma = 0.01; end;
if (isempty(n_shuffle)); n_shuffle = 64; end;
if (isempty(p_set)); p_set = 0.05; end;
if (isempty(n_member_lob)); n_member_lob = 3; end;
if (isempty(p_prev)); p_prev = 0.00; end;
if (isempty(flag_force_create)); flag_force_create = 0; end;
p_use = p_set/(1+2*p_set);
%%%%%%%%;
%n_u = size(E_array_base_(E_array_r_ij_,:),1); n_g = size(E_array_base_(E_array_r_ij_,E_array_c_ij_),2);
n_u = numel(E_array_r_ij_); n_g = numel(E_array_c_ij_);
rev_flag = 0; A_n_rind_ = {[1:n_u]}; A_n_cind = 1:n_g; Z_n_rind_ = {[]}; T_n_ = {ones(n_u,1)}; T_n_cind = 1;
GLOBAL_TEST_sparse = 0; B_MLT = 0; Ireq = 0;  verbose_flag = 0;
%%%%%%%%;
if (isempty(dir_out));
if (verbose>1); disp(sprintf(' %% dir_out not provided, constructing.')); end;
shuffle_num=0;
test_string = sprintf('%s_%s',prefix_base,xxxcluster_uADZSZDA_xfix_gen_ver1('dexnb_trace_ZRmax',rev_flag,A_n_rind_,Z_n_rind_,T_n_cind,GLOBAL_TEST_sparse,gamma,B_MLT,Ireq,shuffle_num));
if (verbose>1); disp(sprintf(' test_string: %s',test_string)); end;
dir_0in = sprintf('%s/dir_%s',dir_trunk,prefix_base);
if (verbose>1); disp(sprintf(' dir_0in: %s',dir_0in)); end;
if (~exist(dir_0in,'dir')); disp(sprintf(' %% mkdir %s',dir_0in)); mkdir(dir_0in); end;
dir_out = sprintf('%s/dir_%s',dir_0in,test_string); 
if (verbose>1); disp(sprintf(' dir_out: %s',dir_out)); end;
end;%if (isempty(dir_out));
if (verbose>0); disp(sprintf(' %% [entering test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_5]: dir_out %s',dir_out)); end;
if (~exist(dir_out,'dir')); disp(sprintf(' %% mkdir %s',dir_out)); mkdir(dir_out); end;
%%%%%%%%;
fname_base = sprintf('%s/base.mat',dir_out);
if (~exist(fname_base,'file') |  flag_force_create);
if (verbose>1); disp(sprintf(' %% %s not found, creating',fname_base)); end;
save(fname_base ...
,'E_array_r_ij_' ...
,'E_array_c_ij_' ...
,'gamma' ...
,'n_shuffle' ...
,'prefix_base' ...
,'dir_out' ...
,'p_set' ...
);
end;%if (~exist(fname_base,'file') |  flag_force_create);
%%%%%%%%;
dir_trace = sprintf('%s/dir_trace',dir_out); 
%if (~exist(dir_trace,'dir')); disp(sprintf(' %% mkdir %s',dir_trace)); mkdir(dir_trace); end;
dir_jpg = sprintf('%s/dir_jpg',dir_out); 
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if nargout>0; %<-- extract labels. ;
output_label_ = cell(n_u,1); for nu=1:n_u; output_label_{nu} = '0'; end;
lpFmax_label_ = cell(n_u,1); for nu=1:n_u; lpFmax_label_{nu} = ''; end;
lpnext_label_ = cell(n_u,1); for nu=1:n_u; lpnext_label_{nu} = ''; end;
fname_branch = sprintf('%s/branch.mat',dir_out);
load(fname_branch,'type_code');
%%%%%%%%;
if (strcmp(type_code,'stop'));
if (verbose>1); disp(sprintf(' %% %s --> stop',dir_out)); end;
end;%if (strcmp(type_code,'stop'));
%%%%%%%%;
if (strcmp(type_code,'split'));
load(fname_branch ...
,'ZR_' ...
,'branch_A_r_ij_' ...
,'branch_B_r_ij_' ...
,'lp_R_max' ...
,'lp_R_F_max' ...
);
p_branch = exp(-lp_R_F_max); %p_branch = min(exp(-[lp_R_max,lp_R_F_max]));
[~,~,~,tmp_p_opt,tmp_p_emp] = gumbel_fit(max(ZR_(:,2:end),[],1),max(ZR_(:,1)));
p_gbranch = tmp_p_emp; if (p_gbranch<1/n_shuffle); p_gbranch = tmp_p_opt; end;
p_next = 1 - (1-p_prev)*(1-p_gbranch); %<-- q_next = q_prev*q_branch. ;
lpnext = -log(p_next);
if (p_next> p_use);
if (verbose>1); disp(sprintf(' %% %s --> stop (p_branch %0.6f p_gbranch %0.6f --> p_next %0.6f > p_use %0.6f)',dir_out,p_branch,p_gbranch,p_next,p_use)); end;
end;%if (p_next>p_use);
if (p_next<=p_use);
if (verbose>1); disp(sprintf(' %% %s --> split (p_branch %0.6f p_gbranch %0.6f --> p_next %0.6f < p_use %0.6f)',dir_out,p_branch,p_gbranch,p_next,p_use)); end;
dir_A_out = sprintf('%s/dir_A',dir_out); if (~exist(dir_A_out,'dir')); disp(sprintf(' %% mkdir %s',dir_A_out)); mkdir(dir_A_out); end;
E_array_A_r_ij_ = E_array_r_ij_(branch_A_r_ij_); E_array_A_c_ij_ = E_array_c_ij_;
if (verbose>1); disp(sprintf(' %% %s <-- calling (n_sample = %d)',dir_A_out,length(E_array_A_r_ij_))); end;
[ ...
 tmp_output_label_ ...
,tmp_lpFmax_label_ ...
,tmp_lpnext_label_ ...
] = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_5( ...
 dir_trunk ...
,dir_A_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_A_r_ij_ ...
,E_array_A_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,p_next ...
,flag_force_create ...
);
for nu=1:length(branch_A_r_ij_);
tmp_ij = branch_A_r_ij_(nu);
output_label_{tmp_ij} = sprintf('%s%s',sprintf('A'),tmp_output_label_{nu});
lpFmax_label_{tmp_ij} = sprintf('%s%s',sprintf('%+0.16f ',lp_R_F_max),tmp_lpFmax_label_{nu});
lpnext_label_{tmp_ij} = sprintf('%s%s',sprintf('%+0.16f ',lpnext),tmp_lpnext_label_{nu});
end;%for nu=1:length(branch_A_r_ij_);
dir_B_out = sprintf('%s/dir_B',dir_out); if (~exist(dir_B_out,'dir')); disp(sprintf(' %% mkdir %s',dir_B_out)); mkdir(dir_B_out); end;
E_array_B_r_ij_ = E_array_r_ij_(branch_B_r_ij_); E_array_B_c_ij_ = E_array_c_ij_;
if (verbose>1); disp(sprintf(' %% %s <-- calling (n_sample = %d)',dir_B_out,length(E_array_B_r_ij_))); end;
[ ...
 tmp_output_label_ ...
,tmp_lpFmax_label_ ...
,tmp_lpnext_label_ ...
] = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_5( ...
 dir_trunk ...
,dir_B_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_B_r_ij_ ...
,E_array_B_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,p_next ...
,flag_force_create ...
);
for nu=1:length(branch_B_r_ij_);
tmp_ij = branch_B_r_ij_(nu);
output_label_{tmp_ij} = sprintf('%s%s',sprintf('B'),tmp_output_label_{nu});
lpFmax_label_{tmp_ij} = sprintf('%s%s',sprintf('%+0.16f ',lp_R_F_max),tmp_lpFmax_label_{nu});
lpnext_label_{tmp_ij} = sprintf('%s%s',sprintf('%+0.16f ',lpnext),tmp_lpnext_label_{nu});
end;%for nu=1:length(branch_B_r_ij_);
end;%if (p_next<p_use);
end;%if (strcmp(type_code,'split'));
%%%%%%%%;
end;%if nargout>0; %<-- extract labels. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if nargout==0; %<-- calculation. ;

fname_trace_E_ = sprintf('%s/out_trace_.mat',dir_out);
fname_xdrop_E = sprintf('%s/out_xdrop_a.txt',dir_out);
if ( exist(fname_trace_E_,'file') &  exist(fname_xdrop_E,'file') & ~flag_force_create);
disp(sprintf(' %% %s and %s found, not creating',fname_trace_E_,fname_xdrop_E));
end;%if ( exist(fname_trace_E_,'file') &  exist(fname_xdrop_E,'file') & ~flag_force_create);
if (~exist(fname_trace_E_,'file') | ~exist(fname_xdrop_E,'file') |  flag_force_create);
disp(sprintf(' %% %s and %s not found, creating',fname_trace_E_,fname_xdrop_E));
trace_E_ = cell(1+n_shuffle,1);
%%%%%%%%;
nshuffle=0;
tmp_fname_xdrop_E = sprintf('%s/out_xdrop_a.txt',dir_out);
tmp_fname_trace_E = sprintf('%s/out_trace_s%.4d.txt',dir_trace,nshuffle);
if ( exist(tmp_fname_trace_E,'file') &  exist(tmp_fname_xdrop_E,'file') & ~flag_force_create);
out_trace_E_0_ = textread(tmp_fname_trace_E);
out_xdrop_E_0_ = textread(tmp_fname_xdrop_E);
end;%if ( exist(tmp_fname_trace_E,'file') &  exist(tmp_fname_xdrop_E,'file') & ~flag_force_create);
if (~exist(tmp_fname_trace_E,'file') | ~exist(tmp_fname_xdrop_E,'file') |  flag_force_create);
L_array_ = mean_center_0(E_array_base_(E_array_r_ij_,E_array_c_ij_));
[out_xdrop_E_0_,out_trace_E_0_] = dexcluster_nonbinary_AAAA_ver0(L_array_,gamma);
fp=fopen(tmp_fname_xdrop_E,'w'); fprintf(fp,'%d %d\n',transpose(out_xdrop_E_0_)); fclose(fp);
clear L_array_;
end;%if (~exist(tmp_fname_trace_E,'file') | ~exist(tmp_fname_xdrop_E,'file') |  flag_force_create);
trace_E_{1+nshuffle} = out_trace_E_0_;
clear out_xdrop_E_0_ out_trace_E_0_ ;
%%%%%%%%;
for nshuffle=1:n_shuffle;
tmp_fname_trace_E = sprintf('%s/out_trace_s%.4d.txt',dir_trace,nshuffle);
if ( exist(tmp_fname_trace_E,'file') & ~flag_force_create);
out_trace_E_0_ = textread(tmp_fname_trace_E);
end;%if ( exist(tmp_fname_trace_E,'file') & ~flag_force_create);
if (~exist(tmp_fname_trace_E,'file') |  flag_force_create);
if (verbose>2); disp(sprintf(' %% %d/%d: not found, creating',nshuffle,n_shuffle)); end;
rng(nshuffle); [tmp_Q_,~] = qr(randn(n_u));
L_array_ = mean_center_0(tmp_Q_*mean_center_0(E_array_base_(E_array_r_ij_,E_array_c_ij_)));
[out_xdrop_E_0_,out_trace_E_0_] = dexcluster_nonbinary_AAAA_ver0(L_array_,gamma);
clear tmp_Q_ L_array_;
end;%if (~exist(tmp_fname_trace_E,'file') |  flag_force_create);
trace_E_{1+nshuffle} = out_trace_E_0_;
clear out_xdrop_E_0_ out_trace_E_0_ ;
end;%for nshuffle=1:n_shuffle;
%%%%%%%%;
save(fname_trace_E_ ...
,'n_shuffle','trace_E_' ...
);
for nshuffle=0:n_shuffle;
tmp_fname_trace_E = sprintf('%s/out_trace_s%.4d.txt',dir_trace,nshuffle);
if ( exist(tmp_fname_trace_E,'file')); disp(sprintf(' %% found %s, deleting',tmp_fname_trace_E)); delete(tmp_fname_trace_E); end;
end;%for nshuffle=0:n_shuffle;
tmp_dir_trace = dir_trace; try; rmdir(tmp_dir_trace); catch; end;%try;
%%%%%%%%;
end;%if (~exist(fname_trace_E_,'file') | ~exist(fname_xdrop_E,'file') |  flag_force_create);
  
%%%%%%%%;
fname_branch = sprintf('%s/branch.mat',dir_out);
if ( exist(fname_branch,'file') & ~flag_force_create);
load(fname_branch);
end;%if ( exist(fname_branch,'file') & ~flag_force_create);
if (~exist(fname_branch,'file') |  flag_force_create);
if (verbose>1); disp(sprintf(' %% %s not found, creating',fname_branch)); end;
%%%%%%%%;
tmp_fname_xdrop_E = sprintf('%s/out_xdrop_a.txt',dir_out);
fcheck(tmp_fname_xdrop_E);
out_xdrop_E_0_ = textread(tmp_fname_xdrop_E);
rdrop_E_0_ = 1+out_xdrop_E_0_(find(out_xdrop_E_0_(:,1)>-1),1); rkeep_E_0_ = rdrop_E_0_(end:-1:1);
cdrop_E_0_ = 1+out_xdrop_E_0_(find(out_xdrop_E_0_(:,2)>-1),2); ckeep_E_0_ = cdrop_E_0_(end:-1:1);
%%%%%%%%;
load(fname_trace_E_);
trace_s0000_ = trace_E_{1};
QR_ = zeros(size(trace_s0000_,1),1+n_shuffle);
QC_ = zeros(size(trace_s0000_,1),1+n_shuffle);
for nshuffle=0:n_shuffle;
tmp_d = min(length(trace_s0000_(:,1)),length(trace_E_{1+nshuffle}(:,1)));
QR_(1:tmp_d,1+nshuffle) = trace_E_{1+nshuffle}(1:tmp_d,4);
QC_(1:tmp_d,1+nshuffle) = trace_E_{1+nshuffle}(1:tmp_d,5);
end;%for nshuffle=1:n_shuffle;
n_iter_ = trace_s0000_(:,1); n_iteration = n_iter_(end);
r_rem_ = trace_s0000_(:,2);
c_rem_ = trace_s0000_(:,3);
%%%%%%%%;
QR_avg_ = mean(QR_(:,2:end),2); QR_std_ = std(QR_(:,2:end),[],2);
ZR_ = (QR_ - repmat(QR_avg_,1,1+n_shuffle))./repmat(QR_std_,1,1+n_shuffle);
QC_avg_ = mean(QC_(:,2:end),2); QC_std_ = std(QC_(:,2:end),[],2);
ZC_ = (QC_ - repmat(QC_avg_,1,1+n_shuffle))./repmat(QC_std_,1,1+n_shuffle);
%%%%%%%%;
lp_R_ = z_to_p_0(ZR_); lp_C_ = z_to_p_0(ZC_);
[lp_R_max,ZR_ij] = Z_imax(0,-lp_R_(:,1),0);
r_rem_E_cut = r_rem_(ZR_ij); c_rem_E_cut = c_rem_(ZR_ij); n_rank = min(c_rem_E_cut,2);
L_array_ = mean_center_0(E_array_base_(E_array_r_ij_,E_array_c_ij_));
[tmp_U_,~,tmp_V_] = svds( L_array_(rkeep_E_0_(1:r_rem_E_cut),ckeep_E_0_(1:c_rem_E_cut)) ,n_rank);
tmp_E1_pos_ = L_array_(:,ckeep_E_0_(1:c_rem_E_cut))*tmp_V_; %<-- original. ;
E1_ZRmax_label_ = zeros(n_u,1); E1_ZRmax_label_(rkeep_E_0_(1:r_rem_E_cut)) = 1;
clear L_array_;
%%%%%%%%;
if n_iteration>1;
fname_fig = sprintf('%s/base_ZRmax',dir_jpg);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',sprintf('%s.jpg',fname_fig))); end;
figure(1);
clf;
subplot(1,2,1); hold on;
plot(n_iter_,-lp_R_(:,2:end),'k-','LineWidth',1);
plot(n_iter_,-lp_R_(:,1),'r-','LineWidth',2);
plot(n_iter_(ZR_ij),-lp_R_(ZR_ij,1),'go');
hold off;
xlim([0,n_iteration]);
xlabel('iteration'); ylabel('-log(p(ZR))'); title(sprintf('%s: -log(p(ZR)) vs iteration',dir_out),'Interpreter','none');
subplot(1,2,2); hold on;
scatter(tmp_E1_pos_(:,1),tmp_E1_pos_(:,n_rank),25,E1_ZRmax_label_); colormap('lines');
hold off;
xlabel('pc1');ylabel('pc2'); title('scatterplot');
figbig;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if n_iteration>1;
%%%%%%%%;
if n_iteration>1;
fname_fig = sprintf('%s/base_QZ',dir_jpg);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',sprintf('%s.jpg',fname_fig))); end;
figure(1);
clf;
subplot(1,2,1); hold on;
plot(n_iter_,ZR_(:,2:end),'k-','LineWidth',1);
plot(n_iter_,ZR_(:,1),'r-','LineWidth',2);
hold off;
xlim([0,n_iteration]);
xlabel('iteration'); ylabel('ZR'); title(sprintf('%s: ZR vs iteration',dir_out),'Interpreter','none');
subplot(1,2,2); hold on;
plot(n_u-r_rem_,ZR_(:,2:end),'k-','LineWidth',1);
plot(n_u-r_rem_,ZR_(:,1),'r-','LineWidth',2);
hold off;
xlim([0,n_u]);
xlabel('samples remaining'); ylabel('ZR'); title(sprintf('%s: ZR vs sample',dir_out),'Interpreter','none');
figbig;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if n_iteration>1;
%%%%%%%%;
flag_split = 0;
p_branch = 1;
%%%%%%%%;
if (lp_R_max< -log(p_use)); 
if (verbose>1); disp(sprintf(' %% lp_R_max %0.6f; not significant, stopping ;',lp_R_max));  end;
end;%if (lp_R_max< -log(p_use));  
%%%%%%%%;
if (lp_R_max>=-log(p_use)); 
if (verbose>1); disp(sprintf(' %% lp_R_max %0.6f; statistically significant, ZR_ij %d. ;',lp_R_max,ZR_ij)); end;
[lp_R_max,ZR_ij] = Z_imax(0,-lp_R_(:,1),-log(p_use));
if (verbose>1); disp(sprintf(' %% redefining: lp_R_max %0.6f; statistically significant, ZR_ij %d. ;',lp_R_max,ZR_ij)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now recalculate traces, dropping only rows. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
gamma_rdrop = 0;
fname_trace_F_ = sprintf('%s/out_trace_F_.mat',dir_out);
fname_xdrop_F = sprintf('%s/out_xdrop_F_a.txt',dir_out);
if ( exist(fname_trace_F_,'file') &  exist(fname_xdrop_F,'file') & ~flag_force_create);
disp(sprintf(' %% %s and %s found, not creating',fname_trace_F_,fname_xdrop_F));
end;%if ( exist(fname_trace_F_,'file') &  exist(fname_xdrop_F,'file') & ~flag_force_create);
if (~exist(fname_trace_F_,'file') | ~exist(fname_xdrop_F,'file') |  flag_force_create);
disp(sprintf(' %% %s and %s not found, creating',fname_trace_F_,fname_xdrop_F));
trace_F_ = cell(1+n_shuffle,1);
%%%%%%%%;
nshuffle=0;
tmp_fname_xdrop_F = sprintf('%s/out_xdrop_F_a.txt',dir_out);
tmp_fname_trace_F = sprintf('%s/out_trace_F_s%.4d.txt',dir_trace,nshuffle);
if ( exist(tmp_fname_trace_F,'file') &  exist(tmp_fname_xdrop_F,'file'));
out_trace_F_0_ = textread(tmp_fname_trace_F);
out_xdrop_F_0_ = textread(tmp_fname_xdrop_F);
end;%if ( exist(tmp_fname_trace_F,'file') &  exist(tmp_fname_xdrop_F,'file'));
if (~exist(tmp_fname_trace_F,'file') | ~exist(tmp_fname_xdrop_F,'file'));
L_array_ = mean_center_0(E_array_base_(E_array_r_ij_,E_array_c_ij_(ckeep_E_0_(1:c_rem_E_cut))));
[out_xdrop_F_0_,out_trace_F_0_] = dexcluster_nonbinary_AAAA_rdrop_ver0(L_array_,gamma_rdrop);
fp=fopen(tmp_fname_xdrop_F,'w'); fprintf(fp,'%d %d\n',transpose(out_xdrop_F_0_)); fclose(fp);
clear L_array_;
end;%if (~exist(tmp_fname_trace_F,'file') | ~exist(tmp_fname_xdrop_F,'file'));
trace_F_{1+nshuffle} = out_trace_F_0_;
clear out_xdrop_F_0_ out_trace_F_0_ ;
%%%%%%%%;
for nshuffle=1:n_shuffle;
tmp_fname_trace_F = sprintf('%s/out_trace_F_s%.4d.txt',dir_trace,nshuffle);
if ( exist(tmp_fname_trace_F,'file'));
out_trace_F_0_ = textread(tmp_fname_trace_F);
end;%if ( exist(tmp_fname_trace_F,'file'));
if (~exist(tmp_fname_trace_F,'file'));
if (verbose>2); disp(sprintf(' %% %d/%d: not found, creating',nshuffle,n_shuffle)); end;
rng(nshuffle); [tmp_Q_,~] = qr(randn(n_u));
L_array_ = mean_center_0(tmp_Q_*mean_center_0(E_array_base_(E_array_r_ij_,E_array_c_ij_(ckeep_E_0_(1:c_rem_E_cut)))));
[out_xdrop_F_0_,out_trace_F_0_] = dexcluster_nonbinary_AAAA_rdrop_ver0(L_array_,gamma_rdrop);
clear tmp_Q_ L_array_;
end;%if (~exist(tmp_fname_trace_F,'file'));
trace_F_{1+nshuffle} = out_trace_F_0_;
clear out_xdrop_F_0_ out_trace_F_0_ ;
end;%for nshuffle=1:n_shuffle;
%%%%%%%%;
save(fname_trace_F_ ...
,'n_shuffle','trace_F_' ...
);
for nshuffle=0:n_shuffle;
tmp_fname_trace_F = sprintf('%s/out_trace_F_s%.4d.txt',dir_trace,nshuffle);
if ( exist(tmp_fname_trace_F,'file')); disp(sprintf(' %% found %s, deleting',tmp_fname_trace_F)); delete(tmp_fname_trace_F); end;
end;%for nshuffle=0:n_shuffle;
tmp_dir_trace = dir_trace; try; rmdir(tmp_dir_trace); catch; end;%try;
%%%%%%%%;
end;%if (~exist(fname_trace_F_,'file') | ~exist(fname_xdrop_F,'file') |  flag_force_create);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% finished recalculating traces. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now extract out_xdrop_F. ;
%%%%%%%%;
tmp_fname_xdrop_F = sprintf('%s/out_xdrop_F_a.txt',dir_out);
fcheck(tmp_fname_xdrop_F); out_xdrop_F_0_ = textread(tmp_fname_xdrop_F);
rdrop_F_0_ = 1+out_xdrop_F_0_(find(out_xdrop_F_0_(:,1)>-1),1); rkeep_F_0_ = rdrop_F_0_(end:-1:1);
%%%%%%%%;
% Now extract trace_F. ;
%%%%%%%%;
load(fname_trace_F_);
trace_F_s0000_ = trace_F_{1};
QR_F_ = zeros(size(trace_F_s0000_,1),1+n_shuffle);
QC_F_ = zeros(size(trace_F_s0000_,1),1+n_shuffle);
for nshuffle=0:n_shuffle;
tmp_d = min(length(trace_F_s0000_(:,1)),length(trace_F_{1+nshuffle}(:,1)));
QR_F_(1:tmp_d,1+nshuffle) = trace_F_{1+nshuffle}(1:tmp_d,4);
QC_F_(1:tmp_d,1+nshuffle) = trace_F_{1+nshuffle}(1:tmp_d,5);
end;%for nshuffle=1:n_shuffle;
n_iter_F_ = trace_F_s0000_(:,1); n_iteration_F = n_iter_F_(end);
r_rem_F_ = trace_F_s0000_(:,2);
%%%%%%%%;
% Now calculate ZR_F. ;
%%%%%%%%;
QR_F_avg_ = mean(QR_F_(:,2:end),2); QR_F_std_ = std(QR_F_(:,2:end),[],2);
ZR_F_ = (QR_F_ - repmat(QR_F_avg_,1,1+n_shuffle))./repmat(QR_F_std_,1,1+n_shuffle);
QC_F_avg_ = mean(QC_F_(:,2:end),2); QC_F_std_ = std(QC_F_(:,2:end),[],2);
ZC_F_ = (QC_F_ - repmat(QC_F_avg_,1,1+n_shuffle))./repmat(QC_F_std_,1,1+n_shuffle);
%%%%%%%%;
% Now define F1_ZRmax_label_. ;
%%%%%%%%;
lp_R_F_ = z_to_p_0(ZR_F_); lp_R_F_ij_ = find( (r_rem_F_>n_member_lob) & ((n_u - r_rem_F_)>n_member_lob) );
[lp_R_F_max,ZR_F_ij] = Z_imax(0,-lp_R_F_(lp_R_F_ij_,1),-log(p_use)*1); ZR_F_ij = lp_R_F_ij_(ZR_F_ij); %<-- threshold at -log(p_use);
if (lp_R_F_max<=-log(p_use)); if (verbose>0); disp(sprintf(' %% Warning: lp_R_F_max %0.2f <=-log(p_use) %0.2f',lp_R_F_max,-log(p_use))); end; end;
if (lp_R_F_max> -log(p_use)); 
[lp_R_F_max,ZR_F_ij] = Z_imax(0,-lp_R_F_(lp_R_F_ij_,1),-log(p_use)); ZR_F_ij = lp_R_F_ij_(ZR_F_ij);
if (verbose>1); disp(sprintf(' %% redefining: lp_R_F_max %0.6f; statistically significant, ZR_F_ij %d',lp_R_F_max,ZR_F_ij)); end;
end;%if (lp_R_F_max> -log(p_use)); 
r_rem_F_cut = r_rem_F_(ZR_F_ij); n_rank = min(n_g,2);
L_array_ = mean_center_0(E_array_base_(E_array_r_ij_,E_array_c_ij_));
[tmp_U_,~,tmp_V_] = svds( L_array_(rkeep_F_0_(1:r_rem_F_cut),ckeep_E_0_(1:c_rem_E_cut)) ,n_rank);
tmp_F1_pos_ = L_array_(:,ckeep_E_0_(1:c_rem_E_cut))*tmp_V_; %<-- original. ;
F1_ZRmax_label_ = zeros(n_u,1); F1_ZRmax_label_(rkeep_F_0_(1:r_rem_F_cut)) = 1;
if (verbose>1); disp(sprintf(' %% cutting at iteration %d, [rem:%d,cut:%d] --> numel(F1_ZRmax_label_) %d',ZR_F_ij,r_rem_F_cut,n_u - r_rem_F_cut,numel(F1_ZRmax_label_))); end;
clear L_array_;
%%%%%%%%;
if ~isempty(ZR_F_ij) & ~isempty(ZR_ij) & n_iteration>1;
fname_fig = sprintf('%s/base_ZRmax_F',dir_jpg);
if (~exist(sprintf('%s.jpg',fname_fig),'file') |  flag_force_create);
if (verbose); disp(sprintf(' %% %s not found, creating',sprintf('%s.jpg',fname_fig))); end;
figure(1);
clf;
subplot(2,2,1); hold on;
plot(n_iter_,-lp_R_(:,2:end),'k-','LineWidth',1);
plot(n_iter_,-lp_R_(:,1),'r-','LineWidth',2);
plot(n_iter_(ZR_ij),-lp_R_(ZR_ij,1),'go');
hold off;
xlim([0,n_iteration]);
xlabel('iteration'); ylabel('-log(p(ZR))'); title(sprintf('%s: -log(p(ZR_E)) vs iteration',dir_out),'Interpreter','none');
subplot(2,2,3); hold on;
plot(n_iter_F_,-lp_R_F_(:,2:end),'k-','LineWidth',1);
plot(n_iter_F_,-lp_R_F_(:,1),'r-','LineWidth',2);
plot(n_iter_F_(ZR_F_ij),-lp_R_F_(ZR_F_ij,1),'go');
hold off;
xlim([0,n_iteration]);
xlabel('iteration'); ylabel('-log(p(ZR))'); title(sprintf('%s: -log(p(ZR_F)) vs iteration',dir_out),'Interpreter','none');
subplot(2,2,2); hold on;
scatter(tmp_E1_pos_(:,1),tmp_E1_pos_(:,n_rank),25,E1_ZRmax_label_); colormap('lines');
hold off;
xlabel('pc1');ylabel('pc2'); title('scatterplot E');
subplot(2,2,4); hold on;
scatter(tmp_F1_pos_(:,1),tmp_F1_pos_(:,n_rank),25,F1_ZRmax_label_); colormap('lines');
hold off;
xlabel('pc1');ylabel('pc2'); title('scatterplot F');
figbig;
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file') |  flag_force_create);
end;%if if ~isempty(ZR_F_ij) & ~isempty(ZR_ij) & n_iteration>1;
%%%%%%%%;
p_brach = 1;
if (~isempty(lp_R_F_max)); p_branch = exp(-lp_R_F_max); end; %<-- p_branch = min(exp(-[lp_R_max,lp_R_F_max]));
p_next = 1 - (1-p_prev)*(1-p_branch); %<-- q_next = q_prev*q_branch. ;
tmp_u_ = unique(F1_ZRmax_label_);
flag_split = (length(tmp_u_)==2);
flag_split = flag_split & (p_next<=p_use); %<-- only split if significant. ;
if (flag_split); flag_split = flag_split & length(find(F1_ZRmax_label_==tmp_u_(1)))>n_member_lob; end;
if (flag_split); flag_split = flag_split & length(find(F1_ZRmax_label_==tmp_u_(2)))>n_member_lob; end;
%%%%%%%%;
if ~flag_split;
type_code = 'stop';
branch_A_r_ij_ = [];
branch_B_r_ij_ = [];
else; 
type_code = 'split';
branch_A_r_ij_ = find(F1_ZRmax_label_==tmp_u_(1));
branch_B_r_ij_ = find(F1_ZRmax_label_==tmp_u_(2));
end;%if flag_split;
%%%%%%%%;
save(fname_branch ...
,'trace_E_' ...
,'trace_F_' ...
,'QR_' ...
,'QC_' ...
,'ZR_' ...
,'ZC_' ...
,'QR_F_' ...
,'ZR_F_' ...
,'lp_R_' ...
,'lp_R_max' ...
,'lp_R_F_' ...
,'lp_R_F_max' ...
,'p_set' ...
,'p_branch' ...
,'type_code' ...
,'ZR_ij' ...
,'ZR_F_ij' ...
,'branch_A_r_ij_' ...
,'branch_B_r_ij_' ...
,'E1_ZRmax_label_' ...
,'F1_ZRmax_label_' ...
);
%%%%%%%%;
end;%if (lp_R_max>=-log(p_use)); 
%%%%%%%%;
if ((lp_R_max< -log(p_use)) | ~flag_split);
type_code = 'stop';
branch_A_r_ij_ = [];
branch_B_r_ij_ = [];
if (~exist(fname_branch,'file'))
save(fname_branch ...
,'trace_E_' ...
,'QR_' ...
,'QC_' ...
,'ZR_' ...
,'ZC_' ...
,'lp_R_' ...
,'lp_R_max' ...
,'p_set' ...
,'p_branch' ...
,'type_code' ...
,'ZR_ij' ...
,'branch_A_r_ij_' ...
,'branch_B_r_ij_' ...
,'E1_ZRmax_label_' ...
);
end;%if (~exist(fname_branch,'file') |  flag_force_create); 
if ( exist(fname_branch,'file') &  flag_force_create);
save(fname_branch ...
,'trace_E_' ...
,'QR_' ...
,'QC_' ...
,'ZR_' ...
,'ZC_' ...
,'lp_R_' ...
,'lp_R_max' ...
,'p_set' ...
,'p_branch' ...
,'type_code' ...
,'ZR_ij' ...
,'branch_A_r_ij_' ...
,'branch_B_r_ij_' ...
,'E1_ZRmax_label_' ...
,'-append' ...
);
end;%if ( exist(fname_branch,'file') &  flag_force_create);
end;%if ((lp_R_max< -log(p_use)) | ~flag_split);
%%%%%%%%;
end;%if (~exist(fname_branch,'file') |  flag_force_create);

fname_branch = sprintf('%s/branch.mat',dir_out);
pause(1.0);
load(fname_branch);
[~,~,~,tmp_p_opt,tmp_p_emp] = gumbel_fit(max(ZR_(:,2:end),[],1),max(ZR_(:,1)));
p_gbranch = tmp_p_emp; if (p_gbranch<1/n_shuffle); p_gbranch = tmp_p_opt; end;
p_next = 1 - (1-p_prev)*(1-p_gbranch); %<-- q_next = q_prev*q_branch. ;
if strcmp(type_code,'stop');
if (verbose>1); disp(sprintf(' %% %s --> stop',dir_out)); end;
end;%if strcmp(type_code,'stop');
if strcmp(type_code,'split');
if (verbose>1); disp(sprintf(' %% %s --> split (lp_R_max %0.2f lp_R_F_max %0.2f --> p_branch %0.6f p_gbranch %0.6f --> p_next %0.6f)',dir_out,lp_R_max,lp_R_F_max,exp(-lp_R_F_max),p_gbranch,p_next)); end;
dir_A_out = sprintf('%s/dir_A',dir_out); if (~exist(dir_A_out,'dir')); disp(sprintf(' %% mkdir %s',dir_A_out)); mkdir(dir_A_out); end;
E_array_A_r_ij_ = E_array_r_ij_(branch_A_r_ij_); E_array_A_c_ij_ = E_array_c_ij_;
if (verbose>1); disp(sprintf(' %% %s <-- calling (n_sample = %d)',dir_A_out,length(E_array_A_r_ij_))); end;
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_5( ...
 dir_trunk ...
,dir_A_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_A_r_ij_ ...
,E_array_A_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,p_next ...
,flag_force_create ...
);
dir_B_out = sprintf('%s/dir_B',dir_out); if (~exist(dir_B_out,'dir')); disp(sprintf(' %% mkdir %s',dir_B_out)); mkdir(dir_B_out); end;
E_array_B_r_ij_ = E_array_r_ij_(branch_B_r_ij_); E_array_B_c_ij_ = E_array_c_ij_;
if (verbose>1); disp(sprintf(' %% %s <-- calling (n_sample = %d)',dir_B_out,length(E_array_B_r_ij_))); end;
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_5( ...
 dir_trunk ...
,dir_B_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_B_r_ij_ ...
,E_array_B_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,p_next ...
,flag_force_create ...
);
end;%if strcmp(type_code,'split');

end;%if nargout==0; %<-- calculation. ;
