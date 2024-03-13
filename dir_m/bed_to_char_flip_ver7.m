function ...
[ ...
 parameter ...
] = ...
bed_to_char_flip_ver7( ...
 parameter ...
,n_study ...
,study_name_of_branch_s_ ...
,study_name_without_extension_s_ ...
,n_mds_p ...
,n_mds_v ...
,mds_pv__ ...
,mds_fidandiid_p_ ...
,n_famex ...
,famex_fidandiid_p_ ...
);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_study=[]; end; na=na+1;
if (nargin<1+na); study_name_of_branch_s_=[]; end; na=na+1;
if (nargin<1+na); study_name_without_extension_s_=[]; end; na=na+1;
if (nargin<1+na); n_mds_p=[]; end; na=na+1;
if (nargin<1+na); n_mds_v=[]; end; na=na+1;
if (nargin<1+na); mds_pv__=[]; end; na=na+1;
if (nargin<1+na); mds_fidandiid_p_=[]; end; na=na+1;
if (nargin<1+na); n_famex=[]; end; na=na+1;
if (nargin<1+na); famex_fidandiid_p_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'str_output_prefix'); parameter.str_output_prefix=''; end;
if ~isfield(parameter,'dir_trunk'); parameter.dir_trunk=pwd; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
if ~isfield(parameter,'flag_stop'); parameter.flag_stop=0; end;
if ~isfield(parameter,'flag_crop'); parameter.flag_crop=0; end;
if ~isfield(parameter,'flag_bim_ext'); parameter.flag_bim_ext=1; end;
if ~isfield(parameter,'flag_fam_ext'); parameter.flag_fam_ext=1; end;
if ~isfield(parameter,'ent_cutoff'); parameter.ent_cutoff=0.03; end;
if ~isfield(parameter,'maf_cutoff'); parameter.maf_cutoff=0.01; end;
if ~isfield(parameter,'mss_cutoff'); parameter.mss_cutoff=0.02; end;
if ~isfield(parameter,'flag_force_create'); parameter.flag_force_create=0; end;
str_output_prefix = parameter.str_output_prefix;
dir_trunk = parameter.dir_trunk;
flag_verbose = parameter.flag_verbose;
flag_stop = parameter.flag_stop;
flag_crop = parameter.flag_crop;
flag_bim_ext = parameter.flag_bim_ext;
flag_fam_ext = parameter.flag_fam_ext;
ent_cutoff = parameter.ent_cutoff;
maf_cutoff = parameter.maf_cutoff;
mss_cutoff = parameter.mss_cutoff;
flag_force_create = parameter.flag_force_create;
%%%%%%%%;
parameter.n_study = n_study; parameter.n_bin = n_study; %<-- synonyms. ;
parameter.study_name_of_branch_s_ = study_name_of_branch_s_;
parameter.study_name_without_extension_s_ = study_name_without_extension_s_;
parameter.n_mds_p = n_mds_p;
parameter.n_mds_v = n_mds_v;
parameter.n_famex = n_famex;

if (flag_verbose>0); disp(sprintf(' %% [entering_bed_to_char_flip_ver7]')); end;

%%%%%%%%;
% flag_stop: [4=after mds, 3=after bim&fam, 2=after masks, 1=after An&At] ;
%%%%%%%%;
%bitj=16; %<-- not necessary. ;
%bit8=8; %<-- not necessary;
bit4=4; %<-- used for n_patient_bed. ;
nkhr=22; 
male_tag=1;
fema_tag=2;
case_tag=2;
ctrl_tag=1;
%%%%%%%%%%%%%%%%;
% 0: nor <-- both off (homozygous GG). ;
% 1: mss <-- missing. ;
% 2: xor <-- one on one off (heterozygous AG). ;
% 3: and <-- both on (homozygous AA). ;
% 5: not <-- do not use. ;
%%%%%%%%%%%%%%%%;
snp_and_tag=uint8(3);
snp_xor_tag=uint8(2);
snp_nor_tag=uint8(0);
snp_mss_tag=uint8(1);
snp_not_tag=uint8(5);

if (flag_crop==0);
str_tmp = sprintf('%s_maf%.2d',str_output_prefix,floor(100*maf_cutoff));
end;%if (flag_crop==0);
if (flag_crop==1);
str_tmp = sprintf('%s_maf%.2d_crop',str_output_prefix,floor(100*maf_cutoff));
end;%if (flag_crop==0);
str_output_prefix_local = sprintf('%s_',str_tmp);
dir_out = sprintf('%s/dir_%s',dir_trunk,str_tmp); if ~exist(dir_out,'dir'); mkdir(dir_out); end;
if (flag_verbose>0); disp(sprintf(' %% str_output_prefix_local %s ;\n %% dir_out %s ;\n',str_output_prefix_local,dir_out)); end;

flag_exist_all = 1;
tmpchar = sprintf('%s/%schar_bim.ext',dir_out,str_output_prefix_local);
if (flag_verbose>0); if exist(tmpchar,'file'); disp(sprintf(' %% %s yes found',tmpchar)); else; disp(sprintf(' %% %s not found',tmpchar)); end; end;
flag_exist_all = flag_exist_all & exist(tmpchar,'file');
tmpchar = sprintf('%s/%schar_fam.ext',dir_out,str_output_prefix_local);
if (flag_verbose>0); if exist(tmpchar,'file'); disp(sprintf(' %% %s yes found',tmpchar)); else; disp(sprintf(' %% %s not found',tmpchar)); end; end;
flag_exist_all = flag_exist_all & exist(tmpchar,'file');
tmpchar = sprintf('%s/%sA_full_n.char',dir_out,str_output_prefix_local);
if (flag_verbose>0); if exist(tmpchar,'file'); disp(sprintf(' %% %s yes found',tmpchar)); else; disp(sprintf(' %% %s not found',tmpchar)); end; end;
flag_exist_all = flag_exist_all & exist(tmpchar,'file');
tmpchar = sprintf('%s/%sA_full_t.char',dir_out,str_output_prefix_local);
if (flag_verbose>0); if exist(tmpchar,'file'); disp(sprintf(' %% %s yes found',tmpchar)); else; disp(sprintf(' %% %s not found',tmpchar)); end; end;
flag_exist_all = flag_exist_all & exist(tmpchar,'file');
for nstudy=0:n_study-1;
tmpchar = sprintf('%s/%sA_%.2d_n.char',dir_out,str_output_prefix_local,1+nstudy);
if (flag_verbose>0); if exist(tmpchar,'file'); disp(sprintf(' %% %s yes found',tmpchar)); else; disp(sprintf(' %% %s not found',tmpchar)); end; end;
flag_exist_all = flag_exist_all & exist(tmpchar,'file');
tmpchar = sprintf('%s/%sA_%.2d_t.char',dir_out,str_output_prefix_local,1+nstudy);
if (flag_verbose>0); if exist(tmpchar,'file'); disp(sprintf(' %% %s yes found',tmpchar)); else; disp(sprintf(' %% %s not found',tmpchar)); end; end;
flag_exist_all = flag_exist_all & exist(tmpchar,'file');
end;%for nstudy=0:n_study-1;

if ( flag_exist_all);
if (flag_verbose>0); disp(sprintf(' %% all files found in bed_to_char_flip_ver7;')); end;
end;%if ( flag_exist_all);
if (~flag_exist_all);
if (flag_verbose>0); disp(sprintf(' %% some files missing in bed_to_char_flip_ver7;')); end;
end;%if (~flag_exist_all);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ( flag_force_create | ~flag_exist_all);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_crop);
n_snp_crop_s_ = 1*1024 + [1:n_study]; %<-- reduce number of snps for each study. ;
end;%if (flag_crop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Loading data ;')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn_fam_s_ = cell(n_study,1); fp_fam_s_ = cell(n_study,1); 
fn_bim_s_ = cell(n_study,1); fp_bim_s_ = cell(n_study,1); 
fn_bed_s_ = cell(n_study,1); fp_bed_s_ = cell(n_study,1); 
fam_s_ = cell(n_study,1); fam_fidandiid_s_ = cell(n_study,1);
bim_s_ = cell(n_study,1); khrdist_s_ = cell(n_study,1);
n_patient_fam_s_ = zeros(n_study,1); n_patient_bed_s_ = zeros(n_study,1); n_patient_keep_s_ = zeros(n_study,1);
rij_male_s_ = cell(n_study,1);
rij_fema_s_ = cell(n_study,1);
rij_case_s_ = cell(n_study,1);
rij_ctrl_s_ = cell(n_study,1);
rij_keep_s_ = cell(n_study,1);
rindex_keep_s_ = cell(n_study,1);
n_and_s_ = cell(n_study,1); n_xor_s_ = cell(n_study,1); n_nor_s_ = cell(n_study,1); n_mss_s_ = cell(n_study,1);
m_and_s_ = cell(n_study,1); m_xor_s_ = cell(n_study,1); m_nor_s_ = cell(n_study,1); m_mss_s_ = cell(n_study,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Reading fam ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for ns=0:n_study-1;
if (flag_verbose>0); disp(sprintf(' %% ns %d: Reading fam file %s/%s.fam (patients) ;',ns,study_name_of_branch_s_{1+ns},study_name_without_extension_s_{1+ns})); end;
fn_fam_s_{1+ns} = sprintf('%s/%s/%s.fam',dir_trunk,study_name_of_branch_s_{1+ns},study_name_without_extension_s_{1+ns});
fp_fam_s_{1+ns} = fopen(fn_fam_s_{1+ns}); fam_s_{1+ns} = textscan(fp_fam_s_{1+ns},'%s %s %s %s %d %d'); fclose(fp_fam_s_{1+ns});
n_patient_fam_s_(1+ns) = numel(fam_s_{1+ns}{1+0}); 
tmp_rem = mod(n_patient_fam_s_(1+ns),bit4); tmp_quo = (n_patient_fam_s_(1+ns)-tmp_rem)/bit4;
n_patient_bed_s_(1+ns) = tmp_quo + (tmp_rem>0);
rij_male_s_{1+ns} = find(fam_s_{1+ns}{1+4}==male_tag);rij_fema_s_{1+ns} = find(fam_s_{1+ns}{1+4}==fema_tag);
rij_case_s_{1+ns} = find(fam_s_{1+ns}{1+5}==case_tag);rij_ctrl_s_{1+ns} = find(fam_s_{1+ns}{1+5}==ctrl_tag);
for np=1:n_patient_fam_s_(1+ns);
fam_fidandiid_s_{1+ns}{np} = sprintf('%s%s%s',fam_s_{1+ns}{1+0}{np},'&',fam_s_{1+ns}{1+1}{np});
end;%for np=1:n_patient_fam_s_(1+ns);
[~,rij_keep_s_{1+ns}] = setdiff(fam_fidandiid_s_{1+ns},famex_fidandiid_p_,'stable'); %<-- ignores empty famex_fidandiid_p_. ;
n_patient_keep_s_(1+ns) = numel(rij_keep_s_{1+ns});
if (flag_verbose>0);
disp(sprintf(' %% found total: n_patient %d (%d males, %d femas, %d cases, %d ctrls) --> %d' ...
	     ,n_patient_fam_s_(1+ns) ...
	     ,numel(rij_male_s_{1+ns}) ...
	     ,numel(rij_fema_s_{1+ns}) ...
	     ,numel(rij_case_s_{1+ns}) ...
	     ,numel(rij_ctrl_s_{1+ns}) ...
	     ,n_patient_bed_s_(1+ns) ...
			     ));
end;% if (flag_verbose>0);
if (flag_verbose>0); 
disp(sprintf(' %% will keep:   n_patient %d (%d males, %d femas, %d cases, %d ctrls) --> %d' ...
	     ,numel(rij_keep_s_{1+ns}) ...
	     ,numel(intersect(rij_keep_s_{1+ns},rij_male_s_{1+ns})) ...
	     ,numel(intersect(rij_keep_s_{1+ns},rij_fema_s_{1+ns})) ...
	     ,numel(intersect(rij_keep_s_{1+ns},rij_case_s_{1+ns})) ...
	     ,numel(intersect(rij_keep_s_{1+ns},rij_ctrl_s_{1+ns})) ...
	     ,n_patient_bed_s_(1+ns) ...
			     ));
end;% if (flag_verbose>0);
rindex_keep_s_{1+ns} = rij_keep_s_{1+ns}-1;
end;%for ns=0:n_study-1;
n_patient_fam_csum_s_ = cumsum([0;n_patient_fam_s_(:)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Reading bim ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_snp_s_ = zeros(n_study,1);
for ns=0:n_study-1;
if (flag_verbose>0); disp(sprintf(' %% ns %d: Reading bim file %s/%s.bim (snp names) ;',ns,study_name_of_branch_s_{1+ns},study_name_without_extension_s_{1+ns})); end;
fn_bim_s_{1+ns} = sprintf('%s/%s/%s.bim',dir_trunk,study_name_of_branch_s_{1+ns},study_name_without_extension_s_{1+ns});
fp_bim_s_{1+ns} = fopen(fn_bim_s_{1+ns}); bim_s_{1+ns} = textscan(fp_bim_s_{1+ns},'%d %s %d %d %s %s'); fclose(fp_bim_s_{1+ns});
n_snp_s_(1+ns) = numel(bim_s_{1+ns}{1+0});
if (flag_crop);
for nl=1:6; bim_s_{1+ns}{nl} = bim_s_{1+ns}{nl}(1:n_snp_crop_s_(1+ns),:); end;
n_snp_s_(1+ns) = n_snp_crop_s_(1+ns); 
end%if (flag_crop);
tmp=zeros(nkhr,1); for nk=1:nkhr; tmp(nk) = numel(find(bim_s_{1+ns}{1+0}==nk)); end;%for nk=1:nkhr;
khrdist_s_{1+ns} = tmp;
if (flag_verbose>0); disp(sprintf(' %% n_snp_s_(%d) --> chromosome distribution: [%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d]',n_snp_s_(1+ns),khrdist_s_{1+ns})); end;
end;%for ns=0:n_study-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Calculating intersection (i.e., cap) of snps ;')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snp_0id_s_ = cell(n_study,1); snp_pos_s_ = cell(n_study,1); 
for ns=0:n_study-1; snp_0id_s_{1+ns} = bim_s_{1+ns}{1+1}; snp_pos_s_{1+ns} = bim_s_{1+ns}{1+3}; end;%for ns=0:n_study-1;
snp_al1_ = cell(n_study,1); snp_al2_ = cell(n_study,1);
for ns=0:n_study-1; snp_al1_{1+ns} = bim_s_{1+ns}{1+4}; snp_al2_{1+ns} = bim_s_{1+ns}{1+5}; end;%for ns=0:n_study-1;
[snp_cap_] = intersectall(snp_0id_s_); n_snp_cap = numel(snp_cap_);
ij_snp_cap_sl__ = cell(n_study,1);
for ns=0:n_study-1; [~,~,ij_snp_cap_sl__{1+ns}] = intersect(snp_cap_,snp_0id_s_{1+ns},'stable'); end;%for ns=0:n_study-1;
if (flag_verbose>0); disp(sprintf(' %% n_snp_cap %d',n_snp_cap)); end;
disp_flag=0;
if disp_flag;
figure;cla; 
for ns1=0:n_study-1;
for ns2=ns1+0:n_study-1;
subplot(1,2,1); hold on;
plot(1:n_snp_cap,snp_pos_s_{1+ns1}(ij_snp_cap_sl__{1+ns1}),'ro',1:n_snp_cap,snp_pos_s_{1+ns2}(ij_snp_cap_sl__{1+ns2}),'b.'); xlabel('snp index number'); ylabel('snp pos'); title('snp ij vs pos');
subplot(1,2,2); hold on;
plot(snp_pos_s_{1+ns1}(ij_snp_cap_sl__{1+ns1}),snp_pos_s_{1+ns2}(ij_snp_cap_sl__{1+ns2}),'.'); xlabel('snp pos'); ylabel('snp pos'); title('pos vs pos');
end;%for ns2=ns1+0:n_study-1;
end;%for ns1=0:n_study-1;
end;%if disp_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% now checking to see if the snps in each study are unique or appear multiple times. ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for ns=0:n_study-1; 
if (flag_verbose>0); disp(sprintf('study %d: numel(unique(snp_0id_s_{1+ns}))-numel(unique(snp_0id_s_{1+ns})) = %d',1+ns,numel(unique(snp_0id_s_{1+ns}))-numel(snp_0id_s_{1+ns}))); end; 
end;%for ns=0:n_study-1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% now checking to see if the snps in the intersection are consistent across studies. ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
snp_sg1_ls__ = zeros(n_snp_cap,n_study); snp_sg2_ls__ = zeros(n_snp_cap,n_study);
for ns=0:n_study-1;
tmp_sg1_ = cast(cell2mat(snp_al1_{1+ns}(ij_snp_cap_sl__{1+ns})),'double')-'A'+1;
tmp_sg2_ = cast(cell2mat(snp_al2_{1+ns}(ij_snp_cap_sl__{1+ns})),'double')-'A'+1;
% convert T to A;
tmp_ij = find(tmp_sg1_=='T'-'A'+1); tmp_sg1_(tmp_ij) = 'A'-'A'+1;
% convert C to G;
tmp_ij = find(tmp_sg1_=='C'-'A'+1); tmp_sg1_(tmp_ij) = 'G'-'A'+1;
% convert A to -1 and G to +1;
tmp_sg1_ = 2*(tmp_sg1_=='G'-'A'+1)-1;
snp_sg1_ls__(:,1+ns) = tmp_sg1_;
% convert T to A;
tmp_ij = find(tmp_sg2_=='T'-'A'+1); tmp_sg2_(tmp_ij) = 'A'-'A'+1;
% convert C to G;
tmp_ij = find(tmp_sg2_=='C'-'A'+1); tmp_sg2_(tmp_ij) = 'G'-'A'+1;
% convert A to -1 and G to +1;
tmp_sg2_ = 2*(tmp_sg2_=='G'-'A'+1)-1;
snp_sg2_ls__(:,1+ns) = tmp_sg2_;
end;%for ns=0:n_study-1;
% Now assert that al1 and al2 are complements. ;
if (flag_verbose>0); disp(sprintf(' %% norm(sg1+sg2) %0.16f',norm(snp_sg1_ls__+snp_sg2_ls__,'fro'))); end;
% Now compute hamming distance;
flag_disp=0;
if flag_disp;
disp((transpose(snp_sg1_ls__)*snp_sg1_ls__ + n_snp_cap)/2/n_snp_cap);
imagesc(((transpose(snp_sg1_ls__)*snp_sg1_ls__ + n_snp_cap)/2/n_snp_cap),[0,1]); colorbar;
set(gca,'XTick',1:n_study,'XTickLabel',study_name_without_extension_s_,'ticklabelinterpreter','none'); xtickangle(90);
set(gca,'YTick',1:n_study,'YTickLabel',study_name_without_extension_s_,'ticklabelinterpreter','none');
end;%if flag_disp;
if (flag_verbose>0); disp(sprintf(' %% fraction of snps with mismatch: %d/%d = %0.4f',numel(find(abs(sum(snp_sg1_ls__,2))~=n_study)),numel(snp_sg1_ls__),numel(find(abs(sum(snp_sg1_ls__,2))~=n_study))/numel(snp_sg1_ls__))); end;
if (flag_verbose>0); disp(sprintf(' %% fraction of mismatches (within %d-x-%d data-array) %d/%d = %0.4f',n_snp_cap,n_study,sum((n_study-abs(sum(snp_sg1_ls__(:,:),2)))/2),(numel(snp_sg1_ls__(:,1))*n_study),sum((n_study-abs(sum(snp_sg1_ls__(:,:),2)))/2)/(numel(snp_sg1_ls__(:,1))*n_study))); end;
if (flag_verbose>0); disp(sprintf(' %% number of sign flips per study: ')); disp(num2str(sum(snp_sg1_ls__==-1))); end;
if (flag_verbose>0); disp(sprintf(' %% number of sign keeps per study: ')); disp(num2str(sum(snp_sg1_ls__==+1))); end;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Reading bed. ;')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ns=0:n_study-1;
if (flag_verbose>0); disp(sprintf(' %% ns %d: Reading bed file %s/%s.bed (binarized snp data) ;',ns,study_name_of_branch_s_{1+ns},study_name_without_extension_s_{1+ns})); end;
fn_bed_s_{1+ns} = sprintf('%s/%s/%s.bed',dir_trunk,study_name_of_branch_s_{1+ns},study_name_without_extension_s_{1+ns});
fp_bed_s_{1+ns} = fopen(fn_bed_s_{1+ns}); 
key1 = fread(fp_bed_s_{1+ns},1,'uint8'); key2 = fread(fp_bed_s_{1+ns},1,'uint8'); key3 = fread(fp_bed_s_{1+ns},1,'uint8'); assert(key1==108); assert(key2==27); assert(key3==1);
flag_ver=0;
if flag_ver==0;
if (flag_verbose>0); disp(sprintf(' %% version 0: simultaneous read (assuming sufficient ram)...')); end; tic;
% These are the sums -- over patients -- of the and, xor, nor and mss symbols for each snp. ;
n_and_s_{1+ns} = zeros(1,n_snp_s_(1+ns));n_xor_s_{1+ns} = zeros(1,n_snp_s_(1+ns));n_nor_s_{1+ns} = zeros(1,n_snp_s_(1+ns));n_mss_s_{1+ns} = zeros(1,n_snp_s_(1+ns));
% These are the sums -- over snps -- of the and, xor, nor and mss symbols for each patient. ; only m_mss_s_{1+ns} is used. ;
%m_and_s_{1+ns} = zeros(1,4*n_patient_bed_s_(1+ns));m_xor_s_{1+ns} = zeros(1,4*n_patient_bed_s_(1+ns));m_nor_s_{1+ns} = zeros(1,4*n_patient_bed_s_(1+ns));
m_mss_s_{1+ns} = zeros(1,4*n_patient_bed_s_(1+ns));
%%%%%%%%%%%%%%%%;
% uchar_11111111 in [0,..,255];
% uchar_11000000 in [0,1,2,3];
% uchar_00110000 in [0,1,2,3];
% uchar_00001100 in [0,1,2,3];
% uchar_00000011 in [0,1,2,3];
%%%%%%%%%%%%%%%%;
uchar_11111111 = transpose(fread(fp_bed_s_{1+ns},[n_patient_bed_s_(1+ns),n_snp_s_(1+ns)],'uint8=>uint8'));
t_tmp = toc; if (flag_verbose>0);disp(sprintf(' %% read input: time %0.2f seconds',t_tmp)); end; 
if (flag_verbose>0); disp(sprintf(' %% reading 4 allele pairs per byte... ')); end; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% these commands use mod, which is very slow. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%uchar_00000011 = mod(uchar_11111111,  4);
%uchar_00001111 = mod(uchar_11111111, 16); uchar_00001100 = idivide(uchar_00001111 - uchar_00000011, 4) ; 
%uchar_00111111 = mod(uchar_11111111, 64); uchar_00110000 = idivide(uchar_00111111 - uchar_00001111,16) ; clear uchar_00001111;
%                                          uchar_11000000 = idivide(uchar_11111111 - uchar_00111111,64) ; clear uchar_00111111;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% these commands use bitshifts, which are much faster. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
uchar_00000011 = bitshift(bitshift(uchar_11111111,+6),-6);
uchar_00001100 = bitshift(bitshift(bitshift(uchar_11111111,-2),+6),-6);
uchar_00110000 = bitshift(bitshift(bitshift(uchar_11111111,-4),+6),-6);
uchar_11000000 = bitshift(uchar_11111111,-6);
clear uchar_11111111;
t_tmp = toc; if (flag_verbose>0);disp(sprintf(' %% bitshift: time %0.2f seconds',t_tmp)); end; 
if (flag_verbose>0); disp(sprintf(' %% reading 4 tags per allele pair... ')); end; tic;
%%%%%%%%%%%%%%%%;
% removing patients from famex. ;
%%%%%%%%%%%%%%%%;
stride1 = 1:4:4*n_patient_bed_s_(1+ns); [~,tmp_ij] = setdiff(stride1,1+rindex_keep_s_{1+ns},'stable'); uchar_00000011(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride2 = 2:4:4*n_patient_bed_s_(1+ns); [~,tmp_ij] = setdiff(stride2,1+rindex_keep_s_{1+ns},'stable'); uchar_00001100(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride3 = 3:4:4*n_patient_bed_s_(1+ns); [~,tmp_ij] = setdiff(stride3,1+rindex_keep_s_{1+ns},'stable'); uchar_00110000(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride4 = 4:4:4*n_patient_bed_s_(1+ns); [~,tmp_ij] = setdiff(stride4,1+rindex_keep_s_{1+ns},'stable'); uchar_11000000(:,tmp_ij)=snp_not_tag; clear tmp_ij;
clear uchar_tmp;
%%%%%%%%%%%%%%%%;
% define ;
%%%%%%%%%%%%%%%%;
uchar_tmp = uchar_00000011==snp_and_tag; sum1_and_00000011 = sum(uchar_tmp,1); sum2_and_00000011 = sum(uchar_tmp,2);
uchar_tmp = uchar_00000011==snp_xor_tag; sum1_xor_00000011 = sum(uchar_tmp,1); sum2_xor_00000011 = sum(uchar_tmp,2);
uchar_tmp = uchar_00000011==snp_nor_tag; sum1_nor_00000011 = sum(uchar_tmp,1); sum2_nor_00000011 = sum(uchar_tmp,2);
uchar_tmp = uchar_00000011==snp_mss_tag; sum1_mss_00000011 = sum(uchar_tmp,1); sum2_mss_00000011 = sum(uchar_tmp,2);
clear uchar_tmp;
uchar_tmp = uchar_00001100==snp_and_tag; sum1_and_00001100 = sum(uchar_tmp,1); sum2_and_00001100 = sum(uchar_tmp,2);
uchar_tmp = uchar_00001100==snp_xor_tag; sum1_xor_00001100 = sum(uchar_tmp,1); sum2_xor_00001100 = sum(uchar_tmp,2);
uchar_tmp = uchar_00001100==snp_nor_tag; sum1_nor_00001100 = sum(uchar_tmp,1); sum2_nor_00001100 = sum(uchar_tmp,2);
uchar_tmp = uchar_00001100==snp_mss_tag; sum1_mss_00001100 = sum(uchar_tmp,1); sum2_mss_00001100 = sum(uchar_tmp,2);
clear uchar_tmp;
uchar_tmp = uchar_00110000==snp_and_tag; sum1_and_00110000 = sum(uchar_tmp,1); sum2_and_00110000 = sum(uchar_tmp,2);
uchar_tmp = uchar_00110000==snp_xor_tag; sum1_xor_00110000 = sum(uchar_tmp,1); sum2_xor_00110000 = sum(uchar_tmp,2);
uchar_tmp = uchar_00110000==snp_nor_tag; sum1_nor_00110000 = sum(uchar_tmp,1); sum2_nor_00110000 = sum(uchar_tmp,2);
uchar_tmp = uchar_00110000==snp_mss_tag; sum1_mss_00110000 = sum(uchar_tmp,1); sum2_mss_00110000 = sum(uchar_tmp,2);
clear uchar_tmp;
uchar_tmp = uchar_11000000==snp_and_tag; sum1_and_11000000 = sum(uchar_tmp,1); sum2_and_11000000 = sum(uchar_tmp,2);
uchar_tmp = uchar_11000000==snp_xor_tag; sum1_xor_11000000 = sum(uchar_tmp,1); sum2_xor_11000000 = sum(uchar_tmp,2);
uchar_tmp = uchar_11000000==snp_nor_tag; sum1_nor_11000000 = sum(uchar_tmp,1); sum2_nor_11000000 = sum(uchar_tmp,2);
uchar_tmp = uchar_11000000==snp_mss_tag; sum1_mss_11000000 = sum(uchar_tmp,1); sum2_mss_11000000 = sum(uchar_tmp,2);
clear uchar_tmp;
clear uchar_00000011;clear uchar_00001100;clear uchar_00110000;clear uchar_11000000;
t_tmp = toc; if (flag_verbose>0);disp(sprintf(' %% sums: time %0.2f seconds',t_tmp)); end; 
if (flag_verbose>0); disp(sprintf(' %% totalling sum2...')); end; tic;
n_and_s_{1+ns} = sum2_and_00000011 + sum2_and_00001100 + sum2_and_00110000 + sum2_and_11000000;
n_xor_s_{1+ns} = sum2_xor_00000011 + sum2_xor_00001100 + sum2_xor_00110000 + sum2_xor_11000000;
n_nor_s_{1+ns} = sum2_nor_00000011 + sum2_nor_00001100 + sum2_nor_00110000 + sum2_nor_11000000;
n_mss_s_{1+ns} = sum2_mss_00000011 + sum2_mss_00001100 + sum2_mss_00110000 + sum2_mss_11000000;
t_tmp = toc; if (flag_verbose>0);disp(sprintf(' %% n_xxx: time %0.2f seconds',t_tmp)); end; 
if (flag_verbose>0); disp(sprintf(' %% totalling sum1...')); end; tic;
%m_and_s_{1+ns}(stride1) = m_and_s_{1+ns}(stride1) + sum1_and_00000011;
%m_and_s_{1+ns}(stride2) = m_and_s_{1+ns}(stride2) + sum1_and_00001100;
%m_and_s_{1+ns}(stride3) = m_and_s_{1+ns}(stride3) + sum1_and_00110000;
%m_and_s_{1+ns}(stride4) = m_and_s_{1+ns}(stride4) + sum1_and_11000000;
%m_xor_s_{1+ns}(stride1) = m_xor_s_{1+ns}(stride1) + sum1_xor_00000011;
%m_xor_s_{1+ns}(stride2) = m_xor_s_{1+ns}(stride2) + sum1_xor_00001100;
%m_xor_s_{1+ns}(stride3) = m_xor_s_{1+ns}(stride3) + sum1_xor_00110000;
%m_xor_s_{1+ns}(stride4) = m_xor_s_{1+ns}(stride4) + sum1_xor_11000000;
%m_nor_s_{1+ns}(stride1) = m_nor_s_{1+ns}(stride1) + sum1_nor_00000011;
%m_nor_s_{1+ns}(stride2) = m_nor_s_{1+ns}(stride2) + sum1_nor_00001100;
%m_nor_s_{1+ns}(stride3) = m_nor_s_{1+ns}(stride3) + sum1_nor_00110000;
%m_nor_s_{1+ns}(stride4) = m_nor_s_{1+ns}(stride4) + sum1_nor_11000000;
m_mss_s_{1+ns}(stride1) = m_mss_s_{1+ns}(stride1) + sum1_mss_00000011;
m_mss_s_{1+ns}(stride2) = m_mss_s_{1+ns}(stride2) + sum1_mss_00001100;
m_mss_s_{1+ns}(stride3) = m_mss_s_{1+ns}(stride3) + sum1_mss_00110000;
m_mss_s_{1+ns}(stride4) = m_mss_s_{1+ns}(stride4) + sum1_mss_11000000;
t_tmp = toc; if (flag_verbose>0);disp(sprintf(' %% m_xxx: time %0.2f seconds',t_tmp)); end; tic;
clear sum1_and_00000011;clear sum1_and_00001100;clear sum1_and_00110000;clear sum1_and_11000000;
clear sum1_xor_00000011;clear sum1_xor_00001100;clear sum1_xor_00110000;clear sum1_xor_11000000;
clear sum1_nor_00000011;clear sum1_nor_00001100;clear sum1_nor_00110000;clear sum1_nor_11000000;
clear sum1_mss_00000011;clear sum1_mss_00001100;clear sum1_mss_00110000;clear sum1_mss_11000000;
clear sum2_and_00000011;clear sum2_and_00001100;clear sum2_and_00110000;clear sum2_and_11000000;
clear sum2_xor_00000011;clear sum2_xor_00001100;clear sum2_xor_00110000;clear sum2_xor_11000000;
clear sum2_nor_00000011;clear sum2_nor_00001100;clear sum2_nor_00110000;clear sum2_nor_11000000;
clear sum2_mss_00000011;clear sum2_mss_00001100;clear sum2_mss_00110000;clear sum2_mss_11000000;
clear stride1 stride2 stride3 stride4;
end;%if flag_ver==0;
if flag_ver==1;
if (flag_verbose>0); disp(sprintf(' %% Warning! version 1: sequential read not implemented')); end;
disp(' %% returning'); return;
end;%if flag_ver==1;
fclose(fp_bed_s_{1+ns});
end;%for ns=0:n_study-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Calculating relative-entropy (i.e., KL-divergence) for each snp within each study. ;')); end;
if (flag_verbose>0); disp(sprintf(' %% If we assume that a snp has minor allele frequency q and major allele frequency p, ;')); end;
if (flag_verbose>0); disp(sprintf(' %% then the distribution of (n_and,n_xor,n_nor) should be proportional to: ;')); end;
if (flag_verbose>0); disp(sprintf(' %% (p*p, 2*p*q, q*q), respectively. ;')); end;
if (flag_verbose>0); disp(sprintf(' %% If the true distribution is given by (frq_and,frq_xor,frq_nor) = (n_and,n_xor,n_nor)/n_patient_keep, ;')); end;
if (flag_verbose>0); disp(sprintf(' %% then the relative entropy between the true distribution and the expected distribution is: ;')); end;
if (flag_verbose>0); disp(sprintf(' %% I = frq_and*log(frq_and/p^2) + frq_xor*log(frq_xor/(2*p*q)) + frq_nor*log(frq_nor/q^2). ;')); end;
if (flag_verbose>0); disp(sprintf(' %% The derivative of I with respect to p is: ;')); end;
if (flag_verbose>0); disp(sprintf(' %% dI/dp = -2*frq_and/p - (q-p)*frq_xor/(p*q) + 2*frq_nor/q. ;')); end;
if (flag_verbose>0); disp(sprintf(' %% Setting this to 0, we see that the optimal p and q are given by: ;')); end;
if (flag_verbose>0); disp(sprintf(' %% 2*p*frq_nor - 2*q*frq_and = (q-p)*frq_xor, ;')); end;
if (flag_verbose>0); disp(sprintf(' %% p*(2*frq_nor+2*frq_and+2*frq_xor) = 2*frq_and + frq_xor, ;')); end;
if (flag_verbose>0); disp(sprintf(' %% or simply: ;')); end;
if (flag_verbose>0); disp(sprintf(' %% p_opt = frq_and + 0.5*frq_xor, ')); end;
if (flag_verbose>0); disp(sprintf(' %% q_opt = frq_nor + 0.5*frq_xor. ,')); end;
if (flag_verbose>0); disp(sprintf(' %% I_opt = frq_and*log(frq_and/p_opt^2) + frq_xor*log(frq_xor/(2*p_opt*q_opt)) + frq_nor*log(frq_nor/q_opt^2). ;')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ns=0:n_study-1;
frq_and_ = n_and_s_{1+ns}/n_patient_keep_s_(1+ns); frq_xor_ = n_xor_s_{1+ns}/n_patient_keep_s_(1+ns); frq_nor_ = n_nor_s_{1+ns}/n_patient_keep_s_(1+ns);
p_opt_s_{1+ns} = frq_and_ + 0.5*frq_xor_; q_opt_s_{1+ns} = frq_nor_ + 0.5*frq_xor_;
tmp_and_ = frq_and_.*log(frq_and_./(p_opt_s_{1+ns}.^2)); tmp_and_(find(~isfinite(tmp_and_)))=0;
tmp_xor_ = frq_xor_.*log(frq_xor_./(2*p_opt_s_{1+ns}.*q_opt_s_{1+ns})); tmp_xor_(find(~isfinite(tmp_xor_)))=0;
tmp_nor_ = frq_nor_.*log(frq_nor_./(q_opt_s_{1+ns}.^2)); tmp_nor_(find(~isfinite(tmp_nor_)))=0;
I_opt_s_{1+ns} = tmp_and_ + tmp_xor_ + tmp_nor_ ;
clear tmp_and_ tmp_xor_ tmp_nor_ ;
if (flag_verbose>0); disp(sprintf(' %% ns %d: [05 25 50 75 95] percentile: ',ns)); end;
tmp = [prctile(p_opt_s_{1+ns}, 5) , prctile(p_opt_s_{1+ns},25) , prctile(p_opt_s_{1+ns},50) , prctile(p_opt_s_{1+ns},75) , prctile(p_opt_s_{1+ns},95) ] ; 
if (flag_verbose>0); disp(sprintf(' %% p_opt_s_ [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp)); end;
tmp = [prctile(q_opt_s_{1+ns}, 5) , prctile(q_opt_s_{1+ns},25) , prctile(q_opt_s_{1+ns},50) , prctile(q_opt_s_{1+ns},75) , prctile(q_opt_s_{1+ns},95) ] ; 
if (flag_verbose>0); disp(sprintf(' %% q_opt_s_ [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp)); end;
tmp = [prctile(I_opt_s_{1+ns}, 5) , prctile(I_opt_s_{1+ns},25) , prctile(I_opt_s_{1+ns},50) , prctile(I_opt_s_{1+ns},75) , prctile(I_opt_s_{1+ns},95) ] ; 
if (flag_verbose>0); disp(sprintf(' %% I_opt_s_ [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp)); end;
end;%for ns=0:n_study-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Calculating p_opt_s_, q_opt_s_ and entropy I_opt_s_ for each snp across studies. ;')); end;
if (flag_verbose>0); disp(sprintf(' %% Note that we must ''flip'' the and-vs-nor values for snp-study combinations that have snp_sg1_ls__(1+nsnp,1+ns)==-1. ;')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_and_tot_l_ = zeros(n_snp_cap,1); n_xor_tot_l_ = zeros(n_snp_cap,1); n_nor_tot_l_ = zeros(n_snp_cap,1); n_mss_tot_l_ = zeros(n_snp_cap,1); I_opt_max_l_ = zeros(n_snp_cap,1);
for ns=0:n_study-1; 
n_and_tot_l_ = n_and_tot_l_ + n_and_s_{1+ns}(ij_snp_cap_sl__{1+ns}).*(snp_sg1_ls__(:,1+ns)==+1) + n_nor_s_{1+ns}(ij_snp_cap_sl__{1+ns}).*(snp_sg1_ls__(:,1+ns)==-1) ;
n_xor_tot_l_ = n_xor_tot_l_ + n_xor_s_{1+ns}(ij_snp_cap_sl__{1+ns}); 
n_nor_tot_l_ = n_nor_tot_l_ + n_nor_s_{1+ns}(ij_snp_cap_sl__{1+ns}).*(snp_sg1_ls__(:,1+ns)==+1) + n_and_s_{1+ns}(ij_snp_cap_sl__{1+ns}).*(snp_sg1_ls__(:,1+ns)==-1) ;
n_mss_tot_l_ = n_mss_tot_l_ + n_mss_s_{1+ns}(ij_snp_cap_sl__{1+ns}); 
I_opt_max_l_ = max(I_opt_max_l_,I_opt_s_{1+ns}(ij_snp_cap_sl__{1+ns}));
end;%for ns=0:n_study-1;
frq_and_tot_l_ = n_and_tot_l_ / sum(n_patient_keep_s_); frq_xor_tot_l_ = n_xor_tot_l_ / sum(n_patient_keep_s_); frq_nor_tot_l_ = n_nor_tot_l_ / sum(n_patient_keep_s_); frq_mss_tot_l_ = n_mss_tot_l_ / sum(n_patient_keep_s_);
fr2_and_tot_l_ = frq_and_tot_l_; fr2_xor_tot_l_ = frq_xor_tot_l_; fr2_nor_tot_l_ = frq_nor_tot_l_; fr2_mode_l_ = zeros(n_snp_cap,1);
% account for eventually replacing missing data for each snp with mode taken across the patients (for that snp). ;
tmp_ij = find(frq_and_tot_l_ >frq_xor_tot_l_ & frq_and_tot_l_ >frq_nor_tot_l_); fr2_mode_l_(tmp_ij) = snp_and_tag; fr2_and_tot_l_(tmp_ij) = fr2_and_tot_l_(tmp_ij) + frq_mss_tot_l_(tmp_ij);
tmp_ij = find(frq_xor_tot_l_>=frq_and_tot_l_ & frq_xor_tot_l_>=frq_nor_tot_l_); fr2_mode_l_(tmp_ij) = snp_xor_tag; fr2_xor_tot_l_(tmp_ij) = fr2_xor_tot_l_(tmp_ij) + frq_mss_tot_l_(tmp_ij);
tmp_ij = find(frq_nor_tot_l_>=frq_and_tot_l_ & frq_nor_tot_l_ >frq_xor_tot_l_); fr2_mode_l_(tmp_ij) = snp_nor_tag; fr2_nor_tot_l_(tmp_ij) = fr2_nor_tot_l_(tmp_ij) + frq_mss_tot_l_(tmp_ij);
clear tmp_ij;
p_opt_tot_l_ = fr2_and_tot_l_ + 0.5*fr2_xor_tot_l_;
q_opt_tot_l_ = fr2_nor_tot_l_ + 0.5*fr2_xor_tot_l_;
tmp_nor_ = fr2_nor_tot_l_.*log(fr2_nor_tot_l_./(q_opt_tot_l_.^2)); tmp_nor_(find(~isfinite(tmp_nor_)))=0;
tmp_xor_ = fr2_xor_tot_l_.*log(fr2_xor_tot_l_./(2*p_opt_tot_l_.*q_opt_tot_l_)); tmp_xor_(find(~isfinite(tmp_xor_)))=0;
tmp_and_ = fr2_and_tot_l_.*log(fr2_and_tot_l_./(p_opt_tot_l_.^2)); tmp_and_(find(~isfinite(tmp_and_)))=0;
I_opt_tot_l_ = tmp_nor_ + tmp_xor_ + tmp_and_ ;
clear tmp_nor_ tmp_xor_ tmp_and_ ;
if (flag_verbose>0); disp(sprintf(' %% all studies: [05 25 50 75 95] percentile: ')); end;
tmp = [prctile(p_opt_tot_l_, 5) , prctile(p_opt_tot_l_,25) , prctile(p_opt_tot_l_,50) , prctile(p_opt_tot_l_,75) , prctile(p_opt_tot_l_,95) ] ; 
if (flag_verbose>0); disp(sprintf(' %% p_opt_tot_l_ [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp)); end;
tmp = [prctile(q_opt_tot_l_, 5) , prctile(q_opt_tot_l_,25) , prctile(q_opt_tot_l_,50) , prctile(q_opt_tot_l_,75) , prctile(q_opt_tot_l_,95) ] ; 
if (flag_verbose>0); disp(sprintf(' %% q_opt_tot_l_ [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp)); end;
tmp = [prctile(I_opt_tot_l_, 5) , prctile(I_opt_tot_l_,25) , prctile(I_opt_tot_l_,50) , prctile(I_opt_tot_l_,75) , prctile(I_opt_tot_l_,95) ] ; 
if (flag_verbose>0); disp(sprintf(' %% I_opt_tot_l_ [%0.2f,%0.2f,%0.2f,%0.2f,%0.2f]',tmp)); end;

%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Now we estimate the sum (i.e., the dosage). ')); end;
%%%%%%%%%%%%%%%%;
fr2_sum_tot_l_ = 0*fr2_nor_tot_l_ + 1*fr2_xor_tot_l_ + 2*fr2_and_tot_l_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Calculating subset of col-intersection which satisfies cutoffs. ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% This calculation does not involve patients which are excluded via famex. ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ij_snp_cap_sub_ = find(I_opt_tot_l_ < ent_cutoff & I_opt_max_l_ < ent_cutoff & min(p_opt_tot_l_,q_opt_tot_l_) > maf_cutoff & frq_mss_tot_l_ < mss_cutoff); n_snp_cap_sub = numel(ij_snp_cap_sub_);
if (flag_verbose>0); disp(sprintf(' %% retaining %d/%d snps',n_snp_cap_sub,n_snp_cap)); end;
ij_snp_cap_sub_sl__ = cell(n_study,1);
for ns=0:n_study-1; ij_snp_cap_sub_sl__{1+ns} = ij_snp_cap_sl__{1+ns}(ij_snp_cap_sub_); end;%for ns=0:n_study-1;
disp_flag=0;
if disp_flag;
figure;cla; 
for ns1=0:n_study-1;
for ns2=ns1+0:n_study-1;
subplot(1,2,1); hold on;
plot(1:n_snp_cap_sub,snp_pos_s_{1+ns1}(ij_snp_cap_sub_sl__{1+ns1}),'ro',1:n_snp_cap_sub,snp_pos_s_{1+ns2}(ij_snp_cap_sub_sl__{1+ns2}),'b.'); xlabel('snp index number'); ylabel('snp pos'); title('snp cij vs pos');
subplot(1,2,2); hold on;
plot(snp_pos_s_{1+ns1}(ij_snp_cap_sub_sl__{1+ns1}),snp_pos_s_{1+ns2}(ij_snp_cap_sub_sl__{1+ns2}),'.'); xlabel('snp pos'); ylabel('snp pos'); title('pos vs pos');
end;%for ns2=ns1+0:n_study-1;
end;%for ns1=0:n_study-1;
end;%if disp_flag;

%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Counting frequencies of each column (i.e., each retained snp across studies) ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Note that this calculation involves fr2_and_tot_l_ and fr2_nor_tot_l_, ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% both of which have been ''flipped'' based on the modification of n_and_tot_l_ ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% and n_nor_tot_l_ using snp_sg1_ls__. ; ')); end;
%%%%%%%%%%%%%%%%;
brk=3; col_fr2 = n_study+1; col_sort_frwd = n_study+2; col_srt_back = n_study+3;
cij_array__ = zeros(3*n_snp_cap_sub,n_study+brk);
cij_sum_array__ = zeros(1*n_snp_cap_sub,n_study+brk); %<-- used for dosage. ;
cij_and = 0*n_snp_cap_sub + transpose(1:n_snp_cap_sub);
cij_xor = 1*n_snp_cap_sub + transpose(1:n_snp_cap_sub);
cij_nor = 2*n_snp_cap_sub + transpose(1:n_snp_cap_sub);
cij_sum = 0*n_snp_cap_sub + transpose(1:n_snp_cap_sub);
for ns=0:n_study-1; 
cij_array__(cij_and,1+ns) = ij_snp_cap_sub_sl__{1+ns}; 
cij_array__(cij_xor,1+ns) = ij_snp_cap_sub_sl__{1+ns}; 
cij_array__(cij_nor,1+ns) = ij_snp_cap_sub_sl__{1+ns}; 
cij_sum_array__(cij_sum,1+ns) = ij_snp_cap_sub_sl__{1+ns}; 
end;% for ns=0:n_study-1;
cij_array__(cij_and,col_fr2) = fr2_and_tot_l_(ij_snp_cap_sub_);
cij_array__(cij_xor,col_fr2) = fr2_xor_tot_l_(ij_snp_cap_sub_);
cij_array__(cij_nor,col_fr2) = fr2_nor_tot_l_(ij_snp_cap_sub_);
[~,tmp_r] = sort(cij_array__(:,col_fr2),'ascend');
cij_array__(:,col_sort_frwd) = tmp_r;
[~,tmp_s] = sort(tmp_r,'ascend');
cij_array__(:,col_srt_back) = tmp_s;
cij_and_rpl_l_ = cij_array__(cij_and,col_srt_back);
cij_xor_rpl_l_ = cij_array__(cij_xor,col_srt_back);
cij_nor_rpl_l_ = cij_array__(cij_nor,col_srt_back);
clear tmp_r tmp_s ;
%%%%;
cij_sum_array__(cij_sum,col_fr2) = fr2_sum_tot_l_(ij_snp_cap_sub_);
[~,tmp_r] = sort(cij_sum_array__(:,col_fr2),'ascend');
cij_sum_array__(:,col_sort_frwd) = tmp_r;
[~,tmp_s] = sort(tmp_r,'ascend');
cij_sum_array__(:,col_srt_back) = tmp_s;
cij_sum_rpl_l_ = cij_sum_array__(cij_sum,col_srt_back);
clear tmp_r tmp_s ;
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Calculating subset of rows which satisfies cutoffs. ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% This calculation does involve all the snps when calculating missing fraction, ;')); end;
if (flag_verbose>0); disp(sprintf(' %% including those which have been excluded above. ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rij_cap_sp__ = cell(n_study,1);
rindex_cap_sp__ = cell(n_study,1);
n_rij_cap_s_ = zeros(n_study,1);
for ns=0:n_study-1; 
rij_tmp_ = zeros(1,n_patient_fam_s_(1+ns)); rij_tmp_(rij_keep_s_{1+ns})=1;
rij_cap_sp__{1+ns} = find((m_mss_s_{1+ns}(1,1:n_patient_fam_s_(1+ns))/n_snp_s_(1+ns) < mss_cutoff) & rij_tmp_); 
n_rij_cap_s_(1+ns) = numel(rij_cap_sp__{1+ns});
rindex_cap_sp__{1+ns} = rij_cap_sp__{1+ns}-1;
if (flag_verbose>0); disp(sprintf(' %% ns %d: retaining %d/%d/%d patients',ns,n_rij_cap_s_(1+ns),n_patient_keep_s_(1+ns),n_patient_fam_s_(1+ns))); end;
end;%for ns=0:n_study-1;
clear rij_tmp_ ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Grabbing mds_pv__ components ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrows_T_full=0; ncols_Tn = 1 + size(mds_pv__,2);
for ns=0:n_study-1;
tmp_fidandiid_p_ = cell(n_rij_cap_s_(1+ns),1);
for np=0:n_rij_cap_s_(1+ns)-1; 
tmp_fidandiid_p_{1+np} = sprintf('%s%s%s',fam_s_{1+ns}{1+0}{1+rindex_cap_sp__{1+ns}(1+np)},'&',fam_s_{1+ns}{1+1}{1+rindex_cap_sp__{1+ns}(1+np)}); 
end;%for np=0:n_rij_cap_s_(1+ns)-1;
if ~isempty(mds_fidandiid_p_);
[~,tmp_ij_r_,tmp_ij_m_] = intersect(tmp_fidandiid_p_,mds_fidandiid_p_,'stable');
end;%if ~isempty(mds_fidandiid_p_);
if  isempty(mds_fidandiid_p_);
tmp_ij_r_ = 1:n_rij_cap_s_(1+ns);
tmp_ij_m_ = n_patient_fam_csum_s_(1+ns) + rij_cap_sp__{1+ns}; %<-- assume mds_pv__ is aligned with original stack of fam files. ;
end;%if  isempty(mds_fidandiid_p_);
if numel(tmp_ij_r_)<n_rij_cap_s_(1+ns); 
if (flag_verbose>0); disp(sprintf(' %% Warning! patients not found in mds_fidandiid_p_: ns %d, tmp_ij_r_ length %d/%d',ns,numel(tmp_ij_r_),n_rij_cap_s_(1+ns))); end; 
rij_cap_sp__{1+ns} = rij_cap_sp__{1+ns}(tmp_ij_r_);
rindex_cap_sp__{1+ns} = rindex_cap_sp__{1+ns}(tmp_ij_r_);
else; if (flag_verbose>0); disp(sprintf(' %% ns %d, tmp_ij_r_ %d == n_rij_cap_s_(1+ns) %d',ns,numel(tmp_ij_r_),n_rij_cap_s_(1+ns))); end;
end;%if numel(tmp_ij_r_)<n_rij_cap_s_(1+ns); 
if (size(mds_pv__,1)<max(tmp_ij_m_)); disp(sprintf(' %% Warning, mds_pv__ is size (%d,%d), too small in bed_to_char_flip_ver7, returning',size(mds_pv__))); return; end;
T_{1+ns} = mds_pv__(tmp_ij_m_,:);
nrows_T_full = nrows_T_full + numel(tmp_ij_r_);
end;%for ns=0:n_study-1;
clear tmp_ij_r_ tmp_ij_m_ ;
T_full_ = zeros(nrows_T_full,size(mds_pv__,2));
nrows_T_sum=0;
for ns=0:n_study-1;
T_full_(nrows_T_sum + (1:size(T_{1+ns},1)),1:size(T_{1+ns},2)) = T_{1+ns};
nrows_T_sum = nrows_T_sum + size(T_{1+ns},1);
end;%for ns=0:n_study-1;
T_med_ = median(T_full_);
T_full_ = T_full_ - repmat(T_med_,nrows_T_sum,1);
for ns=0:n_study-1;
T_{1+ns} = T_{1+ns} - repmat(T_med_,n_rij_cap_s_(1+ns),1);
end;%for ns=0:n_study-1;
mc_T = ones(ncols_Tn,1);
tmpchar = sprintf('%s/%smc_T.char',dir_out,str_output_prefix_local); 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(mc_T(:)>0,tmpchar); 
end;%if (flag_force_create | ~exist(tmpchar,'file'));
tmpchar = sprintf('%s/%sT_full_n.char',dir_out,str_output_prefix_local);
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress([ones(nrows_T_sum,1) , (T_full_>0)],tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
tmpchar = sprintf('%s/%sT_full_t.char',dir_out,str_output_prefix_local);
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(transpose([ones(nrows_T_sum,1) , (T_full_>0)]),tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
for ns=0:n_study-1;
tmpchar = sprintf('%s/%sT_%.2d_n.char',dir_out,str_output_prefix_local,1+ns);
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress([ones(n_rij_cap_s_(1+ns),1) , (T_{1+ns}>0)],tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
tmpchar = sprintf('%s/%sT_%.2d_t.char',dir_out,str_output_prefix_local,1+ns);
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(transpose([ones(n_rij_cap_s_(1+ns),1) , (T_{1+ns}>0)]),tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
end;%for ns=0:n_study-1;

if (size(T_full_,2)>=2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% creating mds_m2rx ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_m2rz_med_ = transpose(adi_median(transpose(T_full_(:,1:2)),1e-9,max(0,flag_verbose-1)));
if (flag_verbose>0); disp(sprintf('T_m2rz_med_: %0.6f,%0.6f',T_m2rz_med_)); end;
T_m2rz_full_ = T_full_(:,1:2) - repmat(T_m2rz_med_,nrows_T_sum,1);
for ns=0:n_study-1;
T_m2rz_{1+ns} = T_{1+ns}(:,1:2) - repmat(T_m2rz_med_,n_rij_cap_s_(1+ns),1);
end;%for ns=0:n_study-1;
n_z=4;
for nz=1:n_z;
T_m2rx_full_ = zeros(size(T_m2rz_full_,1),size(T_m2rz_full_,2)*nz);
T_m2rx_full_(:,1:2) = T_m2rz_full_;
for ns=0:n_study-1;
T_m2rx_{1+ns} = zeros(size(T_m2rz_{1+ns},1),size(T_m2rz_{1+ns},2)*nz);
T_m2rx_{1+ns}(:,1:2) = T_m2rz_{1+ns};
end;%for ns=0:n_study-1;
w_tmp = pi/2/nz; R_tmp = [cos(w_tmp),+sin(w_tmp);-sin(w_tmp),cos(w_tmp)];
for nx=2:nz;
T_m2rx_full_(:,(1:2) + (nx-1)*2) = T_m2rx_full_(:,(1:2) + (nx-2)*2)*R_tmp;
for ns=0:n_study-1;
T_m2rx_{1+ns}(:,(1:2) + (nx-1)*2) = T_m2rx_{1+ns}(:,(1:2) + (nx-2)*2)*R_tmp;
end;%for ns=0:n_study-1;
end;%for nx=2:nz;
mc_T_m2rx = ones(1 + 2*nz,1);
tmpchar = sprintf('%s/%smc_T_m2r%d.char',dir_out,str_output_prefix_local,nz); if (flag_verbose>0); disp(sprintf(' %% Writing %s of size %d sum %d',tmpchar,numel(mc_T_m2rx),sum(mc_T_m2rx>0))); end; 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(mc_T_m2rx(:)>0,tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
tmpchar = sprintf('%s/%sT_m2r%d_full_n.char',dir_out,str_output_prefix_local,nz); if (flag_verbose>0); disp(sprintf(' %% Writing %s of size %d-x%d',tmpchar,size(T_m2rx_full_,1),1+size(T_m2rx_full_,2))); end; 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress([ones(nrows_T_sum,1) , (T_m2rx_full_>0)],tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
tmpchar = sprintf('%s/%sT_m2r%d_full_t.char',dir_out,str_output_prefix_local,nz); if (flag_verbose>0); disp(sprintf(' %% Writing %s of size %d-x%d',tmpchar,1+size(T_m2rx_full_,2),size(T_m2rx_full_,1))); end; 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(transpose([ones(nrows_T_sum,1) , (T_m2rx_full_>0)]),tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
save(sprintf('%s/%sT_m2r%d_full_n.mat',dir_out,str_output_prefix_local,nz),'T_m2rz_med_','T_m2rz_full_','T_m2rx_full_','w_tmp','R_tmp');
for ns=0:n_study-1;
tmpchar = sprintf('%s/%sT_m2r%d_%.2d_n.char',dir_out,str_output_prefix_local,nz,1+ns); if (flag_verbose>0); disp(sprintf(' %% Writing %s of size %d-x%d',tmpchar,size(T_m2rx_{1+ns},1),1+size(T_m2rx_{1+ns},2))); end; 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress([ones(n_rij_cap_s_(1+ns),1) , (T_m2rx_{1+ns}>0)],tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
tmpchar = sprintf('%s/%sT_m2r%d_%.2d_t.char',dir_out,str_output_prefix_local,nz,1+ns); if (flag_verbose>0); disp(sprintf(' %% Writing %s of size %d-x%d',tmpchar,1+size(T_m2rx_{1+ns},2),size(T_m2rx_{1+ns},1))); end; 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(transpose([ones(n_rij_cap_s_(1+ns),1) , (T_m2rx_{1+ns}>0)]),tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
end;%for ns=0:n_study-1;
save(sprintf('%s/%sT_m2r%d_xx_n.mat',dir_out,str_output_prefix_local,nz),'T_m2rz_med_','T_m2rz_','T_m2rx_','w_tmp','R_tmp');
tmp_M = size(T_m2rz_full_,1); tmp_NT = size(T_m2rx_full_,2);
tmp_b_ = 2*(T_m2rx_full_>0)-1;
tmp_d = sum(sum(tmp_b_(:,:),1).^2)/(tmp_M*tmp_M*tmp_NT);
tmp_c_ = zeros(tmp_M,1);
for nr=0:tmp_M-1;
u_ = T_m2rz_full_(1+nr,:); u_ = u_/norm(u_,'fro');
tmp_ij = find((T_m2rz_full_*transpose(u_))>0); tmp_m = numel(tmp_ij);
tmp_c_(1+nr) = sum(sum(tmp_b_(tmp_ij,:),1).^2)/(tmp_m*tmp_m*tmp_NT);
end;%for nr=0:tmp_M-1;
tmp_kappa_squared = mean(tmp_c_);
tmpchar = sprintf('%s/%sT_m2r%d_kappa_squared.txt',dir_out,str_output_prefix_local,nz);
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
tmp_fp = fopen(tmpchar,'w'); fprintf(tmp_fp,'%0.16f\n',tmp_kappa_squared); fclose(tmp_fp);
if (flag_verbose>0); disp(sprintf(' %% nx %d: kappa^2 = %0.6f (versus %0.6f) --> writing to %s',nz,tmp_kappa_squared,tmp_d,tmpchar)); end;
end;%if (flag_force_create | ~exist(tmpchar,'file'));
end;%for nz=1:n_z;
%%%%%%%%;
end;%if (size(T_full_,2)>=2);

if flag_stop>3;
disp(' %% returning after mds'); return; 
end;%if flag_stop>3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% defining ncols_An, nrows_An_s_ and nrows_An_full_ ;')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ncols_An = 3*n_snp_cap_sub; %<-- This is for the .b16 file, which has 3 columns per snp. ;
ncols_An = 1*n_snp_cap_sub; %<-- Now we only keep 1 column per snp, corresponding to the dosage. ;
%ncols_An_extend = mod(bitj - mod(ncols_An,bitj),bitj); %<-- not needed. ;
%lcols_An = (ncols_An + ncols_An_extend)/bit8; %<-- not needed. ;
nrows_An_s_ = zeros(1,n_study);
for ns=0:n_study-1;
nrows_An_s_(1+ns) = n_rij_cap_s_(1+ns);
%nrows_An_extend_s_(1+ns) = mod(bitj - mod(nrows_An_s_(1+ns),bitj),bitj); %<-- not needed. ;
%lrows_An_s_(1+ns) = (nrows_An_s_(1+ns) + nrows_An_extend_s_(1+ns))/bit8; %<-- not needed. ;
if (flag_verbose>0); disp(sprintf(' %% nrows_An_s_(1+%d): %d',ns,nrows_An_s_(1+ns))); end;
end;%for ns=0:n_study-1;
nrows_An_full = sum(nrows_An_s_);
%nrows_An_full_extend = mod(bitj - mod(nrows_An_full,bitj),bitj); %<-- not needed. ;
%lrows_An_full = (nrows_An_full + nrows_An_full_extend)/bit8; %<-- not needed. ;
if (flag_verbose>0); disp(sprintf(' %% nrows_An_full: %d',nrows_An_full)); end;
nrows_An_csum_s_ = cumsum([0,nrows_An_s_]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Writing extended bim and fam files (i.e., char_bim.ext and char_fam.ext). ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Bim: ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% This has the same structure as a typical bim file, ; ')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Chromosome code (either an integer, or ''X''/''Y''/''XY''/''MT''; ''0'' indicates unknown) or name ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Variant identifier ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Position in morgans or centimorgans (safe to use dummy value of ''0'') ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Base-pair coordinate (1-based; limited to 2^31-2) ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Allele 1 (corresponding to clear bits in .bed; usually minor) ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Allele 2 (corresponding to set bits in .bed; usually major) ;')); end;
if (flag_verbose>0); disp(sprintf(' %% with five extra fields. ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% field 6+1: allele type. ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% field 6+2: total relative entropy across (nonexcluded) patients. ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% field 6+3: frequency of allele type (i.e., sparsity of column). ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% field 6+4: frequency of missing data for that snp. ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% field 6+5: minor-allele-frequency of that snp. ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Moreover, the allele-ordering should be sorted in terms of column-sparsity. ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_bim_ext;
tmp_bim = cell(6 + 5);
tmp_bim{ 1} = zeros(ncols_An,1);
tmp_bim{ 2} = cell(ncols_An,1);
tmp_bim{ 3} = zeros(ncols_An,1);
tmp_bim{ 4} = zeros(ncols_An,1);
tmp_bim{ 5} = cell(ncols_An,1);
tmp_bim{ 6} = cell(ncols_An,1);
tmp_bim{ 7} = cell(ncols_An,1);
tmp_bim{ 8} = zeros(ncols_An,1);
tmp_bim{ 9} = zeros(ncols_An,1);
tmp_bim{10} = zeros(ncols_An,1);
tmp_bim{11} = zeros(ncols_An,1);
for nl3=0:n_snp_cap_sub-1; 
nl2 = ij_snp_cap_sub_(1+nl3)-1;
ns=0; nla = ij_snp_cap_sub_sl__{1+ns}(1+nl3)-1; nlb = ij_snp_cap_sl__{1+ns}(1+nl2)-1; assert(nla==nlb); nl1=nla;
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_bim{ 1}(cij_and_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+0}(1+nl1);
tmp_bim{ 2}(cij_and_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+1}(1+nl1);
tmp_bim{ 3}(cij_and_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+2}(1+nl1);
tmp_bim{ 4}(cij_and_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+3}(1+nl1);
tmp_bim{ 5}(cij_and_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+4}(1+nl1);
tmp_bim{ 6}(cij_and_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+5}(1+nl1);
tmp_bim{ 7}(cij_and_rpl_l_(1+nl3)) = {'and(AA)'};
tmp_bim{ 8}(cij_and_rpl_l_(1+nl3)) = I_opt_tot_l_(1+nl2);
tmp_bim{ 9}(cij_and_rpl_l_(1+nl3)) = fr2_and_tot_l_(1+nl2);
tmp_bim{10}(cij_and_rpl_l_(1+nl3)) = frq_mss_tot_l_(1+nl2);
tmp_bim{11}(cij_and_rpl_l_(1+nl3)) = min(p_opt_tot_l_(1+nl2),q_opt_tot_l_(1+nl2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_bim{ 1}(cij_xor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+0}(1+nl1);
tmp_bim{ 2}(cij_xor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+1}(1+nl1);
tmp_bim{ 3}(cij_xor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+2}(1+nl1);
tmp_bim{ 4}(cij_xor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+3}(1+nl1);
tmp_bim{ 5}(cij_xor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+4}(1+nl1);
tmp_bim{ 6}(cij_xor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+5}(1+nl1);
tmp_bim{ 7}(cij_xor_rpl_l_(1+nl3)) = {'xor(AG)'};
tmp_bim{ 8}(cij_xor_rpl_l_(1+nl3)) = I_opt_tot_l_(1+nl2);
%tmp_bim{ 9}(cij_xor_rpl_l_(1+nl3)) = fr2_and_tot_l_(1+nl2) + fr2_xor_tot_l_(1+nl2);
tmp_bim{ 9}(cij_xor_rpl_l_(1+nl3)) = fr2_xor_tot_l_(1+nl2);
tmp_bim{10}(cij_xor_rpl_l_(1+nl3)) = frq_mss_tot_l_(1+nl2);
tmp_bim{11}(cij_xor_rpl_l_(1+nl3)) = min(p_opt_tot_l_(1+nl2),q_opt_tot_l_(1+nl2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_bim{ 1}(cij_nor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+0}(1+nl1);
tmp_bim{ 2}(cij_nor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+1}(1+nl1);
tmp_bim{ 3}(cij_nor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+2}(1+nl1);
tmp_bim{ 4}(cij_nor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+3}(1+nl1);
tmp_bim{ 5}(cij_nor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+4}(1+nl1);
tmp_bim{ 6}(cij_nor_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+5}(1+nl1);
tmp_bim{ 7}(cij_nor_rpl_l_(1+nl3)) = {'nor(GG)'};
tmp_bim{ 8}(cij_nor_rpl_l_(1+nl3)) = I_opt_tot_l_(1+nl2);
tmp_bim{ 9}(cij_nor_rpl_l_(1+nl3)) = fr2_nor_tot_l_(1+nl2);
tmp_bim{10}(cij_nor_rpl_l_(1+nl3)) = frq_mss_tot_l_(1+nl2);
tmp_bim{11}(cij_nor_rpl_l_(1+nl3)) = min(p_opt_tot_l_(1+nl2),q_opt_tot_l_(1+nl2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_bim{ 1}(cij_sum_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+0}(1+nl1);
tmp_bim{ 2}(cij_sum_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+1}(1+nl1);
tmp_bim{ 3}(cij_sum_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+2}(1+nl1);
tmp_bim{ 4}(cij_sum_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+3}(1+nl1);
tmp_bim{ 5}(cij_sum_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+4}(1+nl1);
tmp_bim{ 6}(cij_sum_rpl_l_(1+nl3)) = bim_s_{1+ns}{1+5}(1+nl1);
tmp_bim{ 7}(cij_sum_rpl_l_(1+nl3)) = {'sum(0*GG+1*AG+2*AA)'};
tmp_bim{ 8}(cij_sum_rpl_l_(1+nl3)) = I_opt_tot_l_(1+nl2);
tmp_bim{ 9}(cij_sum_rpl_l_(1+nl3)) = fr2_sum_tot_l_(1+nl2);
tmp_bim{10}(cij_sum_rpl_l_(1+nl3)) = frq_mss_tot_l_(1+nl2);
tmp_bim{11}(cij_sum_rpl_l_(1+nl3)) = min(p_opt_tot_l_(1+nl2),q_opt_tot_l_(1+nl2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nl3=0:n_snp_cap_sub-1;
tmpchar = sprintf('%s/%schar_bim.ext',dir_out,str_output_prefix_local);
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
fid = fopen(tmpchar,'w');
for nc=0:ncols_An-1;
fprintf(fid ...
,'%d\t%s\t%d\t%d\t%s\t%s\t%s\t%f\t%f\t%f\t%f\n' ...
,tmp_bim{ 1}(1+nc) ...
,tmp_bim{ 2}{1+nc} ...
,tmp_bim{ 3}(1+nc) ...
,tmp_bim{ 4}(1+nc) ...
,tmp_bim{ 5}{1+nc} ...
,tmp_bim{ 6}{1+nc} ...
,tmp_bim{ 7}{1+nc} ...
,tmp_bim{ 8}(1+nc) ...
,tmp_bim{ 9}(1+nc) ...
,tmp_bim{10}(1+nc) ...
,tmp_bim{11}(1+nc) ...
);
end;%for nc=0:ncols_An-1;
fclose(fid);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
end;% if flag_bim_ext;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Fam: ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% This has the same structure as a typical fam file, ; ')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Family ID (''FID'') ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Within-family ID (''IID''; cannot be ''0'') ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Within-family ID of father (''0'' if father isn''t in dataset) ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Within-family ID of mother (''0'' if mother isn''t in dataset) ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Sex code (''1'' = male, ''2'' = female, ''0'' = unknown) ;')); end;
if (flag_verbose>0); disp(sprintf(' %%\t Phenotype value (''1'' = control, ''2'' = case, ''-9''/''0''/non-numeric = missing data if case/control) ;')); end;
if (flag_verbose>0); disp(sprintf(' %% with one extra field. ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% field 6+1: directory and filename from which fam file was drawn . ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_fam_ext;
tmp_fam = cell(6 + 1);
tmp_fam{1+0} = cell(nrows_An_full,1);
tmp_fam{1+1} = cell(nrows_An_full,1);
tmp_fam{1+2} = cell(nrows_An_full,1);
tmp_fam{1+3} = cell(nrows_An_full,1);
tmp_fam{1+4} = zeros(nrows_An_full,1);
tmp_fam{1+5} = zeros(nrows_An_full,1);
tmp_fam{1+6} = cell(nrows_An_full,1);
nr=0;
for ns=0:n_study-1;
for np2=0:nrows_An_s_(1+ns)-1;
np1 = rindex_cap_sp__{1+ns}(1+np2);
tmp_fam{1+0}(1+nr) = fam_s_{1+ns}{1+0}(1+np1);
tmp_fam{1+1}(1+nr) = fam_s_{1+ns}{1+1}(1+np1);
tmp_fam{1+2}(1+nr) = fam_s_{1+ns}{1+2}(1+np1);
tmp_fam{1+3}(1+nr) = fam_s_{1+ns}{1+3}(1+np1);
tmp_fam{1+4}(1+nr) = fam_s_{1+ns}{1+4}(1+np1);
tmp_fam{1+5}(1+nr) = fam_s_{1+ns}{1+5}(1+np1);
tmp_fam{1+6}(1+nr) = {sprintf('%s/%s/%s.fam',dir_trunk,study_name_of_branch_s_{1+ns},study_name_without_extension_s_{1+ns})};
nr = nr+1;
end;%for np2=0:nrows_An_s_(1+ns)-1;
end;%for ns=0:n_study-1;
tmpchar = sprintf('%s/%schar_fam.ext',dir_out,str_output_prefix_local);
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
fid = fopen(tmpchar,'w');
for nr=0:nrows_An_full-1;
fprintf( ...
 fid ...
,'%s\t%s\t%s\t%s\t%d\t%d\t%s\n' ...
,tmp_fam{1+0}{1+nr} ...
,tmp_fam{1+1}{1+nr} ...
,tmp_fam{1+2}{1+nr} ...
,tmp_fam{1+3}{1+nr} ...
,tmp_fam{1+4}(1+nr) ...
,tmp_fam{1+5}(1+nr) ...
,tmp_fam{1+6}{1+nr} ...
);
end;%for nr=0:nrows_An_full-1;
fclose(fid);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
end;%if flag_fam_ext;

if flag_stop>2;
disp(' %% returning after char_bim.ext and char_fam.ext'); return; 
end;%if flag_stop>2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Writing mc and mr masks ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mc_A = ones(ncols_An,1);
tmpchar = sprintf('%s/%smc_A.char',dir_out,str_output_prefix_local); 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(mc_A(:)>0,tmpchar); 
end;%if (flag_force_create | ~exist(tmpchar,'file'));
mr_A_full = zeros(nrows_An_full,1); mr_Z_full = zeros(nrows_An_full,1); 
for ns=0:n_study-1;
mr_A_{1+ns} = zeros(nrows_An_s_(1+ns),1); mr_Z_{1+ns} = zeros(nrows_An_s_(1+ns),1);
tmp_ij = find(fam_s_{1+ns}{1+5}(1+rindex_cap_sp__{1+ns})==case_tag);
mr_A_{1+ns}(tmp_ij) = 1; mr_A_full(nrows_An_csum_s_(1+ns) + tmp_ij) = 1;
tmp_ij = find(fam_s_{1+ns}{1+5}(1+rindex_cap_sp__{1+ns})==ctrl_tag);
mr_Z_{1+ns}(tmp_ij) = 1; mr_Z_full(nrows_An_csum_s_(1+ns) + tmp_ij) = 1;
end;%for ns=0:n_study-1;
tmpchar = sprintf('%s/%smr_A_full.char',dir_out,str_output_prefix_local);if (flag_verbose>0); disp(sprintf(' %% Writing %s of size %d sum %d',tmpchar,numel(mr_A_full),sum(mr_A_full>0))); end; 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(mr_A_full(:)>0,tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
pause(1);
tmpchar = sprintf('%s/%smr_Z_full.char',dir_out,str_output_prefix_local);if (flag_verbose>0); disp(sprintf(' %% Writing %s of size %d sum %d',tmpchar,numel(mr_Z_full),sum(mr_Z_full>0))); end; 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(mr_Z_full(:)>0,tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
pause(1);
for ns=0:n_study-1;
tmpchar = sprintf('%s/%smr_A_%.2d.char',dir_out,str_output_prefix_local,1+ns);if (flag_verbose>0); disp(sprintf(' %% Writing %s of size %d sum %d',tmpchar,numel(mr_A_{1+ns}),sum(mr_A_{1+ns}>0))); end; 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(mr_A_{1+ns}(:)>0,tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
pause(1);
tmpchar = sprintf('%s/%smr_Z_%.2d.char',dir_out,str_output_prefix_local,1+ns);if (flag_verbose>0); disp(sprintf(' %% Writing %s of size %d sum %d',tmpchar,numel(mr_Z_{1+ns}),sum(mr_Z_{1+ns}>0))); end; 
if (flag_force_create | ~exist(tmpchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpchar)); end;
char_compress(mr_Z_{1+ns}(:)>0,tmpchar);
end;%if (flag_force_create | ~exist(tmpchar,'file'));
pause(1);
end;%for ns=0:n_study-1;

if flag_stop>1;
disp(' %% returning after masks'); return; 
end;%if flag_stop>1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Initializing An_ and At_ ;')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear A_n_ A_n_full A_t_ A_t_full;
An_full = zeros(nrows_An_full,ncols_An,'uint8');
At_full = zeros(ncols_An,nrows_An_full,'uint8');
for ns=0:n_study-1;
An_{1+ns} = zeros(nrows_An_s_(1+ns),ncols_An,'uint8');
At_{1+ns} = zeros(ncols_An,nrows_An_s_(1+ns),'uint8');
end;%for ns=0:n_study-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% Writing An_ and At_ ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ns=0:n_study-1;
if (flag_verbose>0); disp(sprintf(' %% ns %d: %s/%s',ns,study_name_of_branch_s_{1+ns},study_name_without_extension_s_{1+ns})); end;
tmp_ij_row_ = 1 + nrows_An_csum_s_(1+ns) + (0:nrows_An_s_(1+ns)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reading bed file (binarized snp data) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_bed_s_{1+ns} = fopen(fn_bed_s_{1+ns}); 
key1 = fread(fp_bed_s_{1+ns},1,'uint8'); key2 = fread(fp_bed_s_{1+ns},1,'uint8'); key3 = fread(fp_bed_s_{1+ns},1,'uint8'); assert(key1==108); assert(key2==27); assert(key3==1);
flag_ver=0;
if flag_ver==0;
if (flag_verbose>0); disp(sprintf(' %% version 0: simultaneous read (assuming sufficient ram)... ')); end;tic;
uchar_11111111 = transpose(fread(fp_bed_s_{1+ns},[n_patient_bed_s_(1+ns),n_snp_s_(1+ns)],'uint8=>uint8')); fclose(fp_bed_s_{1+ns});
t_tmp = toc; if (flag_verbose>0); disp(sprintf(' %% read input: time %0.2f seconds',t_tmp)); end; 
if (flag_verbose>0); disp(sprintf(' %% reading 4 allele pairs per byte... ')); end; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% these commands use mod, which is very slow. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%uchar_00000011 = mod(uchar_11111111,  4);
%uchar_00001111 = mod(uchar_11111111, 16); uchar_00001100 = idivide(uchar_00001111 - uchar_00000011, 4) ; 
%uchar_00111111 = mod(uchar_11111111, 64); uchar_00110000 = idivide(uchar_00111111 - uchar_00001111,16) ; clear uchar_00001111;
%                                          uchar_11000000 = idivide(uchar_11111111 - uchar_00111111,64) ; clear uchar_00111111;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% these commands use bitshifts, which are much faster. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
uchar_00000011 = bitshift(bitshift(uchar_11111111,+6),-6);
uchar_00001100 = bitshift(bitshift(bitshift(uchar_11111111,-2),+6),-6);
uchar_00110000 = bitshift(bitshift(bitshift(uchar_11111111,-4),+6),-6);
uchar_11000000 = bitshift(uchar_11111111,-6);
clear uchar_11111111;
stride1 = 1:4:4*n_patient_bed_s_(1+ns); [~,tmp_ij] = setdiff(stride1,1+rindex_keep_s_{1+ns},'stable'); uchar_00000011(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride2 = 2:4:4*n_patient_bed_s_(1+ns); [~,tmp_ij] = setdiff(stride2,1+rindex_keep_s_{1+ns},'stable'); uchar_00001100(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride3 = 3:4:4*n_patient_bed_s_(1+ns); [~,tmp_ij] = setdiff(stride3,1+rindex_keep_s_{1+ns},'stable'); uchar_00110000(:,tmp_ij)=snp_not_tag; clear tmp_ij;
stride4 = 4:4:4*n_patient_bed_s_(1+ns); [~,tmp_ij] = setdiff(stride4,1+rindex_keep_s_{1+ns},'stable'); uchar_11000000(:,tmp_ij)=snp_not_tag; clear tmp_ij;
if (flag_verbose>0); disp(sprintf(' %% writing binary files... ')); end; t_start = tic; t_sub = tic;
%br_u = zeros(1,bit8,'uint8'); br_u = cast(2.^(bit8-1:-1:0),'uint8'); %<-- not necessary. ;
%br_d = zeros(1,bit8,'double'); br_d = cast(2.^(bit8-1:-1:0),'double'); %<-- not necessary. ;
for nl3=0:n_snp_cap_sub-1; 
% nl3 refers to the index within 0:n_snp_cap_sub-1. ;
% cl3_and refers to the column index for the and-type (within b16 file). ;
% cl3_xor refers to the column index for the xor-type (within b16 file). ;
% cl3_nor refers to the column index for the nor-type (within b16 file). ;
%bl3_and = mod((cij_and_rpl_l_(1+nl3)-1),bit8); al3_and = floor((cij_and_rpl_l_(1+nl3)-1)/bit8); cl3_and = (cij_and_rpl_l_(1+nl3)-1);
%bl3_xor = mod((cij_xor_rpl_l_(1+nl3)-1),bit8); al3_xor = floor((cij_xor_rpl_l_(1+nl3)-1)/bit8); cl3_xor = (cij_xor_rpl_l_(1+nl3)-1);
%bl3_nor = mod((cij_nor_rpl_l_(1+nl3)-1),bit8); al3_nor = floor((cij_nor_rpl_l_(1+nl3)-1)/bit8); cl3_nor = (cij_nor_rpl_l_(1+nl3)-1);
cl3_sum = (cij_sum_rpl_l_(1+nl3)-1);
% nl2 refers to the index listing from snp_cap_ ;
nl2 = ij_snp_cap_sub_(1+nl3)-1;
% nla and/or nlb refers to the original index within the bed file. ;
nla = ij_snp_cap_sub_sl__{1+ns}(1+nl3)-1; nlb = ij_snp_cap_sl__{1+ns}(1+nl2)-1; assert(nla==nlb); nl1=nla;
if (mod(nl3,1024)==0); 
t_tmp = toc(t_sub); t_tot_tmp = t_tmp/(1+nl3)*n_snp_cap_sub; 
if (flag_verbose>0); disp(sprintf(' %% nl3 %d/%d nl2 %d/%d nl1 %d/%d; time %0.2f seconds; estimated total time %0.2f seconds (%0.2f hours)',nl3,n_snp_cap_sub,nl2,n_snp_cap,nl1,n_snp_s_(1+ns),t_tmp,t_tot_tmp,t_tot_tmp/3600)); end; 
end;%if (mod(nl3,1024)==0); 
% data_snp_p_ holds the tags for the patients (columns) associated with snp nl1. ;
data_snp_p_ = zeros(1,4*n_patient_bed_s_(1+ns));
data_snp_p_(stride1) = uchar_00000011(1+nl1,:);
data_snp_p_(stride2) = uchar_00001100(1+nl1,:);
data_snp_p_(stride3) = uchar_00110000(1+nl1,:);
data_snp_p_(stride4) = uchar_11000000(1+nl1,:);
% truncate to exactly n_patient_fam_s_(1+ns) entries. ;
data_snp_p_ = data_snp_p_(1,1:n_patient_fam_s_(1+ns));
% We replace the missing data with the mode below, stored in fr2_mode_l_(1+nl2). ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if (snp_sg1_ls__(1+nl2,1+ns)==+1);
data_snp_and_p_ = ( (data_snp_p_==snp_and_tag) | ( (fr2_mode_l_(1+nl2)==snp_and_tag) & (data_snp_p_==snp_mss_tag) ) ) ; 
end;%if (snp_sg1_ls__(1+nl2,1+ns)==+1);
if (snp_sg1_ls__(1+nl2,1+ns)==-1);
data_snp_and_p_ = ( (data_snp_p_==snp_nor_tag) | ( (fr2_mode_l_(1+nl2)==snp_and_tag) & (data_snp_p_==snp_mss_tag) ) ) ; 
end;%if (snp_sg1_ls__(1+nl2,1+ns)==-1);
data_snp_and_sub_p_ = data_snp_and_p_(1,1+rindex_cap_sp__{1+ns}); data_snp_and_bin_u_p_ = cast(data_snp_and_sub_p_(1,:),'uint8');
data_snp_and_sub_full_p_ = zeros(1,nrows_An_full); data_snp_and_sub_full_p_(1,tmp_ij_row_) = data_snp_and_sub_p_;
%data_snp_and_bin_d_p_ = cast(data_snp_and_sub_p_(1,:),'double'); 
%data_snp_and_full_bin_d_p_ = cast(data_snp_and_sub_full_p_(1,:),'double'); 
%An_full(1:lrows_An_full,1+cl3_and) = An_full(1:lrows_An_full,1+cl3_and) + cast(transpose(br_d*reshape([data_snp_and_full_bin_d_p_,zeros(1,nrows_An_full_extend)],bit8,lrows_An_full)),'uint8');
%At_full(1+al3_and,tmp_ij_row_) = At_full(1+al3_and,tmp_ij_row_) + br_u(1+bl3_and)*data_snp_and_bin_u_p_;
%An_{1+ns}(1:lrows_An_s_(1+ns),1+cl3_and) = An_{1+ns}(1:lrows_An_s_(1+ns),1+cl3_and) + cast(transpose(br_d*reshape([data_snp_and_bin_d_p_,zeros(1,nrows_An_extend_s_(1+ns))],bit8,lrows_An_s_(1+ns))),'uint8');
%At_{1+ns}(1+al3_and,(1:nrows_An_s_(1+ns))) = At_{1+ns}(1+al3_and,(1:nrows_An_s_(1+ns))) + br_u(1+bl3_and)*data_snp_and_bin_u_p_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
data_snp_xor_p_ = ( (data_snp_p_==snp_xor_tag) | ( (fr2_mode_l_(1+nl2)==snp_xor_tag) & (data_snp_p_==snp_mss_tag) ) ) ;
data_snp_xor_sub_p_ = data_snp_xor_p_(1,1+rindex_cap_sp__{1+ns}); data_snp_xor_bin_u_p_ = cast(data_snp_xor_sub_p_(1,:),'uint8');
data_snp_xor_sub_full_p_ = zeros(1,nrows_An_full); data_snp_xor_sub_full_p_(1,tmp_ij_row_) = data_snp_xor_sub_p_;
%data_snp_xor_bin_d_p_ = cast(data_snp_xor_sub_p_(1,:),'double'); 
%data_snp_xor_full_bin_d_p_ = cast(data_snp_xor_sub_full_p_(1,:),'double'); 
%An_full(1:lrows_An_full,1+cl3_xor) = An_full(1:lrows_An_full,1+cl3_xor) + cast(transpose(br_d*reshape([data_snp_xor_full_bin_d_p_,zeros(1,nrows_An_full_extend)],bit8,lrows_An_full)),'uint8');
%At_full(1+al3_xor,tmp_ij_row_) = At_full(1+al3_xor,tmp_ij_row_) + br_u(1+bl3_xor)*data_snp_xor_bin_u_p_;
%An_{1+ns}(1:lrows_An_s_(1+ns),1+cl3_xor) = An_{1+ns}(1:lrows_An_s_(1+ns),1+cl3_xor) + cast(transpose(br_d*reshape([data_snp_xor_bin_d_p_,zeros(1,nrows_An_extend_s_(1+ns))],bit8,lrows_An_s_(1+ns))),'uint8');
%At_{1+ns}(1+al3_xor,(1:nrows_An_s_(1+ns))) = At_{1+ns}(1+al3_xor,(1:nrows_An_s_(1+ns))) + br_u(1+bl3_xor)*data_snp_xor_bin_u_p_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if (snp_sg1_ls__(1+nl2,1+ns)==+1);
data_snp_nor_p_ = ( (data_snp_p_==snp_nor_tag) | ( (fr2_mode_l_(1+nl2)==snp_nor_tag) & (data_snp_p_==snp_mss_tag) ) ) ; 
end;%if (snp_sg1_ls__(1+nl2,1+ns)==+1);
if (snp_sg1_ls__(1+nl2,1+ns)==-1);
data_snp_nor_p_ = ( (data_snp_p_==snp_and_tag) | ( (fr2_mode_l_(1+nl2)==snp_nor_tag) & (data_snp_p_==snp_mss_tag) ) ) ; 
end;%if (snp_sg1_ls__(1+nl2,1+ns)==-1);
data_snp_nor_sub_p_ = data_snp_nor_p_(1,1+rindex_cap_sp__{1+ns}); data_snp_nor_bin_u_p_ = cast(data_snp_nor_sub_p_(1,:),'uint8');
data_snp_nor_sub_full_p_ = zeros(1,nrows_An_full); data_snp_nor_sub_full_p_(1,tmp_ij_row_) = data_snp_nor_sub_p_;
%data_snp_nor_bin_d_p_ = cast(data_snp_nor_sub_p_(1,:),'double'); 
%data_snp_nor_full_bin_d_p_ = cast(data_snp_nor_sub_full_p_(1,:),'double'); 
%An_full(1:lrows_An_full,1+cl3_nor) = An_full(1:lrows_An_full,1+cl3_nor) + cast(transpose(br_d*reshape([data_snp_nor_full_bin_d_p_,zeros(1,nrows_An_full_extend)],bit8,lrows_An_full)),'uint8');
%At_full(1+al3_nor,tmp_ij_row_) = At_full(1+al3_nor,tmp_ij_row_) + br_u(1+bl3_nor)*data_snp_nor_bin_u_p_;
%An_{1+ns}(1:lrows_An_s_(1+ns),1+cl3_nor) = An_{1+ns}(1:lrows_An_s_(1+ns),1+cl3_nor) + cast(transpose(br_d*reshape([data_snp_nor_bin_d_p_,zeros(1,nrows_An_extend_s_(1+ns))],bit8,lrows_An_s_(1+ns))),'uint8');
%At_{1+ns}(1+al3_nor,(1:nrows_An_s_(1+ns))) = At_{1+ns}(1+al3_nor,(1:nrows_An_s_(1+ns))) + br_u(1+bl3_nor)*data_snp_nor_bin_u_p_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
An_full(:,1+cl3_sum) = An_full(:,1+cl3_sum) + reshape(cast(0*data_snp_nor_sub_full_p_ + 1*data_snp_xor_sub_full_p_ + 2*data_snp_and_sub_full_p_,'uint8'),[nrows_An_full,1]);
At_full(1+cl3_sum,:) = At_full(1+cl3_sum,:) + reshape(cast(0*data_snp_nor_sub_full_p_ + 1*data_snp_xor_sub_full_p_ + 2*data_snp_and_sub_full_p_,'uint8'),[1,nrows_An_full]);
An_{1+ns}(:,1+cl3_sum) = An_{1+ns}(:,1+cl3_sum) + reshape(cast(0*data_snp_nor_sub_p_ + 1*data_snp_xor_sub_p_ + 2*data_snp_and_sub_p_,'uint8'),[nrows_An_s_(1+ns)],1);
At_{1+ns}(1+cl3_sum,:) = At_{1+ns}(1+cl3_sum,:) + reshape(cast(0*data_snp_nor_sub_p_ + 1*data_snp_xor_sub_p_ + 2*data_snp_and_sub_p_,'uint8'),[1,nrows_An_s_(1+ns)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
end;%for nl3=0:n_snp_cap_sub-1;
t_tmp = toc(t_start); if (flag_verbose>0); disp(sprintf(' %% bin: time %0.2f seconds',t_tmp)); end; 
clear uchar_00000011;clear uchar_00001100;clear uchar_00110000;clear uchar_11000000;
clear stride1 stride2 stride3 stride4;
end;%if flag_ver==0;
end;%for ns=0:n_study-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% ; ')); end;
if (flag_verbose>0); disp(sprintf(' %% writing A_n.char and A_t.char ; ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
tmpAnchar = sprintf('%s/%sA_full_n.char',dir_out,str_output_prefix_local);
if (flag_force_create | ~exist(tmpAnchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpAnchar)); end;
fid = fopen(tmpAnchar,'w');
%fwrite(fid,bitj,'int');
fwrite(fid,nrows_An_full,'int');
fwrite(fid,ncols_An,'int');
fwrite(fid,An_full(:),'uint8');
fclose(fid);
end;%if (flag_force_create | ~exist(tmpAnchar,'file'));
tmpAtchar = sprintf('%s/%sA_full_t.char',dir_out,str_output_prefix_local);
if (flag_force_create | ~exist(tmpAtchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpAtchar)); end;
fid = fopen(tmpAtchar,'w');
%fwrite(fid,bitj,'int');
fwrite(fid,ncols_An,'int');
fwrite(fid,nrows_An_full,'int');
fwrite(fid,At_full(:),'uint8');
fclose(fid);
end;%if (flag_force_create | ~exist(tmpAtchar,'file'));
for ns=0:n_study-1;
tmpAnchar = sprintf('%s/%sA_%.2d_n.char',dir_out,str_output_prefix_local,1+ns);
if (flag_force_create | ~exist(tmpAnchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpAnchar)); end;
fid = fopen(tmpAnchar,'w');
%fwrite(fid,bitj,'int');
fwrite(fid,nrows_An_s_(1+ns),'int');
fwrite(fid,ncols_An,'int');
fwrite(fid,An_{1+ns}(:),'uint8');
fclose(fid);
end;%if (flag_force_create | ~exist(tmpAnchar,'file'));
tmpAtchar = sprintf('%s/%sA_%.2d_t.char',dir_out,str_output_prefix_local,1+ns);
if (flag_force_create | ~exist(tmpAtchar,'file'));
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',tmpAtchar)); end;
fid = fopen(tmpAtchar,'w');
%fwrite(fid,bitj,'int');
fwrite(fid,ncols_An,'int');
fwrite(fid,nrows_An_s_(1+ns),'int');
fwrite(fid,At_{1+ns}(:),'uint8');
fclose(fid);
end;%if (flag_force_create | ~exist(tmpAtchar,'file'));
end;%for ns=0:n_study-1;
t_tmp = toc; if (flag_verbose>0); disp(sprintf(' %% char: time %0.2f seconds',t_tmp)); end; 

flag_check = flag_crop;
if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose>0); disp(sprintf(' %% checking .char files for consistency ;')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpAnchar = sprintf('%s/%sA_full_n.char',dir_out,str_output_prefix_local);
Bn_full = cast(char_uncompress(tmpAnchar,1:nrows_An_full,1:ncols_An),'double');
tmpAtchar = sprintf('%s/%sA_full_t.char',dir_out,str_output_prefix_local);
Bt_full = cast(char_uncompress(tmpAtchar,1:ncols_An,1:nrows_An_full),'double');
if (flag_verbose>0); disp(sprintf(' %% ns all: error %f',norm(Bn_full-transpose(Bt_full),'fro'))); end;
disp_flag=0;
if disp_flag;
subplot(1,2,1);imagesc(Bn_full,[-1,1]); subplot(1,2,2);imagesc(transpose(Bt_full),[-1,1]);
print('-djpeg',sprintf('%s/%sA_full_x.jpg',dir_out,str_output_prefix_local));
end;%if disp_flag;
for ns=0:n_study-1;
tmpAnchar = sprintf('%s/%sA_%.2d_n.char',dir_out,str_output_prefix_local,1+ns);
Bn_{1+ns} = cast(char_uncompress(tmpAnchar,1:nrows_An_s_(1+ns),1:ncols_An),'double');
tmpAtchar = sprintf('%s/%sA_%.2d_t.char',dir_out,str_output_prefix_local,1+ns);
Bt_{1+ns} = cast(char_uncompress(tmpAtchar,1:ncols_An,1:nrows_An_s_(1+ns)),'double');
if (flag_verbose>0); disp(sprintf(' %% ns %d: error %f',ns,norm(Bn_{1+ns}-transpose(Bt_{1+ns}),'fro'))); end;
disp_flag=0;
if disp_flag;
subplot(1,2,1);imagesc(Bn_{1+ns},[-1,1]); subplot(1,2,2);imagesc(transpose(Bt_{1+ns}),[-1,1]);
print('-djpeg',sprintf('%s/%sA_%.2d_x.jpg',dir_out,str_output_prefix_local,1+ns));
end;%if disp_flag;
end;%for ns=0:n_study-1;
end;%if flag_check;

if flag_stop>0;
disp(' %% returning after An&At');
if (flag_verbose>0); disp(sprintf(' %% [finished_bed_to_char_flip_ver7]')); end;
return; 
end;%if flag_stop>0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (flag_force_create | ~flag_exist_all);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;










