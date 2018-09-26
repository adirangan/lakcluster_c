function lakcluster_test_scramble_0(dir_trunk,N,verbose_flag)
% creating a larger example involving multiple categorical and continuous covariates, as well as cases and controls (but not involving Y). ;
% This only implants a single bicluster, and is simpler than lakcluster_test_uADZSZDA_ver8 (which implants multiple biclusters). ;
% uses lakcluster_ver18 ;
% try: ;
%{
  N=1024;verbose_flag=0;lakcluster_test_scramble_0(pwd,N);
  %}

if (nargin<3); verbose_flag=0; end;
if (nargin<2); N = 1024; end;
if (nargin<1); dir_trunk = pwd; end;
% clear; dir_trunk = pwd; N = 1024; verbose_flag=0;

nbins = 2; 
nrows_M_ = ceil(N/nbins) + 10*(1:nbins) ;
ncols_A = ceil(N/1920 + 2)*1920; pcols_A = ceil(ncols_A/1920); A_p = transpose(linspace(0.2,0.8,pcols_A));
ncols_Y = 0;
nmds = 0; ncols_T = 1;
QR_strategy = 'YnWt condense';
QC_strategy = 'YnWt store one';

test_string = 'lakcluster_test_scramble';
dir_base = sprintf('%s/dir_%s',dir_trunk,test_string); if ~exist(dir_base,'dir'); mkdir(dir_base); end;
in__name = sprintf('%s__in',test_string);
out_name = sprintf('%s_out',test_string);
dir__in = sprintf('%s/dir_%s',dir_base,in__name); if ~exist(dir__in,'dir'); mkdir(dir__in); end;
dir_out = sprintf('%s/dir_%s',dir_base,out_name); if ~exist(dir_out,'dir'); mkdir(dir_out); end; 
delete(sprintf('%s/An_00_scramble_pos.b16')); delete(sprintf('%s/At_00_scramble_pos.b16'));
delete(sprintf('%s/An_01_scramble_pos.b16')); delete(sprintf('%s/At_01_scramble_pos.b16'));
delete(sprintf('%s/An_00_scramble_pre.b16')); delete(sprintf('%s/At_00_scramble_pre.b16'));
delete(sprintf('%s/An_01_scramble_pre.b16')); delete(sprintf('%s/At_01_scramble_pre.b16'));

%path(path,dir__in); path(path,dir_out);


clear Mn Mt mr_A mr_Z mc_A bitj;
bitj = 16;

rng(1);
mc_A = (rand(1,ncols_A)>0.5);
for nb1=0:nbins-1; 
mr_A{1+nb1} = (rand(nrows_M_(1+nb1),1)>0.50);
mr_Z{1+nb1} = (rand(nrows_M_(1+nb1),1)>0.50); 
mr_Z{1+nb1} = mr_Z{1+nb1} - (mr_Z{1+nb1} & mr_A{1+nb1});
end;%for nb1=0:nbins-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate checkered matrices ;
for nb1=0:nbins-1; 
Mn{1+nb1} = -ones(nrows_M_(1+nb1),ncols_A);
Mn{1+nb1}(find(mr_A{1+nb1}),find(mc_A)) = 2*make_checkered_array(sum(mr_A{1+nb1}),sum(mc_A),17 + 50*nb1)-1;
Mn{1+nb1}(find(mr_Z{1+nb1}),find(mc_A)) = 2*make_checkered_array(sum(mr_Z{1+nb1}),sum(mc_A),31 + 75*nb1)-1;
Mt{1+nb1} = transpose(Mn{1+nb1});
end;%for nb1=0:nbins-1; 

flag_disp=0;
if flag_disp;
for nb1=0:nbins-1;
subplot(2,nbins,1+nb1+0*nbins); imagesc(Mn{1+nb1}(find(mr_A{1+nb1}),find(mc_A)));
subplot(2,nbins,1+nb1+1*nbins); imagesc(Mn{1+nb1}(find(mr_Z{1+nb1}),find(mc_A)));
end;%for nb1=0:nbins-1;
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpchar = sprintf('%s/A_p.mda',dir__in);mda_compress(A_p,tmpchar);
tmpAtchar = sprintf('%s/mc_A.b16',dir__in);tutorial_binary_compress(bitj,mc_A(:)>0,tmpAtchar);
for nb1=0:nbins-1;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Mn{1+nb1}>0,tmpAnchar); 
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Mt{1+nb1}>0,tmpAtchar); 
tmpAtchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_A{1+nb1}(:)>0,tmpAtchar);
tmpZtchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_Z{1+nb1}(:)>0,tmpZtchar);
end;%for nb1=0:nbins-1;

for nb1=0:nbins-1;
rij{1+nb1} = find(mr_A{1+nb1}); rij_etc{1+nb1} = find(~mr_A{1+nb1});
rij{1+nb1} = rij{1+nb1}(ceil(end/2):end); 
end;%for nb1=0:nbins-1;
scramble_rij = [];
ic=0;
for nb1=0:nbins-1;
scramble_rij = [scramble_rij;rij{1+nb1}+ic;rij_etc{1+nb1}+ic];
ic = ic+nrows_M_(1+nb1);
end;%for nb1=0:nbins-1;
scramble_cij_etc = find(~mc_A); 
scramble_cij = find(mc_A); scramble_cij = scramble_cij(ceil(end/2):end);
scramble_cij = [scramble_cij , scramble_cij_etc];
scramble_out_xdrop = -ones(length(scramble_rij)+length(scramble_cij),2);
scramble_out_xdrop(1:length(scramble_rij),1) = scramble_rij;
scramble_out_xdrop(length(scramble_rij) + (1:length(scramble_cij)),2) = transpose(scramble_cij);
tmp_p_= randperm(size(scramble_out_xdrop,1));
scramble_out_xdrop = scramble_out_xdrop(tmp_p_,:);
fname = sprintf('%s/scramble_out_xdrop.txt',dir__in);
fp = fopen(fname,'w');
fprintf(fp,'%d %d\n',transpose(scramble_out_xdrop-1));
fclose(fp);

fname__in = sprintf('%s/%s.in',dir__in,in__name);
fp = fopen(fname__in,'w');
fprintf(fp,'GLOBAL_verbose= 0;\n');
fprintf(fp,'GLOBAL_thread_count= 8;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','lakcluster_driver');
fprintf(fp,'GLOBAL_QR_strategy= %s;\n',QR_strategy);
fprintf(fp,'GLOBAL_QC_strategy= %s;\n',QC_strategy);
fprintf(fp,'GLOBAL_NBINS= %d;\n',nbins);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',32);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',0.50);
%fprintf(fp,'GLOBAL_TEST_sparse= 0;\n');
fprintf(fp,'GLOBAL_DIR_XPRE= %s;\n',dir__in);
fprintf(fp,'GLOBAL_A_p_name= %s;\n','A_p.mda');
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_M_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',ncols_A);
fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'mr_A_%.2d.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cind= mc_A.b16;\n');
fprintf(fp,'GLOBAL_Z_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_M_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'mr_Z_%.2d.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'GLOBAL_scramble_num= %d;\n',1);
fprintf(fp,'GLOBAL_scramble_out_xdrop_= %s;\n',sprintf('%s/scramble_out_xdrop.txt',dir__in));
fprintf(fp,'GLOBAL_scramble_rseed_= %d;\n',25);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lakcluster_test_scramble_0.m on %s;\n',date);
fclose(fp);
type(sprintf('%s/%s.in',dir__in,in__name));

call_flag = 1;
if call_flag;
string_command = sprintf('%s/../lakcluster_ver18_dev < %s/%s.in',dir_trunk,dir__in,in__name);
disp(string_command); system(string_command);
end;%if call_flag;

flag_disp=1;
if flag_disp;
% plot original matrices. ;
figure(1);clf;
colormap(colormap_pm(64));
pgap = 16;
nb1=0;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1); tmpAn = tutorial_binary_uncompress(tmpAnchar);
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1); tmpAt = tutorial_binary_uncompress(tmpAtchar);
tmpmrchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1); tmp_A_mr = tutorial_binary_uncompress(tmpmrchar)>0;
tmpmrchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1); tmp_Z_mr = tutorial_binary_uncompress(tmpmrchar)>0;
tmpmcchar = sprintf('%s/mc_A.b16',dir__in); tmp_A_mc = tutorial_binary_uncompress(tmpmcchar)>0;
subplot(2,2,1);imagesc([tmpAn(find(tmp_A_mr),find(tmp_A_mc)) ; zeros(pgap,length(find(tmp_A_mc))) ; tmpAn(find(tmp_Z_mr),find(tmp_A_mc))],[-1,+1]); title('An 00');
subplot(2,2,3);imagesc([tmpAt(find(tmp_A_mc),find(tmp_A_mr)) , zeros(length(find(tmp_A_mc)),pgap) , tmpAt(find(tmp_A_mc),find(tmp_Z_mr))],[-1,+1]); title('At 00');
disp(sprintf(' %% original nb1 %d; norm An - transpose(At): %f',nb1,norm(tmpAn-transpose(tmpAt))));
nb1=1;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1); tmpAn = tutorial_binary_uncompress(tmpAnchar);
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1); tmpAt = tutorial_binary_uncompress(tmpAtchar);
tmpmrchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1); tmp_A_mr = tutorial_binary_uncompress(tmpmrchar)>0;
tmpmrchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1); tmp_Z_mr = tutorial_binary_uncompress(tmpmrchar)>0;
tmpmcchar = sprintf('%s/mc_A.b16',dir__in); tmp_A_mc = tutorial_binary_uncompress(tmpmcchar)>0;
subplot(2,2,2);imagesc([tmpAn(find(tmp_A_mr),find(tmp_A_mc)) ; zeros(pgap,length(find(tmp_A_mc))) ; tmpAn(find(tmp_Z_mr),find(tmp_A_mc))],[-1,+1]); title('An 01');
subplot(2,2,4);imagesc([tmpAt(find(tmp_A_mc),find(tmp_A_mr)) , zeros(length(find(tmp_A_mc)),pgap) , tmpAt(find(tmp_A_mc),find(tmp_Z_mr))],[-1,+1]); title('At 01');
disp(sprintf(' %% original nb1 %d; norm An - transpose(At): %f',nb1,norm(tmpAn-transpose(tmpAt))));
end;%if flag_disp;

flag_disp=1;
if flag_disp;
% plot shuffle pre. ;
tmpchar = sprintf('%s/mc_A.b16',dir__in); mc = tutorial_binary_uncompress(tmpchar)>0;
tmpchar = sprintf('%s/mr_A_00.b16',dir__in); mr_00 = tutorial_binary_uncompress(tmpchar)>0;
tmpchar = sprintf('%s/mr_A_01.b16',dir__in); mr_01 = tutorial_binary_uncompress(tmpchar)>0;
tmpchar = 'dir_lakcluster_test_scramble/dir_lakcluster_test_scramble_out/An_00_scramble_pre.b16';
An_00_pre = tutorial_binary_uncompress(tmpchar);
tmpchar = 'dir_lakcluster_test_scramble/dir_lakcluster_test_scramble_out/An_01_scramble_pre.b16';
An_01_pre = tutorial_binary_uncompress(tmpchar);
tmpchar = 'dir_lakcluster_test_scramble/dir_lakcluster_test_scramble_out/At_00_scramble_pre.b16';
At_00_pre = tutorial_binary_uncompress(tmpchar);
tmpchar = 'dir_lakcluster_test_scramble/dir_lakcluster_test_scramble_out/At_01_scramble_pre.b16';
At_01_pre = tutorial_binary_uncompress(tmpchar);
figure(2);clf;
subplot(2,2,1); imagesc(An_00_pre(find(mr_00),:)); title('An 00 pre');
subplot(2,2,2); imagesc(An_01_pre(find(mr_01),:)); title('An 01 pre');
subplot(2,2,3); imagesc(At_00_pre(find(mc),:)); title('At 00 pre');
subplot(2,2,4); imagesc(At_01_pre(find(mc),:)); title('At 01 pre');
disp(sprintf(' %% pre nb1 00; norm An - transpose(At): %f',norm(An_00_pre(find(mr_00),:) - transpose(At_00_pre(find(mc),:)))));
disp(sprintf(' %% pre nb1 01; norm An - transpose(At): %f',norm(An_01_pre(find(mr_01),:) - transpose(At_01_pre(find(mc),:)))));
end;%if flag_disp;

flag_disp=1;
if flag_disp;
% plot shuffle pos. ;
tmpchar = sprintf('%s/mc_A.b16',dir__in); mc = tutorial_binary_uncompress(tmpchar)>0;
tmpchar = sprintf('%s/mr_A_00.b16',dir__in); mr_00 = tutorial_binary_uncompress(tmpchar)>0;
tmpchar = sprintf('%s/mr_A_01.b16',dir__in); mr_01 = tutorial_binary_uncompress(tmpchar)>0;
tmpchar = 'dir_lakcluster_test_scramble/dir_lakcluster_test_scramble_out/An_00_scramble_pos.b16';
An_00_pos = tutorial_binary_uncompress(tmpchar);
tmpchar = 'dir_lakcluster_test_scramble/dir_lakcluster_test_scramble_out/An_01_scramble_pos.b16';
An_01_pos = tutorial_binary_uncompress(tmpchar);
tmpchar = 'dir_lakcluster_test_scramble/dir_lakcluster_test_scramble_out/At_00_scramble_pos.b16';
At_00_pos = tutorial_binary_uncompress(tmpchar);
tmpchar = 'dir_lakcluster_test_scramble/dir_lakcluster_test_scramble_out/At_01_scramble_pos.b16';
At_01_pos = tutorial_binary_uncompress(tmpchar);
figure(3);clf;
subplot(2,2,1); imagesc(An_00_pos(find(mr_00),:)); title('An 00 pos');
subplot(2,2,2); imagesc(An_01_pos(find(mr_01),:)); title('An 01 pos');
subplot(2,2,3); imagesc(At_00_pos(find(mc),:)); title('At 00 pos');
subplot(2,2,4); imagesc(At_01_pos(find(mc),:)); title('At 01 pos');
disp(sprintf(' %% pos nb1 00; norm An - transpose(At): %f',norm(An_00_pos(find(mr_00),:) - transpose(At_00_pos(find(mc),:)))));
disp(sprintf(' %% pos nb1 01; norm An - transpose(At): %f',norm(An_01_pos(find(mr_01),:) - transpose(At_01_pos(find(mc),:)))));
end;%if flag_disp;


