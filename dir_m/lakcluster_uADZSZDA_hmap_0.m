function lakcluster_uADZSZDA_hmap_0(dir_code,dir_trunk,prefix,M_n_,rev_flag,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind,gamma,B_MLT,Ireq,shuffle_num,verbose_flag);
% plots heatmap ;
% For example, try: ; 
%{

  dir_code = '/data/rangan/dir_bcc/dir_code_122417/';
  dir_trunk = '/data/rangan/dir_bcc/dir_GSE93601/';
  prefix = 'GSE93601_A_t_xf_';
  nbins = 1;
  M_n_ = cell(nbins,1); T_n_ = cell(nbins,1);
  A_n_rind_ = cell(nbins,1);
  Z_n_rind_ = cell(nbins,1);
  M_n_{1} = [ A_tlp_ ; A_tgp_ ; A_tln_ ; A_tgn_ ]; 
  A_n_rind_{1} = 1:sum(M_tlp_ | M_tgp_); 
  Z_n_rind_{1} = sum(M_tlp_ | M_tgp_) + (1:sum(M_tln_ | M_tgn_));
  T_n_{1} = ones(size(M_n_{1},1),1);
  rev_flag = 0;
  A_n_cind = 1:ncols;
  T_n_cind = 1;
  gamma = 0.01;
  B_MLT = 32;
  Ireq = 0;
  shuffle_num = 0;
  verbose_flag = 0;
  lakcluster_uADZSZDA_hmap_0(dir_code,dir_trunk,prefix,M_n_,rev_flag,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind,gamma,B_MLT,Ireq,shuffle_num,verbose_flag);

  %}

nbins=length(M_n_);
[A_n_cols,Y_n_cols,T_n_cols,~,~] = lakcluster_uADZSZDA_check_0(shuffle_num,M_n_,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind);

tmp_string = sprintf('%s_%s',prefix,lakcluster_uADZSZDA_xfix_gen_ver0(rev_flag,gamma,B_MLT,Ireq,0));
dir__in = sprintf('%s/dir_%s',dir_trunk,prefix);
dir_out = sprintf('%s/dir_%s',dir__in,tmp_string); 
if (verbose_flag); disp(sprintf(' tmp_string: %s',tmp_string)); end;%verbose;
if (verbose_flag); disp(sprintf(' dir__in: %s',dir__in)); end;%verbose;
if (verbose_flag); disp(sprintf(' dir_out: %s',dir_out)); end;%verbose;
tmpchar = sprintf('%s/out_xdrop_a.txt',dir_out);

out_xdrop_a = textread(tmpchar);
rij = out_xdrop_a(:,1); rij = rij(find(rij>-1)); rij = rij(end:-1:1); r_ij = 1+rij;
cij = out_xdrop_a(:,2); cij = cij(find(cij>-1)); cij = cij(end:-1:1); c_ij = 1+cij;
M_all_n = []; C_all_n = []; M_X___n = []; C_X___n = []; cov_len = ceil(log2(nbins));
for nb1=0:nbins-1;
M_all_n = [M_all_n ; M_n_{1+nb1}];
if (nbins>1); C_all_n = [C_all_n ; repmat((dec2bin(nb1,cov_len)+0-'0'),size(M_n_{1+nb1},1),1)]; end;
if (rev_flag==0); 
M_X___n = [M_X___n ; M_n_{1+nb1}(Z_n_rind_{1+nb1},:) ]; 
if (nbins>1); C_X___n = [C_X___n ; repmat((dec2bin(nb1,cov_len)+0-'0'),length(Z_n_rind_{1+nb1}),1)]; end;
end;%if (rev_flag==0); 
if (rev_flag==1); 
M_X___n = [M_X___n ; M_n_{1+nb1}(A_n_rind_{1+nb1},:) ]; 
if (nbins>1); C_X___n = [C_X___n ; repmat((dec2bin(nb1,cov_len)+0-'0'),length(A_n_rind_{1+nb1}),1)]; end;
end;%if (rev_flag==1); 
end;%for nb1=0:nbins-1;

pgap = 4; ggap = 128; cgap = 256;
if (nbins>1);
C_rep_n = zeros(cov_len,cov_len*cgap);
for nc=1:cov_len;
C_rep_n(cov_len,cgap*(cov_len-1) + (1:cgap)) = 1;
end;%for nc=1:cov_len;
E_all_n = C_all_n*C_rep_n;
E_X___n = C_X___n*C_rep_n;
end;%if (nbins>1);
cmap_autumn = colormap('autumn'); 
cmap_autumn(1,:) = [1,1,1]; 
cmap_autumn(end-1,:) = 0.85*[1,1,1]; 
cmap_autumn(end-0,:) = 0.00*[1,1,1]; 
clen = size(cmap_autumn,1); 
colormap(cmap_autumn);
cbot=-1.5; ctop=1.5; cvals = linspace(cbot,ctop,clen); 
cbotl=cvals(2); ctopl = cvals(end-2); cdiff = cvals(end)-cvals(end-1); 
%subplot(1,1,1);
if (nbins>1); A11 = ctopl + cdiff + cdiff*E_all_n(r_ij,:); end;
if (nbins>1); A12 = cbot*ones(length(r_ij),ggap); end;
A13 = min(ctopl,max(cbotl,1*M_all_n(r_ij,c_ij)));
if (nbins>1); A21 = cbot*ones(pgap,cov_len*cgap); end;
if (nbins>1); A22 = cbot*ones(pgap,ggap); end;
A23 = cbot*ones(pgap,length(c_ij));
if (nbins>1); A31 = ctopl + cdiff + cdiff*E_X___n(:,:); end;
if (nbins>1); A32 = cbot*ones(size(M_X___n,1),ggap); end;
A33 = min(ctopl,max(cbotl,1*M_X___n(:,c_ij)));
%disp([size(A11),size(A12),size(A13);size(A21),size(A22),size(A23);size(A31),size(A32),size(A33)]);
if (nbins>1); imagesc( [ A11 , A12 , A13 ; A21 , A22 , A23 ; A31 , A32 , A33 ] , [cbot,ctop] ); end;
if (nbins==1); imagesc( [ A13 ; A23 ; A33 ] , [cbot,ctop] ); end;
set(gca,'Xtick',[],'Ytick',[]);axis off; 
set(gcf,'Position',[100,100,1124,868/1.5]);


