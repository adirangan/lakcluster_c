function lakcluster_uADZSZDA_C_hmap_0(dir_code,dir_trunk,prefix,M_n_,rev_flag,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind,C_n_,gamma,B_MLT,Ireq,shuffle_num,verbose_flag);
% plots heatmap ;
% Must provide cell-array of covariate-categories C_n_ corresponding to cell-array M_n_. ;

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
M_all_n = []; C_all_n = []; M_X___n = []; C_X___n = []; 
for nb1=0:nbins-1;
M_all_n = [M_all_n ; M_n_{1+nb1}];
C_all_n = [C_all_n ; C_n_{1+nb1}];
if (rev_flag==0); 
M_X___n = [M_X___n ; M_n_{1+nb1}(Z_n_rind_{1+nb1},:) ]; 
C_X___n = [C_X___n ; C_n_{1+nb1}(Z_n_rind_{1+nb1})];
end;%if (rev_flag==0); 
if (rev_flag==1); 
M_X___n = [M_X___n ; M_n_{1+nb1}(A_n_rind_{1+nb1},:) ]; 
C_X___n = [C_X___n ; C_n_{1+nb1}(A_n_rind_{1+nb1})];
end;%if (rev_flag==1); 
end;%for nb1=0:nbins-1;

cov_len = 1;
pgap = 4; ggap = 128; cgap = 256;
C_rep_n = zeros(cov_len,cov_len*cgap);
for nc=1:cov_len;
C_rep_n(cov_len,cgap*(cov_len-1) + (1:cgap)) = 1;
end;%for nc=1:cov_len;
E_all_n = C_all_n*C_rep_n;
E_X___n = C_X___n*C_rep_n;
cmap_autumn = colormap('autumn'); 
cmap_autumn(1,:) = [1,1,1]; 
cmap_autumn(end-1,:) = 0.85*[1,1,1]; 
cmap_autumn(end-0,:) = 0.00*[1,1,1]; 
clen = size(cmap_autumn,1); 
colormap(cmap_autumn);
cbot=-1.5; ctop=1.5; cvals = linspace(cbot,ctop,clen); 
cbotl=cvals(2); ctopl = cvals(end-2); cdiff = cvals(end)-cvals(end-1); 
%subplot(1,1,1);
A11 = ctopl + cdiff + cdiff*E_all_n(r_ij,:);
A12 = cbot*ones(length(r_ij),ggap);
A13 = min(ctopl,max(cbotl,1*M_all_n(r_ij,c_ij)));
A21 = cbot*ones(pgap,cov_len*cgap);
A22 = cbot*ones(pgap,ggap);
A23 = cbot*ones(pgap,length(c_ij));
A31 = ctopl + cdiff + cdiff*E_X___n(:,:);
A32 = cbot*ones(size(M_X___n,1),ggap);
A33 = min(ctopl,max(cbotl,1*M_X___n(:,c_ij)));
%disp([size(A11),size(A12),size(A13);size(A21),size(A22),size(A23);size(A31),size(A32),size(A33)]);
imagesc( [ A11 , A12 , A13 ; A21 , A22 , A23 ; A31 , A32 , A33 ] , [cbot,ctop] );
set(gca,'Xtick',[],'Ytick',[]);axis off; 
set(gcf,'Position',[100,100,1124,868/1.5]);


