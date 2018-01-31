function hmap_autumn_0(M_n_,A_n_rind_,Z_n_rind_,rev_flag,pgap);
% basic hmap ;
if nargin<5; pgap = 16; end;

nbins=length(M_n_);

M_all_n = []; M_D___n = []; M_X___n = []; 
for nb1=0:nbins-1;
%M_all_n = [M_all_n ; M_n_{1+nb1}];
if (rev_flag==0); 
M_D___n = [M_D___n ; M_n_{1+nb1}(A_n_rind_{1+nb1},:) ]; 
M_X___n = [M_X___n ; M_n_{1+nb1}(Z_n_rind_{1+nb1},:) ]; 
end;%if (rev_flag==0); 
if (rev_flag==1); 
M_D___n = [M_D___n ; M_n_{1+nb1}(Z_n_rind_{1+nb1},:) ]; 
M_X___n = [M_X___n ; M_n_{1+nb1}(A_n_rind_{1+nb1},:) ]; 
end;%if (rev_flag==1); 
end;%for nb1=0:nbins-1;

cmap_autumn = colormap('autumn'); 
cmap_autumn(1,:) = [1,1,1]; 
cmap_autumn(end-1,:) = 0.85*[1,1,1]; 
cmap_autumn(end-0,:) = 0.00*[1,1,1]; 
clen = size(cmap_autumn,1); 
colormap(cmap_autumn);
cbot=-1.5; ctop=1.5; cvals = linspace(cbot,ctop,clen); 
cbotl=cvals(2); ctopl = cvals(end-2); cdiff = cvals(end)-cvals(end-1); 
%subplot(1,1,1);
A13 = min(ctopl,max(cbotl,1*M_D___n));
A23 = cbot*ones(pgap,size(M_D___n,2));
A33 = min(ctopl,max(cbotl,1*M_X___n));
if (size(M_X___n,1)>0);
imagesc( [ A13 ; A23 ; A33 ] , [cbot,ctop] );
end;%if (size(M_X___n,1)>0);
if (size(M_X___n,1)==0);
imagesc( [ A13 ] , [cbot,ctop] );
end;%if (size(M_X___n,1)==0);
set(gca,'Xtick',[],'Ytick',[]);axis off; 
set(gcf,'Position',[100,100,1124,868/1.5]);
