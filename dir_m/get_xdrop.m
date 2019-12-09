function [rdrop,cdrop] = get_xdrop(lrij,lcij,gamma);
% determines how many rows (rdrop) and columns (cdrop) to remove. ;
% given number of current rows (lrij) and columns (lcij). ;
% if gamma=0, then we set remove a single row ;

if nargin<2; 
na=1; lrij_(na) = 1800; lcij_(na) = 5e4;
while (lrij_(na)>0 & lcij_(na)>0);
[rdrop_(na),cdrop_(na)] = get_xdrop(lrij_(na),lcij_(na));
lrij_(1+na) = lrij_(na) - rdrop_(na);
lcij_(1+na) = lcij_(na) - cdrop_(na);
na = na+1;
end;%while (lrij_(na)>0 & lcij_(na)>0);
clf;
subplot(1,3,1); plot(1:na,lrij_,'r-'); xlabel('iteration'); ylabel('n_r'); xlim([0,na]);
subplot(1,3,2); plot(1:na,lcij_,'b-'); xlabel('iteration'); ylabel('n_c'); xlim([0,na]);
subplot(1,3,3); plot(1:na,lcij_./lrij_,'k-'); xlabel('iteration'); ylabel('n_c/n_r'); xlim([0,na]);
set(gcf,'Position',1+[0,0,1024,512]);
disp('returning'); return;
end;%if nargin<2; 

if nargin<3; gamma = 0.0; end;

verbose=0; 
dbl1=1.000000-0.000001; %<-- same tolerance used by lakcluster. ;
gamma_tmp_row=0;
ammag_tmp_row=0;
gamma_tmp_col=0;
ammag_tmp_col=0;
rdrop=0;cdrop=0;

if (verbose); disp(sprintf(' %% [entering get_xdrop] gamma %0.3f, lrij %d lcij %d',gamma,lrij,lcij)); end;
if (lrij<=2); %<-- drop everything. ;
rdrop = lrij; cdrop = lcij;
end;%if (lrij<=2); %<-- drop everything. ;

if lrij>2;
if (gamma>0);
gamma_tmp_row = max(gamma,(dbl1)/max(1,lrij)); 
ammag_tmp_row = min(1-gamma,(lrij-1)/max(1,lrij));
end;%if (gamma>0);
if (gamma<=0);
gamma_tmp_row = (dbl1)/max(1,lrij); 
ammag_tmp_row = (lrij-1)/max(1,lrij);
end;%if (gamma<=0);
ammag_tmp_col = exp(log(ammag_tmp_row)*log(lcij)/max(1,log(max(1,lrij)))); %<-- assuming log-type balance.; 
gamma_tmp_col = dbl1-ammag_tmp_col; 
rdrop = ceil(gamma_tmp_row*lrij); 
cdrop = ceil(gamma_tmp_col*lcij); 
rdrop = min(rdrop,lrij); 
cdrop = min(cdrop,lcij);
if (verbose>0); disp(sprintf(' %% gamma_tmp_row %0.3f ammag_tmp_row %0.3f rdrop %d/%d gamma_tmp_col %0.3f ammag_tmp_col %0.3f cdrop %d/%d',gamma_tmp_row,ammag_tmp_row,rdrop,lrij,gamma_tmp_col,ammag_tmp_col,cdrop,lcij)); end;
end;%if lrij>2;

if (verbose); disp(sprintf(' %% [finished get_xdrop] rdrop %d cdrop %d',rdrop,cdrop)); end;
