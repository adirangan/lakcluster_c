function [rdrop_,cdrop_,lrij_,lcij_] = get_q3drop_(lrij,lcij,gamma);
% determines how many rows (rdrop) and columns (cdrop) to remove. ;
% given number of current rows (lrij) and columns (lcij). ;
% if gamma=0, then we set remove a single row ;

if nargin<2; 
na=1; lrij_(na) = 1800; lcij_(na) = 5e4;
while (lrij_(na)>0 & lcij_(na)>0);
[rdrop_(na),cdrop_(na)] = get_q3drop(lrij_(na),lcij_(na));
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

na=1; lrij_(na) = lrij; lcij_(na) = lcij;
while (lrij_(na)>0 & lcij_(na)>0);
[rdrop_(na),cdrop_(na)] = get_q3drop(lrij_(na),lcij_(na),gamma);
lrij_(1+na) = lrij_(na) - rdrop_(na);
lcij_(1+na) = lcij_(na) - cdrop_(na);
na = na+1;
end;%while (lrij_(na)>0 & lcij_(na)>0);
