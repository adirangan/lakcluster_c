function [i,j] = eind2sub(siz,ind);
[i,j] = ind2sub(siz,1+ind); i=i-1; j=j-1;
