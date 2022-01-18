function [rdrop_,cdrop_,rkeep_,ckeep_,xdrop_,xkeep_] = lisa_xdrop_ver0(fname);
% extract rows and cols from xdrop. ;
% Note that these indices are updated to range from 1:max (i.e., matlab-convention), ;
% and not from 0:max-1 (i.e., c-convention). ;

fcheck(fname);
xdrop_ = textread(fname); 
rdrop_ = xdrop_(find(xdrop_(:,1)>-1),1)+1; 
cdrop_ = xdrop_(find(xdrop_(:,2)>-1),2)+1;
rkeep_ = rdrop_(end:-1:1);
ckeep_ = cdrop_(end:-1:1);
xkeep_ = xdrop_(end:-1:1,:);
