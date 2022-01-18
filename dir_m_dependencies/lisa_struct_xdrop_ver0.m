function l = lisa_struct_xdrop_ver0(l);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% extract rdrop_a and cdrop_a. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
l.fname_out_xdrop_a = sprintf('%s/out_xdrop_a.txt',l.dir_out_s0); fcheck(l.fname_out_xdrop_a);
l.out_xdrop_a_ = textread(l.fname_out_xdrop_a); l.rdrop_a_ = l.out_xdrop_a_(find(l.out_xdrop_a_(:,1)>-1),1)+1; l.cdrop_a_ = l.out_xdrop_a_(find(l.out_xdrop_a_(:,2)>-1),2)+1;
l.fname_out_xdrop_b = sprintf('%s/out_xdrop_b.txt',l.dir_out_s0); fcheck(l.fname_out_xdrop_b);
l.out_xdrop_b_ = textread(l.fname_out_xdrop_b); l.rdrop_b_ = l.out_xdrop_b_(find(l.out_xdrop_b_(:,1)>-1),1)+1; l.cdrop_b_ = l.out_xdrop_b_(find(l.out_xdrop_b_(:,2)>-1),2)+1;
l.rkeep_a_ = l.rdrop_a_(end:-1:1);
l.ckeep_a_ = l.cdrop_a_(end:-1:1);
