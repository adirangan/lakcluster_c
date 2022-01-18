function l = lisa_struct_xdrop_shuffle_ver0(l,nshuffle);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% extract rdrop_a and cdrop_a from the shuffled trace. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if (nshuffle>0); 
l.fname_out_xdrop_a_shuffle = sprintf('%s/out_xdrop_a_s%.4d.txt',l.dir_out_trace,nshuffle); fcheck(l.fname_out_xdrop_a_shuffle);
 else
l.fname_out_xdrop_a_shuffle = sprintf('%s/out_xdrop_a.txt',l.dir_out_s0); fcheck(l.fname_out_xdrop_a_shuffle);
end;%if (nshuffle>0); 
l.out_xdrop_a_shuffle_ = textread(l.fname_out_xdrop_a_shuffle); 
l.rdrop_a_shuffle_ = l.out_xdrop_a_shuffle_(find(l.out_xdrop_a_shuffle_(:,1)>-1),1)+1; 
l.cdrop_a_shuffle_ = l.out_xdrop_a_shuffle_(find(l.out_xdrop_a_shuffle_(:,2)>-1),2)+1;
l.rkeep_a_shuffle_ = l.rdrop_a_shuffle_(end:-1:1);
l.ckeep_a_shuffle_ = l.cdrop_a_shuffle_(end:-1:1);
