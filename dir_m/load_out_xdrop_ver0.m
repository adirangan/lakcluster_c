function ...
[ ...
 xdrop_ ...
] = ...
load_out_xdrop_ver0( ...
fname ...
);

if (~exist(fname,'file')); disp(sprintf(' %% Warning, %s not found in load_out_xdrop_ver0',fname)); end;

out_xdrop_ = textread(fname);
index_rdrop_ = out_xdrop_(find(out_xdrop_(:,1)>-1),1);
index_cdrop_ = out_xdrop_(find(out_xdrop_(:,2)>-1),2);
index_rkeep_ = index_rdrop_(end:-1:1);
index_ckeep_ = index_cdrop_(end:-1:1);
ij_rdrop_ = index_rdrop_ + 1;
ij_cdrop_ = index_cdrop_ + 1;
ij_rkeep_ = index_rkeep_ + 1;
ij_ckeep_ = index_ckeep_ + 1;

xdrop_ = struct('type','xdrop');
xdrop_.fname = fname;
xdrop_.index_rdrop_ = index_rdrop_;
xdrop_.index_cdrop_ = index_cdrop_;
xdrop_.index_rkeep_ = index_rkeep_;
xdrop_.index_ckeep_ = index_ckeep_;
xdrop_.ij_rdrop_ = ij_rdrop_;
xdrop_.ij_cdrop_ = ij_cdrop_;
xdrop_.ij_rkeep_ = ij_rkeep_;
xdrop_.ij_ckeep_ = ij_ckeep_;
xdrop_.out_xdrop_ = out_xdrop_;






 
