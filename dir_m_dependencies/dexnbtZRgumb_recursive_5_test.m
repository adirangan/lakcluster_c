clear;

n_r = 128;
n_c = 512;
E_base_rc__ = zeros(n_r,n_c);
x_ = linspace(-1,+1,n_r);
y_ = linspace(-1,+1,n_c);
for nc=0:n_c-1;
for nr=0:n_r-1;
x = x_(1+nr);
y = y_(1+nc);
z = (x+y)/4;
E_base_rc__(1+nr,1+nc) = sin(2*pi*x) + cos(2*pi*2*y) + x^2 + y^3 + cos(2*pi*4*z);
end;%for nr=0:n_r-1;
end;%for nc=0:n_c-1;
disp(E_base_rc__([1,2,127,128],[1,2,511,512]));
tmp__ = mean_center_0(E_base_rc__);
disp(tmp__([1,2,127,128],[1,2,511,512]));
%imagesc(E_base_rc__);

dir_trunk = pwd;
dir_out = [];
prefix_base = "test_c";
E_array_base_ = E_base_rc__;
E_array_r_ij_ = 1:n_r;
E_array_c_ij_ = 1:n_c;
gamma = 0;
n_shuffle = 64;
p_set = 0.05;
n_member_lob = 3;
p_prev = [];
flag_force_create = [];

dexnbtZRgumb_recursive_5( ...
 dir_trunk ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,p_prev ...
,flag_force_create ...
);

out_trace_E_ = load('./dir_test_c/dir_test_c_dexnbtZRgumb_g000/out_trace_.mat');
QR_ = out_trace_E_.trace_E_{1}(:,1+3);
disp(QR_([1,2,end-1,end]));
n_iteration = numel(QR_);
QR__ = zeros(n_iteration,1+n_shuffle);
QR__(:,1) = QR_;
for nshuffle=0:n_shuffle-1;
QR__(:,1+nshuffle+1) = out_trace_E_.trace_E_{1+nshuffle+1}(:,1+3);
end;%for nshuffle=0:n_shuffle-1;
ZR__ = zeros(n_iteration,1+n_shuffle);
for niteration=0:n_iteration-1;
QR_avg = mean(QR__(1+niteration,2:end));
QR_std = std(QR__(1+niteration,2:end),1);
ZR__(1+niteration,:) = (QR__(1+niteration,:) - QR_avg)/max(1e-12,QR_std);
end;%for niteration=0:n_iteration-1;
p_use= p_set/(1+2*p_set);
nlp_ZR_ = -z_to_lp(ZR__(:,1));
[zone_max,zone_max_index] = Z_imax_zerobased(0,nlp_ZR_,-log(p_use));
disp([zone_max , zone_max_index]);
