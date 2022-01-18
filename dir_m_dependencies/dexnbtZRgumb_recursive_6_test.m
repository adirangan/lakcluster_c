clear;

if (~libisloaded('halfloop_lib'));
tmp_dir = pwd;
cd('/home/rangan/dir_bcc/dir_halfloop_dev');
loadlibrary('halfloop_lib','halfloop_lib.h','includepath','/home/rangan/dir_bcc/dir_halfloop_dev');
cd(tmp_dir);
end;%if (~libisloaded('halfloop_lib'));

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
prefix_base = "test_c"; system('rm -rf /home/rangan/dir_bcc/dir_jamison/dir_m/dir_test_c/dir_test_c_dexnbtZRgumb_g020');
E_array_base_ = E_base_rc__;
E_array_r_ij_ = 1:n_r;
E_array_c_ij_ = 1:n_c;
gamma = 0.02;
n_shuffle = 64;
p_set = 0.05;
n_member_lob = 3;
p_prev = [];
flag_force_create = 1;

dexnbtZRgumb_recursive_6( ...
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

[ ...
 binary_label_ ...
,output_label_ ...
,nlpbra_label_ ...
,nlpnex_label_ ...
] = ...
dexnbtZRgumb_recursive_6( ...
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
,0*flag_force_create ...
);

verbose=0;
flag_omp_use = 1;
flag_r0drop_vs_rcdrop = 0;
binary_label_ = cast(zeros(n_r,1),'uint64');
[~,~,~,binary_label_] = ...
calllib( ...
 'halfloop_lib' ...
,'halfloop_nonbinary_f_gateway_matlab' ...
,verbose ...
,n_r ...
,n_c ...
,cast(E_array_base_,'single') ...
,flag_r0drop_vs_rcdrop ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,pwd ...
,'test_gateway' ...
,flag_force_create ...
,flag_omp_use ...
,binary_label_ ...
);

str_output_label_ = sprintf(' "%s"',output_label_{1});
for nr=1:n_r-1; str_output_label_ = sprintf('%s , "%s"',str_output_label_,output_label_{1+nr}); end;%for nr=1:n_r-1;
str_output_label_ = sprintf('%s \n',str_output_label_);
disp(str_output_label_);
disp(sprintf('\n'));
str_nlpbra_label_ = sprintf(' "%s"',nlpbra_label_{1});
for nr=1:n_r-1; str_nlpbra_label_ = sprintf('%s , "%s"',str_nlpbra_label_,nlpbra_label_{1+nr}); end;%for nr=1:n_r-1;
str_nlpbra_label_ = sprintf('%s \n',str_nlpbra_label_);
disp(str_nlpbra_label_);
disp(sprintf('\n\n'));
str_nlpnex_label_ = sprintf(' "%s"',nlpnex_label_{1});
for nr=1:n_r-1; str_nlpnex_label_ = sprintf('%s , "%s"',str_nlpnex_label_,nlpnex_label_{1+nr}); end;%for nr=1:n_r-1;
str_nlpnex_label_ = sprintf('%s \n',str_nlpnex_label_);
disp(str_nlpnex_label_);

%%%%%%%%;
% compare outputs. ;
%%%%%%%%;
c_out_xdrop__ = transpose(c_MDA_read_i4('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/xdrop_E__.mda'));
m_out_xdrop__ = textread('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/out_xdrop_a.txt');
disp(sprintf(' %% m_out_xdrop__ vs c_out_xdrop__: %0.16f',fnorm(m_out_xdrop__ - c_out_xdrop__)/fnorm(m_out_xdrop__)));
c_QR_E__ = c_MDA_read_r8('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/QR_E__.mda');
m_out_trace_E_ = load('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/out_trace_.mat');
m_QR_E_ = m_out_trace_E_.trace_E_{1}(:,1+3);
n_iteration = numel(m_QR_E_);
m_QR_E__ = zeros(n_iteration,1+n_shuffle);
m_QR_E__(:,1) = m_QR_E_;
for nshuffle=0:n_shuffle-1;
m_QR_E__(:,1+nshuffle+1) = m_out_trace_E_.trace_E_{1+nshuffle+1}(:,1+3);
end;%for nshuffle=0:n_shuffle-1;
disp(sprintf(' %% m_QR_E__ vs c_QR_E__: %0.16f',fnorm(m_QR_E__ - c_QR_E__)/fnorm(m_QR_E__)));
%%%%%%%%;
c_out_xdrop__ = transpose(c_MDA_read_i4('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/xdrop_F__.mda'));
m_out_xdrop__ = textread('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/out_xdrop_F_a.txt');
disp(sprintf(' %% m_out_xdrop__ vs c_out_xdrop__: %0.16f',fnorm(m_out_xdrop__ - c_out_xdrop__)/fnorm(m_out_xdrop__)));
c_QR_F__ = c_MDA_read_r8('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/QR_F__.mda');
m_out_trace_F_ = load('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/out_trace_F_.mat');
m_QR_F_ = m_out_trace_F_.trace_F_{1}(:,1+3);
n_iteration = numel(m_QR_F_);
m_QR_F__ = zeros(n_iteration,1+n_shuffle);
m_QR_F__(:,1) = m_QR_F_;
for nshuffle=0:n_shuffle-1;
m_QR_F__(:,1+nshuffle+1) = m_out_trace_F_.trace_F_{1+nshuffle+1}(:,1+3);
end;%for nshuffle=0:n_shuffle-1;
disp(sprintf(' %% m_QR_F__ vs c_QR_F__: %0.16f',fnorm(m_QR_F__ - c_QR_F__)/fnorm(m_QR_F__)));
%%%%%%%%;

%%%%%%%%;
% compare outputs. ;
%%%%%%%%;
c_out_xdrop__ = transpose(c_MDA_read_i4('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/A/xdrop_E__.mda'));
m_out_xdrop__ = textread('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/dir_A/out_xdrop_a.txt');
disp(sprintf(' %% m_out_xdrop__ vs c_out_xdrop__: %0.16f',fnorm(m_out_xdrop__ - c_out_xdrop__)/fnorm(m_out_xdrop__)));
c_QR_E__ = c_MDA_read_r8('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/A/QR_E__.mda');
m_out_trace_E_ = load('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/dir_A/out_trace_.mat');
m_QR_E_ = m_out_trace_E_.trace_E_{1}(:,1+3);
n_iteration = numel(m_QR_E_);
m_QR_E__ = zeros(n_iteration,1+n_shuffle);
m_QR_E__(:,1) = m_QR_E_;
for nshuffle=0:n_shuffle-1;
m_QR_E__(:,1+nshuffle+1) = m_out_trace_E_.trace_E_{1+nshuffle+1}(:,1+3);
end;%for nshuffle=0:n_shuffle-1;
disp(sprintf(' %% m_QR_E__ vs c_QR_E__: %0.16f',fnorm(m_QR_E__ - c_QR_E__)/fnorm(m_QR_E__)));
%%%%%%%%;
c_out_xdrop__ = transpose(c_MDA_read_i4('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/B/xdrop_E__.mda'));
m_out_xdrop__ = textread('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/dir_B/out_xdrop_a.txt');
disp(sprintf(' %% m_out_xdrop__ vs c_out_xdrop__: %0.16f',fnorm(m_out_xdrop__ - c_out_xdrop__)/fnorm(m_out_xdrop__)));
c_QR_E__ = c_MDA_read_r8('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/B/QR_E__.mda');
m_out_trace_E_ = load('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/dir_B/out_trace_.mat');
m_QR_E_ = m_out_trace_E_.trace_E_{1}(:,1+3);
n_iteration = numel(m_QR_E_);
m_QR_E__ = zeros(n_iteration,1+n_shuffle);
m_QR_E__(:,1) = m_QR_E_;
for nshuffle=0:n_shuffle-1;
m_QR_E__(:,1+nshuffle+1) = m_out_trace_E_.trace_E_{1+nshuffle+1}(:,1+3);
end;%for nshuffle=0:n_shuffle-1;
disp(sprintf(' %% m_QR_E__ vs c_QR_E__: %0.16f',fnorm(m_QR_E__ - c_QR_E__)/fnorm(m_QR_E__)));
%%%%%%%%;
c_out_xdrop__ = transpose(c_MDA_read_i4('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/B/xdrop_F__.mda'));
m_out_xdrop__ = textread('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/dir_B/out_xdrop_F_a.txt');
disp(sprintf(' %% m_out_xdrop__ vs c_out_xdrop__: %0.16f',fnorm(m_out_xdrop__ - c_out_xdrop__)/fnorm(m_out_xdrop__)));
c_QR_F__ = c_MDA_read_r8('/home/rangan/dir_bcc/dir_halfloop_dev/dir_test_rc/dir_test_rc_g020_n03/B/QR_F__.mda');
m_out_trace_F_ = load('./dir_test_c/dir_test_c_dexnbtZRgumb_g020/dir_B/out_trace_F_.mat');
m_QR_F_ = m_out_trace_F_.trace_F_{1}(:,1+3);
disp(num2str(m_QR_F_([1,2,3,end-2,end-1,end-0])));
n_iteration = numel(m_QR_F_);
m_QR_F__ = zeros(n_iteration,1+n_shuffle);
m_QR_F__(:,1) = m_QR_F_;
for nshuffle=0:n_shuffle-1;
m_QR_F__(:,1+nshuffle+1) = m_out_trace_F_.trace_F_{1+nshuffle+1}(:,1+3);
end;%for nshuffle=0:n_shuffle-1;
disp(sprintf(' %% m_QR_F__ vs c_QR_F__: %0.16f',fnorm(m_QR_F__ - c_QR_F__)/fnorm(m_QR_F__)));
%%%%%%%%;

