%%%%%%%%;
% Warning, does not pass flag_r0drop_vs_rcdrop. ;
%%%%%%%%;
function ...
[ ...
parameter ...
,output_label_ ...
,nlpbra_label_ ...
,nlpnex_label_ ...
] = ...
halfloop_nonbinary_recursive_gateway_c( ...
 parameter ...
,str_code ...
,dir_trunk ...
,E_base_rc__ ...
,prefix_base ...
);

if (nargin<1);
%%%%%%%%;
disp(sprintf(' %% testing halfloop_nonbinary_recursive_gateway_c'));
M =  178; N =  2e3; n_cluster = 1; n_rank =  2; rseed = 0;
snr_ = transpose(linspace(0.5,1.0,1+16)); n_snr = numel(snr_);
nlpv_ = zeros(n_snr,1);
%%%%%%%%;
for nsnr=0:n_snr-1;
snr = snr_(1+nsnr);
disp(sprintf(' %% nsnr %d/%d snr %0.2f',nsnr,n_snr,snr));
rng(rseed);
[A_n_,label_A_,n_label_A_,pf_,pi_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(1);clf;figmed;figbeach();
imagesc(A_n_(pi_{1},pi_{2})); axis image; axisnotick;
end;%if flag_plot;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_force_create = 1;
parameter.flag_plot = 0;
str_code = '/home/rangan/dir_bcc/dir_halfloop_dev/halfloop_dev';
dir_trunk = pwd;
E_base_rc__ = A_n_;
prefix_base = 'test';
[ ...
parameter ...
,output_label_ ...
,nlpbra_label_ ...
,nlpnex_label_ ...
] = ...
halfloop_nonbinary_recursive_gateway_c( ...
 parameter ...
,str_code ...
,dir_trunk ...
,E_base_rc__ ...
,prefix_base ...
);
%%%%%%%%;
label_B_ = label_str_to_enum_1(output_label_);
[lpv,lP_0,flag_method,cap_,cup_] = label_to_label_enrichment_quad_4(label_A_,label_B_);
nlpv_(1+nsnr) = -lpv;
end;%for nsnr=0:n_snr-1;
%%%%%%%%;
figure(1);figsml;
plot(snr_,nlpv_,'ko-','LineWidth',2,'MarkerFaceColor','r'); xlabel('snr'); ylabel('nlpv');
title(sprintf('(M %d,N %d), recovery for single planted cluster',M,N),'Interpreter','none');
%%%%%%%%;
disp('returning');return;
end;%if (nargin<1);

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'verbose')); parameter.verbose = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'halfloop_recursion_limit')); parameter.halfloop_recursion_limit = 63; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_r0drop_vs_rcdrop')); parameter.flag_r0drop_vs_rcdrop = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'gamma')); parameter.gamma = 0.01; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_shuffle')); parameter.n_shuffle = 64; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_orth_brute')); parameter.flag_orth_brute = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'p_set')); parameter.p_set = 0.05; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_member_lob')); parameter.n_member_lob = 2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_force_create')); parameter.flag_force_create = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_omp_use')); parameter.flag_omp_use = 1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_plot')); parameter.flag_plot = 1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_replot')); parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
rseed = parameter.rseed;
verbose = parameter.verbose;
halfloop_recursion_limit = parameter.halfloop_recursion_limit;
flag_r0drop_vs_rcdrop = parameter.flag_r0drop_vs_rcdrop;
gamma = parameter.gamma;
n_shuffle = parameter.n_shuffle;
flag_orth_brute = parameter.flag_orth_brute;
p_set = parameter.p_set;
n_member_lob = parameter.n_member_lob;
flag_force_create = parameter.flag_force_create;
flag_omp_use = parameter.flag_omp_use;
flag_plot = parameter.flag_plot;
flag_replot = parameter.flag_replot;

if (verbose); disp(sprintf(' %% [entering halfloop_nonbinary_recursive_gateway_c]')); end;
prefix_base0 = prefix_base;
dir_0in = sprintf('%s/dir_%s',dir_trunk,prefix_base0);
prefix_gamma = ''; if (gamma> 0); prefix_gamma = sprintf('_g%.3d',floor(1000*gamma)); end;
prefix_n_member_lob = ''; if (n_member_lob> 2); prefix_n_member_lob = sprintf('_n%.2d',n_member_lob); end;
prefix_base1 = sprintf('%s%s%s',prefix_base0,prefix_gamma,prefix_n_member_lob);
dir_out = sprintf('%s/dir_%s',dir_0in,prefix_base1);
fname_output_label = sprintf('%s/output_label__.txt',dir_out);
fname_nlpbra_label = sprintf('%s/nlpbra_label__.txt',dir_out);
fname_nlpnex_label = sprintf('%s/nlpnex_label__.txt',dir_out);

flag_skip = exist(fname_output_label,'file') & exist(fname_nlpbra_label,'file') & exist(fname_nlpnex_label,'file');

if (flag_force_create | ~flag_skip);
if (verbose); disp(sprintf(' %% %s not found, creating',fname_output_label)); end;
%%%%%%%%;
% Write mda_r4 file. ;
%%%%%%%%;
fname_mda_r4 = sprintf('%s/%s_E_base_rc__.mda',dir_trunk,prefix_base);
c_MDA_write_r4(fname_mda_r4,full(E_base_rc__));
%%%%%%%%;
% Now loop over recursion_limit. ;
%%%%%%%%;
tmp_halfloop_recursion_limit = 0;
flag_continue = 1;
while (flag_continue);
%%%%%%%%;
% define input file. ;
%%%%%%%%;
fname_0in = sprintf('%s/halfloop_nonbinary_f_gateway_shell_%s_l%d.in',dir_trunk,prefix_base,tmp_halfloop_recursion_limit);
fp = fopen(fname_0in,'w');
fprintf(fp,'GLOBAL_mode= halfloop_nonbinary_f_gateway_shell;\n');
fprintf(fp,'GLOBAL_verbose= %d;\n',verbose);
fprintf(fp,'GLOBAL_flag_omp_use= %d;\n',flag_omp_use);
fprintf(fp,'GLOBAL_halfloop_recursion_limit= %d;\n',tmp_halfloop_recursion_limit);
fprintf(fp,'GLOBAL_E_base_mda_r4= %s;\n',fname_mda_r4);
fprintf(fp,'GLOBAL_flag_r0drop_vs_rcdrop= 0;\n',flag_r0drop_vs_rcdrop);
fprintf(fp,'GLOBAL_gamma= %0.3f;\n',gamma);
fprintf(fp,'GLOBAL_n_shuffle= %d;\n',n_shuffle);
fprintf(fp,'GLOBAL_flag_orth_brute= %d;\n',flag_orth_brute);
fprintf(fp,'GLOBAL_p_set= %0.3f;\n',p_set);
fprintf(fp,'GLOBAL_n_member_lob= %d;\n',n_member_lob);
fprintf(fp,'GLOBAL_dir_trunk= %s;\n',dir_trunk);
fprintf(fp,'GLOBAL_prefix_base= %s;\n',prefix_base0);
fprintf(fp,'GLOBAL_flag_force_create= %d;\n',flag_force_create);
fprintf(fp,'END= 0;\n');
fclose(fp);
%%%%;
if (verbose); type(fname_0in); end;
%%%%%%%%;
% Now run code. ;
%%%%%%%%;
str_command = sprintf('%s < %s;',str_code,fname_0in);
flag_error = system(str_command);
%%%%%%%%;
% Now check for output. ;
%%%%%%%%;
flag_output_exist = exist(fname_output_label,'file') & exist(fname_nlpbra_label,'file') & exist(fname_nlpnex_label,'file');
%%%%%%%%;
% load output. ;
%%%%%%%%;
if flag_output_exist;
if (nargout>1);
fp = fopen(fname_output_label,'r'); output_label_ = textscan(fp,'%s','Delimiter','\n'); fclose(fp); output_label_ = output_label_{1};
fp = fopen(fname_nlpbra_label,'r'); nlpbra_label_ = textscan(fp,'%s','Delimiter','\n'); fclose(fp); nlpbra_label_ = nlpbra_label_{1};
fp = fopen(fname_nlpnex_label,'r'); nlpnex_label_ = textscan(fp,'%s','Delimiter','\n'); fclose(fp); nlpnex_label_ = nlpnex_label_{1};
end;%if (nargout>1);
end;%if flag_output_exist;
%%%%%%%%;
% Now check for success. ;
%%%%%%%%;
flag_success = ~flag_error & flag_output_exist;
tmp_halfloop_recursion_limit = tmp_halfloop_recursion_limit + 1;
flag_continue = flag_success & (tmp_halfloop_recursion_limit<=halfloop_recursion_limit);
%%%%%%%%;
% Now delete input file if we are not at the final recursion;
%%%%%%%%;
if (flag_continue); delete(fname_0in); end;
%%%%%%%%;
end;%while (flag_continue);
%%%%%%%%;
% Now delete mda_r4 only if final recursion was successful. ;
if  (flag_success & (tmp_halfloop_recursion_limit> halfloop_recursion_limit));
delete(fname_mda_r4);
else;
disp(sprintf(' %% Warning, %s not completed',fname_0in));
end;%if  (flag_success & (tmp_halfloop_recursion_limit> halfloop_recursion_limit));
%%%%%%%%;
%%%%%%%%;
end;%if (flag_force_create | ~flag_skip);

%%%%%%%%;
% Now load what we have so far, even if it is not the final recursion. ;
%%%%%%%%;
flag_success = exist(fname_output_label,'file') & exist(fname_nlpbra_label,'file') & exist(fname_nlpnex_label,'file');
if flag_success;
%%%%%%%%;
% load output. ;
%%%%%%%%;
if (nargout>1);
fp = fopen(fname_output_label,'r'); output_label_ = textscan(fp,'%s','Delimiter','\n'); fclose(fp); output_label_ = output_label_{1};
fp = fopen(fname_nlpbra_label,'r'); nlpbra_label_ = textscan(fp,'%s','Delimiter','\n'); fclose(fp); nlpbra_label_ = nlpbra_label_{1};
fp = fopen(fname_nlpnex_label,'r'); nlpnex_label_ = textscan(fp,'%s','Delimiter','\n'); fclose(fp); nlpnex_label_ = nlpnex_label_{1};
end;%if (nargout>1);
%%%%%%%%;
if flag_plot;
fname_fig = sprintf('%s/tree_FIGA',dir_out);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
fp = fopen(fname_output_label,'r'); output_label_ = textscan(fp,'%s','Delimiter','\n'); fclose(fp); output_label_ = output_label_{1};
fp = fopen(fname_nlpbra_label,'r'); nlpbra_label_ = textscan(fp,'%s','Delimiter','\n'); fclose(fp); nlpbra_label_ = nlpbra_label_{1};
fp = fopen(fname_nlpnex_label,'r'); nlpnex_label_ = textscan(fp,'%s','Delimiter','\n'); fclose(fp); nlpnex_label_ = nlpnex_label_{1};
label_plot_recursive_2(output_label_,nlpbra_label_,nlpnex_label_,[],[]);
sgtitle(sprintf('%s',dir_out),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;
%%%%%%%%;
end;%if flag_success;

if (verbose); disp(sprintf(' %% [finished halfloop_nonbinary_recursive_gateway_c]')); end;
