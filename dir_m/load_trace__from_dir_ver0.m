function ...
[ ...
 trace__ ...
] = ...
load_trace__from_dir_ver0( ...
 dir_out_trace ...
);

nshuffle=0;
fname = sprintf('%s/out_trace_s%.4d.txt',dir_out_trace,nshuffle);
[ ...
 niter_s0000_ ...
,r_rem_s0000_ ...
,c_rem_s0000_ ...
,QR_s0000_ ...
,QC_s0000_ ...
,I_rem_s0000_ ...
,trace_s0000_ ...
] = ...
load_trace_ver0( ...
 fname ...
);
n_iter = numel(niter_s0000_);

nshuffle=0;
flag_continue=1;
while flag_continue;
nshuffle = nshuffle + 1;
fname = sprintf('%s/out_trace_s%.4d.txt',dir_out_trace,nshuffle);
flag_exist = exist(fname,'file');
if flag_exist;
flag_match = 1;
[ niter_ ] = load_trace_ver0( fname );
if (numel(niter_)==n_iter); flag_match=1; end;
if (numel(niter_)~=n_iter); flag_match=0; end;
if ~flag_match; flag_continue=0; end;
end;%if flag_exist;
if ~flag_exist; flag_continue=0; end;
end;%while flag_continue;
n_shuffle = nshuffle-1;
if (n_shuffle<1); disp(sprintf(' %% Warning, n_shuffle %d in load_trace__from_dir_ver0',n_shuffle)); end;

trace_st__ = cell(1+n_shuffle,1);
niter_is__ = zeros(n_iter,1+n_shuffle);
r_rem_is__ = zeros(n_iter,1+n_shuffle);
c_rem_is__ = zeros(n_iter,1+n_shuffle);
QR_is__ = zeros(n_iter,1+n_shuffle);
QC_is__ = zeros(n_iter,1+n_shuffle);
I_rem_is__ = zeros(n_iter,1+n_shuffle);
for nshuffle=0:n_shuffle;
fname = sprintf('%s/out_trace_s%.4d.txt',dir_out_trace,nshuffle);
[ ...
 niter_ ...
,r_rem_ ...
,c_rem_ ...
,QR_ ...
,QC_ ...
,I_rem_ ...
,trace_ ...
] = ...
load_trace_ver0( ...
 fname ...
);
niter_is__(:,1+nshuffle) = niter_;
r_rem_is__(:,1+nshuffle) = r_rem_;
c_rem_is__(:,1+nshuffle) = c_rem_;
QR_is__(:,1+nshuffle) = QR_;
QC_is__(:,1+nshuffle) = QC_;
I_rem_is__(:,1+nshuffle) = I_rem_;
trace_st__{1+nshuffle} = trace_;
end;%for nshuffle=0:n_shuffle;
%%%%;
QR_avg_i_ = mean(QR_is__(:,2:end),2);
QR_std_i_ = std(QR_is__(:,2:end),1,2);
ZR_is__ = bsxfun(@rdivide,bsxfun(@minus,QR_is__,QR_avg_i_),max(1e-12,QR_std_i_));
nlpR_is__ = -z_to_lp(ZR_is__);
%%%%;
QC_avg_i_ = mean(QC_is__(:,2:end),2);
QC_std_i_ = std(QC_is__(:,2:end),1,2);
ZC_is__ = bsxfun(@rdivide,bsxfun(@minus,QC_is__,QC_avg_i_),max(1e-12,QC_std_i_));
nlpC_is__ = -z_to_lp(ZC_is__);
%%%%;

%%%%;
trace__ = struct('type','trace__');
trace__.dir_out_trace = dir_out_trace;
trace__.n_shuffle = n_shuffle;
trace__.trace_st__ = trace_st__;
trace__.niter_is__ = niter_is__;
trace__.r_rem_is__ = r_rem_is__;
trace__.c_rem_is__ = c_rem_is__;
trace__.QR_is__ = QR_is__;
trace__.QC_is__ = QC_is__;
trace__.I_rem_is__ = I_rem_is__;
trace__.QR_avg_i_ = QR_avg_i_;
trace__.QR_std_i_ = QR_std_i_;
trace__.ZR_is__ = ZR_is__;
trace__.nlpR_is__ = nlpR_is__;
trace__.QC_avg_i_ = QC_avg_i_;
trace__.QC_std_i_ = QC_std_i_;
trace__.ZC_is__ = ZC_is__;
trace__.nlpC_is__ = nlpC_is__;
%%%%;
trace__.niter_s0000_ = niter_s0000_;
trace__.r_rem_s0000_ = r_rem_s0000_;
trace__.c_rem_s0000_ = c_rem_s0000_;
trace__.QR_s0000_ = QR_s0000_;
trace__.QC_s0000_ = QC_s0000_;
trace__.I_rem_s0000_ = I_rem_s0000_;
trace__.trace_s0000_ = trace_s0000_;
%%%%;
ZR_s0000_ = bsxfun(@rdivide,bsxfun(@minus,QR_s0000_,QR_avg_i_),max(1e-12,QR_std_i_));
nlpR_s0000_ = -z_to_lp(ZR_s0000_);
ZC_s0000_ = bsxfun(@rdivide,bsxfun(@minus,QC_s0000_,QC_avg_i_),max(1e-12,QC_std_i_));
nlpC_s0000_ = -z_to_lp(ZC_s0000_);
trace__.ZR_s0000_ = ZR_s0000_;
trace__.nlpR_s0000_ = nlpR_s0000_;
trace__.ZC_s0000_ = ZC_s0000_;
trace__.nlpC_s0000_ = nlpC_s0000_;
%%%%;



 
