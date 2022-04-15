function ...
[ ...
 mx__ ...
] = ...
load_mx__from_parameter_ver0( ...
 parameter ...
);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if isempty(parameter); parameter = struct('type','parameter'); end;
%%%%%%%%;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'dir_0in'); parameter.dir_0in = pwd; end;
if ~isfield(parameter,'str_prefix'); parameter.str_prefix = 'test'; end;
if ~isfield(parameter,'str_mr_0in'); parameter.str_mr_0in = ''; end;
if ~isfield(parameter,'str_mc_0in'); parameter.str_mc_0in = ''; end;
if ~isfield(parameter,'flag_reverse'); parameter.flag_reverse = 0; end;
if ~isfield(parameter,'n_study'); parameter.n_study = 1; end;
if ~isfield(parameter,'n_bin'); parameter.n_bin = parameter.n_study; end;

n_study = parameter.n_study;

mx__ = struct('type','mx__');
mx__.flag_verbose = parameter.flag_verbose;
mx__.flag_reverse = parameter.flag_reverse;
mx__.n_study = parameter.n_study;
mx__.dir_0in = parameter.dir_0in;
mx__.str_prefix = parameter.str_prefix;
mx__.str_mr_0in = parameter.str_mr_0in;
mx__.str_mc_0in = parameter.str_mc_0in;

mx__.mr_A_default__ = cell(n_study,1);
mx__.mr_Z_default__ = cell(n_study,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% check to make sure that the mr_A_full_ lines up with the mr_A_%0.2d_ ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
mx__.mr_A_default__ = cell(mx__.n_study,1);
mx__.mr_Z_default__ = cell(mx__.n_study,1);
for ns=0:mx__.n_study-1;
tmp_fname = sprintf('%s/%s_mr_A_%0.2d.b16',mx__.dir_0in,mx__.str_prefix,1+ns); fcheck(tmp_fname);
mx__.mr_A_default__{1+ns} = (binary_uncompress(tmp_fname)>0);
tmp_fname = sprintf('%s/%s_mr_Z_%0.2d.b16',mx__.dir_0in,mx__.str_prefix,1+ns); fcheck(tmp_fname);
mx__.mr_Z_default__{1+ns} = (binary_uncompress(tmp_fname)>0);
end;%for ns=0:mx__.n_study-1;
tmp_fname = sprintf('%s/%s_mr_A_full.b16',mx__.dir_0in,mx__.str_prefix); fcheck(tmp_fname);
mx__.mr_A_default_full_ = (binary_uncompress(tmp_fname)>0);
tmp_fname = sprintf('%s/%s_mr_Z_full.b16',mx__.dir_0in,mx__.str_prefix); fcheck(tmp_fname);
mx__.mr_Z_default_full_ = (binary_uncompress(tmp_fname)>0);
mr_A_tmp_ = []; mr_Z_tmp_ = [];
for ns=0:mx__.n_study-1;
mr_A_tmp_ = [mr_A_tmp_;mx__.mr_A_default__{1+ns}];
mr_Z_tmp_ = [mr_Z_tmp_;mx__.mr_Z_default__{1+ns}];
end;% for ns=0:mx__.n_study-1;
%%%%%%%%;
tmp_fname = sprintf('%s/%s_mc_A.b16',mx__.dir_0in,mx__.str_prefix);
mx__.mc_A_default_ = (binary_uncompress(tmp_fname)>0);
%%%%%%%%;
if (mx__.flag_verbose>0); disp(sprintf(' %% default masks: mr_A error %0.16f ; mr_Z error %0.16f ; mc_A %d',fnorm(mr_A_tmp_-mx__.mr_A_default_full_),fnorm(mr_Z_tmp_-mx__.mr_Z_default_full_),sum(mx__.mc_A_default_))); end;
%%%%%%%%;
if length(mx__.str_mr_0in)>0;
mx__.mr_A__ = cell(mx__.n_study,1);
mx__.mr_Z__ = cell(mx__.n_study,1);
for ns=0:mx__.n_study-1;
tmp_fname = sprintf('%s/%s_mr_A_%s_%0.2d.b16',mx__.dir_0in,mx__.str_prefix,mx__.str_mr_0in,1+ns); fcheck(tmp_fname);
mx__.mr_A__{1+ns} = (binary_uncompress(tmp_fname)>0);
tmp_fname = sprintf('%s/%s_mr_Z_%s_%0.2d.b16',mx__.dir_0in,mx__.str_prefix,mx__.str_mr_0in,1+ns); fcheck(tmp_fname);
mx__.mr_Z__{1+ns} = (binary_uncompress(tmp_fname)>0);
end;%for ns=0:mx__.n_study-1;
tmp_fname = sprintf('%s/%s_mr_A_%s_full.b16',mx__.dir_0in,mx__.str_prefix,mx__.str_mr_0in); fcheck(tmp_fname);
mx__.mr_A_full_ = (binary_uncompress(tmp_fname)>0);
tmp_fname = sprintf('%s/%s_mr_Z_%s_full.b16',mx__.dir_0in,mx__.str_prefix,mx__.str_mr_0in); fcheck(tmp_fname);
mx__.mr_Z_full_ = (binary_uncompress(tmp_fname)>0);
mr_A_tmp_ = []; mr_Z_tmp_ = [];
for ns=0:mx__.n_study-1;
mr_A_tmp_ = [mr_A_tmp_;mx__.mr_A__{1+ns}];
mr_Z_tmp_ = [mr_Z_tmp_;mx__.mr_Z__{1+ns}];
end;% for ns=0:mx__.n_study-1;
if (mx__.flag_verbose>-1); disp(sprintf(' %% %s masks: mr_A error %0.16f ; mr_Z error %0.16f',mx__.str_mr_0in,fnorm(mr_A_tmp_-mx__.mr_A_full_),fnorm(mr_Z_tmp_-mx__.mr_Z_full_))); end;
 else;
mx__.mr_A__ = mx__.mr_A_default__;
mx__.mr_Z__ = mx__.mr_Z_default__;
mx__.mr_A_full_ = mx__.mr_A_default_full_;
mx__.mr_Z_full_ = mx__.mr_Z_default_full_;
end;%if length(mx__.str_mr_0in)>0;
%%%%%%%%;
if length(mx__.str_mc_0in)>0;
tmp_fname = sprintf('%s/%s_mc_A_%s.b16',mx__.dir_0in,mx__.str_prefix,mx__.str_mc_0in);
mx__.mc_A_ = (binary_uncompress(tmp_fname)>0);
if (mx__.flag_verbose>-1); disp(sprintf(' %% %s masks: mc_A %d',mx__.str_mr_0in,sum(mx__.mc_A_))); end;
 else;
mx__.mc_A_ = mx__.mc_A_default_;
end;%if length(mx__.str_mc_0in)>0;
%%%%%%%%;
if (mx__.flag_reverse==0); mx__.mr_D_ = mx__.mr_A_full_; mx__.mr_X_ = mx__.mr_Z_full_; end;%if (mx__.flag_reverse==0);
if (mx__.flag_reverse==1); mx__.mr_D_ = mx__.mr_Z_full_; mx__.mr_X_ = mx__.mr_A_full_; end;%if (mx__.flag_reverse==1);
%%%%%%%%;
tmp_fname = sprintf('%s/%s_mc_T.b16',mx__.dir_0in,mx__.str_prefix);
mx__.mc_T_ = (binary_uncompress(tmp_fname)>0);
mx__.mc_T_default_ = (binary_uncompress(tmp_fname)>0);


 
