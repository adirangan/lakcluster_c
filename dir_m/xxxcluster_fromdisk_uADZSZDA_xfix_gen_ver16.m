function ...
[ ...
 parameter ...
,str_output ...
,gamma_d ...
] = ...
xxxcluster_fromdisk_uADZSZDA_xfix_gen_ver16( ...
 parameter ...
);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if isempty(parameter); parameter = struct('type','parameter'); end;
%%%%%%%%;
if ~isfield(parameter,'str_lak_vs_dex'); parameter.str_lak_vs_dex = 'dex'; end;
if ~isfield(parameter,'maf_lo_threshold'); parameter.maf_lo_threshold = 0.01; end;
if ~isfield(parameter,'maf_hi_threshold'); parameter.maf_hi_threshold = 0.5; end;
if ~isfield(parameter,'str_mr_0in'); parameter.str_mr_0in = ''; end;
if ~isfield(parameter,'str_mc_0in'); parameter.str_mc_0in = ''; end;
if ~isfield(parameter,'flag_reverse'); parameter.flag_reverse = 0; end;
if ~isfield(parameter,'n_mds'); parameter.n_mds = 2; end;
if ~isfield(parameter,'ij_mds_use_'); parameter.ij_mds_use_ = [1:2]; end;
if ~isfield(parameter,'n_mds_repl'); parameter.n_mds_repl = 0; end;
if ~isfield(parameter,'gamma'); parameter.gamma = 0; end;
if ~isfield(parameter,'B_MLT'); parameter.B_MLT = 0; end;
if ~isfield(parameter,'Ireq'); parameter.Ireq = 0; end;
if ~isfield(parameter,'n_scramble'); parameter.n_scramble = 0; end;
if ~isfield(parameter,'nshuffle'); parameter.nshuffle = 0; end;
str_lak_vs_dex = parameter.str_lak_vs_dex;
maf_lo_threshold = parameter.maf_lo_threshold;
maf_hi_threshold = parameter.maf_hi_threshold;
str_mr_0in = parameter.str_mr_0in;
str_mc_0in = parameter.str_mc_0in;
flag_reverse = parameter.flag_reverse;
n_mds = parameter.n_mds;
ij_mds_use_ = parameter.ij_mds_use_;
n_mds_repl = parameter.n_mds_repl;
gamma = parameter.gamma;
B_MLT = parameter.B_MLT;
Ireq = parameter.Ireq;
n_scramble = parameter.n_scramble;
nshuffle = parameter.nshuffle;
%%%%%%%%;

if nargin<1;
n_mds = 2; ij_mds_use_ = [1:2];
for str_lak_vs_dex_ = {'dex','lak'};
for maf_lo_threshold = [0.01,0.25];
for maf_hi_threshold = [0.15,0.50];
for str_mr_0in_ = {'','mrtest'};
for str_mc_0in_ = {'','mctest'};
for flag_reverse = [0:1];
for n_mds_repl = [0,3];
for gamma = [0.004 , 0.05];
for B_MLT = [24,34];
for Ireq = [0,5];
for n_scramble = [0,1];
for nshuffle = [0,8];
parameter = struct('type','parameter');
parameter.str_lak_vs_dex = str_lak_vs_dex_{1};
parameter.maf_lo_threshold = maf_lo_threshold;
parameter.maf_hi_threshold = maf_hi_threshold;
parameter.str_mr_0in = str_mr_0in_{1};
parameter.str_mc_0in = str_mc_0in_{1};
parameter.flag_reverse = flag_reverse;
parameter.n_mds = n_mds;
parameter.ij_mds_use_ = ij_mds_use_;
parameter.n_mds_repl = n_mds_repl;
parameter.gamma = gamma;
parameter.B_MLT = B_MLT;
parameter.Ireq = Ireq;
parameter.n_scramble = n_scramble;
parameter.nshuffle = nshuffle;
[~,str_output,~] = xxxcluster_fromdisk_uADZSZDA_xfix_gen_ver16(parameter);
disp(sprintf(' %% %s',str_output));
end;end;end;end;end;end;end;end;end;end;end;end;
disp('returning');return;
end;%if nargin<1;

if max(maf_lo_threshold,maf_hi_threshold)>=0.5; 
str_maf = sprintf('_p%.2d',floor(100*min(maf_lo_threshold,maf_hi_threshold)));
else; 
str_maf = sprintf('_p%.2dq%.2d',floor(100*min(maf_lo_threshold,maf_hi_threshold)),floor(100*max(maf_lo_threshold,maf_hi_threshold))); 
end;%if max(maf_lo_threshold,maf_hi_threshold)>=0.5; 
%%%%;
if (length(str_mr_0in)>0); str_mr = sprintf('_%s',str_mr_0in);
else; str_mr = ''; end;%if (length(str_mr_0in)>0);
%%%%;
if (length(str_mc_0in)>0); str_mc = sprintf('_%s',str_mc_0in);
else; str_mc = ''; end;%if (length(str_mc_0in)>0);
%%%%
if (flag_reverse==1); str_reverse = '_X'; else str_reverse = '_D'; end;
%%%%
mds_tmp_ = zeros(1,n_mds); 
mds_tmp_(ij_mds_use_)=1; 
mds_code = length(ij_mds_use_); 
%mds_code = dot(2.^[0:n_mds-1],mds_tmp_);
if (n_mds_repl<1 | length(ij_mds_use_)<2);
str_mds = sprintf('_m%d',mds_code);
else; str_mds = sprintf('_m%dr%d',mds_code,n_mds_repl); 
end;%if (n_mds_repl<1 | length(ij_mds_use_)<2);
%%%%;
gamma_d = floor(gamma*1000); str_gamma = sprintf('_g%.3d',gamma_d);
%%%%;
if ( (B_MLT<=0) | (B_MLT==32) ); str_B = ''; else; str_B = sprintf('_B%.2d',B_MLT); end;
%%%%;
if (Ireq<=0); str_Ireq = ''; else; str_Ireq = sprintf('_I%d',Ireq); end;
%%%%;
if (n_scramble<=0); str_scramble = ''; else; str_scramble = sprintf('_r%d',n_scramble); end;
%%%%;
if (nshuffle<=0); str_shuffle = ''; else; str_shuffle = sprintf('_s%.4d',nshuffle); end;
%%%%;
str_output = ...
sprintf( ...
'%s%s%s%s%s%s%s%s%s%s%s' ...
,str_lak_vs_dex ...
,str_maf ...
,str_mr ...
,str_mc ...
,str_reverse ...
,str_mds ...
,str_gamma ...
,str_B ...
,str_Ireq ...
,str_scramble ...
,str_shuffle ...
);

