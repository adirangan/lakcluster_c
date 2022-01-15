function ...
[ ...
 uAZm__ ...
,lim_use_ ...
,c_use__ ...
] = ...
imagesc_uAZm_1( ...
 u__ ...
,A__ ...
,Z__ ...
,mr__ ...
,mc__ ...
,lim_0in_ ...
,c_0in__ ...
,r_gap ...
,c_gap ...
);
na=0;
if (nargin<1+na); u__=[]; end; na=na+1;
if (nargin<1+na); A__=[]; end; na=na+1;
if (nargin<1+na); Z__=[]; end; na=na+1;
if (nargin<1+na); mr__=[]; end; na=na+1;
if (nargin<1+na); mc__=[]; end; na=na+1;
if (nargin<1+na); lim_0in_=[]; end; na=na+1;
if (nargin<1+na); c_0in__=[]; end; na=na+1;
if (nargin<1+na); r_gap=[]; end; na=na+1;
if (nargin<1+na); c_gap=[]; end; na=na+1;

if isempty(Z__); Z__ = zeros([0,size(A__,2)]); end;
flag_Z = ~isempty(Z__);
if isempty(u__); u__ = zeros([size(A__,1)+size(Z__,1),0]); end;
flag_u = ~isempty(u__);
if isempty(mr__); mr__ = zeros([size(A__,1)+size(Z__,1),0]); end;
flag_mr = ~isempty(mr__);
if isempty(mc__); mc__ = zeros([size(A__,2),0]); end;
flag_mc = ~isempty(mc__);

if isempty(lim_0in_);
min_A = min(A__,[],'all');
max_A = max(A__,[],'all');
min_Z = min(Z__,[],'all');
max_Z = max(Z__,[],'all');
if ( flag_Z); lim_0in_ = [ min(min_A,min_Z) , max(max_A,max_Z) ]; end;
if (~flag_Z); lim_0in_ = [ min_A , max_A ]; end;
end;%if isempty(lim_0in_);
A__ = min(A__,max(lim_0in_)); A__ = max(A__,min(lim_0in_));
Z__ = min(Z__,max(lim_0in_)); Z__ = max(Z__,min(lim_0in_));

[u_nrows,u_ncols] = size(u__);
[A_nrows,A_ncols] = size(A__);
[Z_nrows,Z_ncols] = size(Z__);
[mr_nrows,mr_ncols] = size(mr__);
[mc_nrows,mc_ncols] = size(mc__);
if isempty(c_0in__); c_0in__ = colormap_80s; end;
if isempty(r_gap); r_gap = ceil((2/64)*(A_nrows + Z_nrows)); end;
if isempty(c_gap); c_gap = ceil((2/64)*A_ncols); end;
u_rep = ceil(2*c_gap/max(1,u_ncols));
mr_rep = ceil(2*c_gap/max(1,mr_ncols));
mc_rep = ceil(2*r_gap/max(1,mc_ncols));

uA__ = u__(1:A_nrows,:);
uA_rep__ = zeros(A_nrows,u_ncols*u_rep);
for nc=0:u_ncols-1; uA_rep__(:,1 + u_rep*nc + [0:u_rep-1]) = repmat(uA__(:,1+nc),[1,u_rep]); end;
uZ__ = u__(A_nrows + [1:Z_nrows],:);
uZ_rep__ = zeros(Z_nrows,u_ncols*u_rep);
for nc=0:u_ncols-1; uZ_rep__(:,1 + u_rep*nc + [0:u_rep-1]) = repmat(uZ__(:,1+nc),[1,u_rep]); end;

mrA__ = mr__(1:A_nrows,:);
mrA_rep__ = zeros(A_nrows,mr_ncols*mr_rep);
for nc=0:mr_ncols-1; mrA_rep__(:,1 + mr_rep*nc + [0:mr_rep-1]) = repmat(mrA__(:,1+nc),[1,mr_rep]); end;
mrZ__ = mr__(A_nrows + [1:Z_nrows],:);
mrZ_rep__ = zeros(Z_nrows,mr_ncols*mr_rep);
for nc=0:mr_ncols-1; mrZ_rep__(:,1 + mr_rep*nc + [0:mr_rep-1]) = repmat(mrZ__(:,1+nc),[1,mr_rep]); end;
mc_rep__ = zeros(A_ncols,mc_ncols*mc_rep);
for nc=0:mc_ncols-1; mc_rep__(:,1 + mc_rep*nc + [0:mc_rep-1]) = repmat(mc__(:,1+nc),[1,mc_rep]); end;

n_c = size(c_0in__,1);
dl = diff(lim_0in_)/(n_c-5);
val_up = min(lim_0in_) - 0.5*dl;
val_un = min(lim_0in_) - 1.5*dl;
val_mp = min(lim_0in_) - 2.5*dl;
val_mn = min(lim_0in_) - 3.5*dl;
val_00 = min(lim_0in_) - 4.5*dl;
lim_use_ = [min(lim_0in_)-5*dl , max(lim_0in_)];

c_use__ = c_0in__;
c_use__(1+0,:) = [1,1,1]; %<-- val_00. ;
c_use__(1+3,:) = [0,1,0]; %<-- val_un. ;
c_use__(1+4,:) = [1,1,0]; %<-- val_up. ;
c_use__(1+1,:) = 0.85*[1,1,1]; %<-- val_mn. ;
c_use__(1+2,:) = 0.50*[1,1,1]; %<-- val_mp. ;

uA_rep__ = val_up*(uA_rep__> 0) + val_un*(uA_rep__<=0);
uZ_rep__ = val_up*(uZ_rep__> 0) + val_un*(uZ_rep__<=0);
mrA_rep__ = val_mp*(mrA_rep__> 0) + val_mn*(mrA_rep__<=0);
mrZ_rep__ = val_mp*(mrZ_rep__> 0) + val_mn*(mrZ_rep__<=0);
mc_rep__ = val_mp*(mc_rep__> 0) + val_mn*(mc_rep__<=0);
mc_rep__ = transpose(mc_rep__);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%     g00__ , g01__ , g02__     , g03__ , mc_rep__ ; ...
%     g10__ , g11__ , g12__     , g13__ , g14__    ; ...
%  uA_rep__ , g21__ , mrA_rep__ , g23__ , A__      ; ...
%     g30__ , g31__ , g32__     , g33__ , g34__    ; ...
%  uZ_rep__ , g41__ , mrZ_rep__ , g43__ , Z__      ; ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_r_g_ = [ ...
;flag_mc*mc_ncols*mc_rep ...
;flag_mc*r_gap ...
;A_nrows ...
;flag_Z*r_gap ...
;Z_nrows ...
];
n_c_g_ = [ ...
;flag_u*u_ncols*u_rep ...
;flag_u*c_gap ...
;flag_mr*mr_ncols*mr_rep ...
;flag_mr*c_gap ...
;A_ncols ...
];
g__ = cell(5,5);
for nr_g=0:5-1;
for nc_g=0:5-1;
g__{1+nr_g,1+nc_g} = val_00 * ones(n_r_g_(1+nr_g),n_c_g_(1+nc_g));
end;%for nc_g=0:5-1;
end;%for nr_g=0:5-1;
%%%%%%%%;
uAZm__ = [ ...
 g__{1+0,1+0} , g__{1+0,1+1} , g__{1+0,1+2} , g__{1+0,1+3} , mc_rep__     ; ...
 g__{1+1,1+0} , g__{1+1,1+1} , g__{1+1,1+2} , g__{1+1,1+3} , g__{1+1,1+4} ; ...
 uA_rep__     , g__{1+2,1+1} , mrA_rep__    , g__{1+2,1+3} , A__          ; ...
 g__{1+3,1+0} , g__{1+3,1+1} , g__{1+3,1+2} , g__{1+3,1+3} , g__{1+3,1+4} ; ...
 uZ_rep__     , g__{1+4,1+1} , mrZ_rep__    , g__{1+4,1+3} , Z__          ; ...
];
%%%%%%%%;

colormap(c_use__);
imagesc(uAZm__,lim_use_);
axisnotick;
