function ...
[ ...
 uAZ__ ...
,lim_use_ ...
,c_use__ ...
] = ...
imagesc_uAZ_0( ...
 u__ ...
,A__ ...
,Z__ ...
,lim_0in_ ...
,c_0in__ ...
,r_gap ...
,c_gap ...
);
na=0;
if (nargin<1+na); u__=[]; end; na=na+1;
if (nargin<1+na); A__=[]; end; na=na+1;
if (nargin<1+na); Z__=[]; end; na=na+1;
if (nargin<1+na); lim_0in_=[]; end; na=na+1;
if (nargin<1+na); c_0in__=[]; end; na=na+1;
if (nargin<1+na); r_gap=[]; end; na=na+1;
if (nargin<1+na); c_gap=[]; end; na=na+1;

if isempty(Z__); Z__ = zeros([0,size(A__,2)]); end;
flag_Z = ~isempty(Z__);
if isempty(u__); u__ = zeros([size(A__,1)+size(Z__,1),0]); end;
flag_u = ~isempty(u__);

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
if isempty(c_0in__); c_0in__ = colormap_80s; end;
if isempty(r_gap); r_gap = ceil((2/64)*(A_nrows + Z_nrows)); end;
if isempty(c_gap); c_gap = ceil((2/64)*A_ncols); end;
u_rep = ceil(c_gap/max(1,u_ncols));

uA__ = u__(1:A_nrows,:);
uA_rep__ = zeros(A_nrows,u_ncols*u_rep);
for nc=0:u_ncols-1; uA_rep__(:,1 + u_rep*nc + [0:u_rep-1]) = repmat(uA__(:,1+nc),[1,u_rep]); end;
uZ__ = u__(A_nrows + [1:Z_nrows],:);
uZ_rep__ = zeros(Z_nrows,u_ncols*u_rep);
for nc=0:u_ncols-1; uZ_rep__(:,1 + u_rep*nc + [0:u_rep-1]) = repmat(uZ__(:,1+nc),[1,u_rep]); end;

n_c = size(c_0in__,1);
dl = diff(lim_0in_)/(n_c-3);
val_up = min(lim_0in_) - 0.5*dl;
val_dn = min(lim_0in_) - 1.5*dl;
val_00 = min(lim_0in_) - 2.5*dl;
lim_use_ = [min(lim_0in_)-3*dl , max(lim_0in_)];

c_use__ = c_0in__;
c_use__(1+0,:) = [1,1,1]; %<-- val_00. ;
c_use__(1+1,:) = [0,1,0]; %<-- val_dn. ;
c_use__(1+2,:) = [1,1,0]; %<-- val_up. ;

uA_rep__ = val_up*(uA_rep__> 0) + val_dn*(uA_rep__<=0);
uZ_rep__ = val_up*(uZ_rep__> 0) + val_dn*(uZ_rep__<=0);

%%%%%%%%;
gn__ = val_00 * ones(A_nrows,c_gap) ;
gw__ = val_00 * ones(r_gap,u_ncols*u_rep) ;
gs__ = val_00 * ones(Z_nrows,c_gap) ;
ge__ = val_00 * ones(r_gap,A_ncols) ;
gc__ = val_00 * ones(r_gap,c_gap) ;
if ( flag_u &  flag_Z);
uAZ__ = [ ...
 uA_rep__ , gn__ ,  A__ ; ...
     gw__ , gc__ , ge__ ; ...
 uZ_rep__ , gs__ ,  Z__ ; ...
	];
end;%if ( flag_u &  flag_Z);
if ( flag_u & ~flag_Z);
uAZ__ = [ ...
 uA_rep__ , gn__ ,  A__ ; ...
	];
end;%if ( flag_u & ~flag_Z);
if (~flag_u &  flag_Z);
uAZ__ = [ ...
 A__ ; ...
 ge__ ; ...
 Z__ ; ...
	];
end;%if (~flag_u &  flag_Z);
if (~flag_u & ~flag_Z);
uAZ__ = [ ...
 A__ ; ...
	];
end;%if (~flag_u & ~flag_Z);
%%%%%%%%;

colormap(c_use__);
imagesc(uAZ__,lim_use_);
axisnotick;
