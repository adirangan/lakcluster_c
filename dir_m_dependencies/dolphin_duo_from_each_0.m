function ...
[ ...
 parameter ...
,a_duo__ ...
,A_duo___ ...
,BB_inv_duo___ ...
,CC_inv_duo___ ...
,index_nvar0_from_nvv_ ...
,index_nvar1_from_nvv_ ...
] = ...
dolphin_duo_from_each_0( ...
 parameter ...
,A_cmb__ ...
,BB_inv_cmb__ ...
,CC_inv_cmb__ ...
,a_each__ ...
,A_each__ ...
,BB_inv_each__ ...
,CC_inv_each__ ...
);

verbose=0;
if verbose; disp(sprintf(' %% [entering dolphin_duo_from_each_0]')); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na);,A_cmb__=[]; end; na=na+1;
if (nargin<1+na);,BB_inv_cmb__=[]; end; na=na+1;
if (nargin<1+na);,CC_inv_cmb__=[]; end; na=na+1;
if (nargin<1+na);,a_each__=[]; end; na=na+1;
if (nargin<1+na);,A_each__=[]; end; na=na+1;
if (nargin<1+na);,BB_inv_each__=[]; end; na=na+1;
if (nargin<1+na);,CC_inv_each__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;

n_var = size(A_cmb__,1);
n_duo = n_var*(n_var-1)/2;

index_nvar0_from_nvv_ = zeros(n_duo,1);
index_nvar1_from_nvv_ = zeros(n_duo,1);
a_duo__ = zeros(2,n_duo);
A_duo___ = zeros(2,2,n_duo);
BB_inv_duo___ = zeros(2,2,n_duo);
CC_inv_duo___ = zeros(2,2,n_duo);
na=0;
for nvar0=0:n_var-1;
if verbose; disp(sprintf(' %% nvar0 %d/%d',nvar0,n_var)); end;
for nvar1=nvar0+1:n_var-1;
index_nvar0_from_nvv_(1+na) = nvar0;
index_nvar1_from_nvv_(1+na) = nvar1;
%%%%;
a_pair_0 = a_each__(1+nvar0,1+nvar1);
a_pair_1 = a_each__(1+nvar1,1+nvar0);
a_pair_ = [ a_pair_0 ; a_pair_1 ];
%%%%;
A_pair_00 = A_each__(1+nvar0,1+nvar1);
A_pair_11 = A_each__(1+nvar1,1+nvar0);
A_pair_01 = A_cmb__(1+nvar0,1+nvar1);
A_pair_10 = A_cmb__(1+nvar1,1+nvar0);
A_pair__ = [ A_pair_00 , A_pair_01 ; A_pair_10 , A_pair_11 ];
%%%%;
BB_inv_pair_00 = BB_inv_each__(1+nvar0,1+nvar1);
BB_inv_pair_11 = BB_inv_each__(1+nvar1,1+nvar0);
BB_inv_pair_01 = BB_inv_cmb__(1+nvar0,1+nvar1);
BB_inv_pair_10 = BB_inv_cmb__(1+nvar1,1+nvar0);
BB_inv_pair__ = [ BB_inv_pair_00 , BB_inv_pair_01 ; BB_inv_pair_10 , BB_inv_pair_11 ];
%%%%;
CC_inv_pair_00 = CC_inv_each__(1+nvar0,1+nvar1);
CC_inv_pair_11 = CC_inv_each__(1+nvar1,1+nvar0);
CC_inv_pair_01 = CC_inv_cmb__(1+nvar0,1+nvar1);
CC_inv_pair_10 = CC_inv_cmb__(1+nvar1,1+nvar0);
CC_inv_pair__ = [ CC_inv_pair_00 , CC_inv_pair_01 ; CC_inv_pair_10 , CC_inv_pair_11 ];
%%%%;
a_duo__(:,1+na) = a_pair_;
A_duo___(:,:,1+na) = A_pair__;
BB_inv_duo___(:,:,1+na) = BB_inv_pair__;
CC_inv_duo___(:,:,1+na) = CC_inv_pair__;
%%%%;
na=na+1;
end;%for nvar1=nvar0+1:n_var-1;
end;%for nvar0=0:n_var-1;


if verbose; disp(sprintf(' %% [finished dolphin_duo_from_each_0]')); end;














  
