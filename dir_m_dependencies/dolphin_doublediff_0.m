function out_ = dolphin_doublediff_0(in0_,in1_);
na=0;
if (nargin<1+na); in0_=[]; end; na=na+1;
if (nargin<1+na); in1_=[]; end; na=na+1;

n_dim = ndims(in0_);
%%%%%%%%;
if (n_dim==1 | size(in0_,2)==1);
n_shuffle = size(in0_,1)-1; assert(size(in1_,1)==1+n_shuffle);
index_shuffle_ = [1:1+n_shuffle-1];
n_shuffle_use = numel(index_shuffle_);
out_ = bsxfun( ...
	       @minus ...
	       ,reshape(in0_(1+index_shuffle_),[1,n_shuffle_use]) ...
	       ,reshape(in1_(1+index_shuffle_),[n_shuffle_use,1]) ...
	       );
out_ = cat( 1 , in0_(1+0)-in1_(1+0) , reshape(out_,[n_shuffle^2,1]) );
end;%if (n_dim==1 | size(in0_,2)==1);
%%%%%%%%;
if (n_dim==2 & size(in0_,2)> 1);
n_var = size(in0_,1); assert(size(in1_,1)==n_var);
n_shuffle = size(in0_,2)-1; assert(size(in1_,2)==1+n_shuffle);
index_shuffle_ = [1:1+n_shuffle-1];
n_shuffle_use = numel(index_shuffle_);
out_ = bsxfun( ...
	       @minus ...
	       ,reshape(in0_(:,1+index_shuffle_),[n_var,1,n_shuffle_use]) ...
	       ,reshape(in1_(:,1+index_shuffle_),[n_var,n_shuffle_use,1]) ...
	       );
out_ = cat( 2 , in0_(:,1+0)-in1_(:,1+0) , reshape(out_,[n_var,n_shuffle^2]) );
end;%if (n_dim==2 & size(in0_,2)> 1);
%%%%%%%%;
if n_dim==3;
n_var = size(in0_,1); assert(size(in0_,2)==n_var); assert(size(in1_,1)==n_var); assert(size(in1_,2)==n_var);
n_shuffle = size(in0_,3)-1; assert(size(in1_,3)==1+n_shuffle);
index_shuffle_ = [1:1+n_shuffle-1];
n_shuffle_use = numel(index_shuffle_);
out_ = bsxfun( ...
	       @minus ...
	       ,reshape(in0_(:,:,1+index_shuffle_),[n_var,n_var,1,n_shuffle_use]) ...
	       ,reshape(in1_(:,:,1+index_shuffle_),[n_var,n_var,n_shuffle_use,1]) ...
	       );
out_ = cat( 3 , in0_(:,:,1+0)-in1_(:,:,1+0) , reshape(out_,[n_var,n_var,n_shuffle^2]) );
end%if n_dim==3;
%%%%%%%%;
