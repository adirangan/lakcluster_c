function bytes_written = mda_write_r8(input_,fname);

if nargin<1;
disp(sprintf(' %% testing mda_write_r8;'));
rng(0);
input_1_ = randn(3,4,5);
fname = 'test_mda_write_r8.mda';
mda_write_r8(input_1_,fname);
input_2_ = mda_read_r8(fname);
disp(sprintf(' %% input_1_: '));
disp(num2str(transpose(input_1_(:))));
disp(sprintf(' %% input_2_: '));
disp(num2str(transpose(input_2_(:))));
disp(sprintf(' %% input_1_ vs input_2_: %0.16f',fnorm(input_1_-input_2_)/max(1e-12,fnorm(input_1_))));
disp('returning'); return;
end;%if nargin<1;

% write mda file readable by mda_read_r8.m ;
d_ = size(input_);
n_d = numel(d_);
dd = prod(d_);
fp = fopen(fname,'w');
if (fp==-1); disp(sprintf(' %% Warning! cannot open %s in mda_write_r8',fname)); end;
bytes_written = 0;
bytes_written = bytes_written + fwrite(fp,n_d,'int32');
bytes_written = bytes_written + fwrite(fp,d_,'int32');
bytes_written = bytes_written + fwrite(fp,reshape(input_,[dd,1]),'double');
fclose(fp);

