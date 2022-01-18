dt_avg = 0.257; %<-- matches d00 data. ;
relative_variation = 0.8973; %<-- matches d00 data. ;

for n_var = [2,3,4]; %<-- pair. ;
for n_T_0in = [128,256,512,1024]; %n_T_0in = 5664; %<-- matches d00 data. ;
snr_ = 2.^[-3:+0.5:+3]; n_snr = numel(snr_);
rseed_ = 0:16; n_rseed = numel(rseed_);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
for nsnr=0:n_snr-1;
snr = snr_(1+nsnr);
dolphin_test_5(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
end;%for nsnr=0:n_snr-1;
end;%for nrseed=0:n_rseed-1;
end;%for n_T_0in = [128,256,512,1024];
end;%for n_var = [2,3,4];
