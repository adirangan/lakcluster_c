function rseed = RSEED_adv1(rseed);
disp(sprintf(' %% Warning, does not loop correctly, do not use.'));
POW2RPOWPLUSRADD = cast(35,'uint32');
POW22RPOWMINUSONE = cast(2147483647,'uint32');
rseed_01 = mod(cast(rseed,'uint32'),POW22RPOWMINUSONE);
rseed_02 = mod(cast(rseed_01+rseed_01,'uint32'),POW22RPOWMINUSONE);
rseed_04 = mod(cast(rseed_02+rseed_02,'uint32'),POW22RPOWMINUSONE);
rseed_08 = mod(cast(rseed_04+rseed_04,'uint32'),POW22RPOWMINUSONE);
rseed_16 = mod(cast(rseed_08+rseed_08,'uint32'),POW22RPOWMINUSONE);
rseed_32 = mod(cast(rseed_16+rseed_16,'uint32'),POW22RPOWMINUSONE);
rseed_35 = mod(cast(rseed_32+rseed_02+rseed_01,'uint32'),POW22RPOWMINUSONE);
rseed = rseed_35;




