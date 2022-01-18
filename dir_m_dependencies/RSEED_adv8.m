function rseed = RSEED_adv8(rseed);
disp(sprintf(' %% [entering RSEED_adv8], %d',rseed));
for nr=0:8-1;
rseed = cast(RSEED_adv1(cast(rseed,'uint32')),'uint32');
disp(sprintf(' %% rseed %d',rseed));
end;%for nr=0:8-1;
disp(sprintf(' %% [finished RSEED_adv8], %d',rseed));

