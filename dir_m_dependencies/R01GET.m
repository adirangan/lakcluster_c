function [val,rseed] = R01GET(rseed);
% RNEXT = (RPREV*((unsigned long int)(pow(2,RPOW)+RADD))%((unsigned long int) pow(2,2*RPOW-1)));
POW2RPOWPLUSRADD = cast(35,'uint32');
POW22RPOWMINUSONE = cast(2147483647,'uint32');
rseed = RSEED_adv1(cast(rseed,'uint32'));
val = cast(rseed,'double')/cast(POW22RPOWMINUSONE,'double');

