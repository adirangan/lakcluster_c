function [Ireq,Icat,mds_used_,mds_repl] = lisa_set_mds_used(n_cov);
% set Ireq, Icat, mds_used_ and mds_repl. ;
if 0;
elseif n_cov==0; Ireq=0; Icat=1; mds_used_=[]; mds_repl=0; % do not correct for mds-components ;
elseif n_cov==1; Ireq=0; Icat=1; mds_used_=[1:2]; mds_repl=1; % only correct for mds-components [1:2], but replicated only once, so that each 'sector' is 90 degrees ;
elseif n_cov==2; Ireq=0; Icat=1; mds_used_=[1:2]; mds_repl=2; % only correct for mds-components [1:2], but replicated twice, so that each 'sector' is 45 degrees ;
elseif n_cov==3; Ireq=Ireq_half; Icat=Icat_full; mds_used_=[]; mds_repl=0; % correct for study, requiring at least Ireq_half, but not for mds-components. ;
elseif n_cov==4; Ireq=Ireq_half; Icat=Icat_full; mds_used_=[1:2]; mds_repl=1; % correct for study, requiring at least Ireq_half, as well as correct for mds-components [1:2], but replicated only once, so that each 'sector' is 90 degrees ;
elseif n_cov==5; Ireq=Ireq_half; Icat=Icat_full; mds_used_=[1:2]; mds_repl=2; % correct for study, requiring at least Ireq_half, as well as correct for mds-components [1:2], replicated twice, so that each 'sector' is 45 degrees ;
 else disp(sprintf(' %% Warning! n_cov %d undefined in lisa_set_mds_used.m',n_cov)); end;
