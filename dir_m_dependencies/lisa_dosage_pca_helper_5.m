function ...
[...
 tmp_dosage_imputed_row_AnV_,tmp_dosage_rounded_row_AnV_,tmp_dosage_genotyp_row_AnV_...
,tmp_dosage_genotyp_row_and_,tmp_dosage_genotyp_row_xor_,tmp_dosage_genotyp_row_nor_...
,tmp_dosage_genotyp_row_and_normalized_,tmp_dosage_genotyp_row_xor_normalized_,tmp_dosage_genotyp_row_nor_normalized_...
,mismatch_d1_imputed_,mismatch_d1_genotyp_...
,mismatch_d2_imputed_,mismatch_d2_genotyp_...
 ] = lisa_dosage_pca_helper_5(...
 tmp_dosage_row_1_,tmp_dosage_row_2_...
,n_dosage_snp,nstudy_arm2...
,mismatch_d1_genotyp__,mismatch_d2_genotyp__....
,bo1_,bu1_,bo2_,bu2_,dud_,eu1_,eu2_...
,allele1_dosage_,allele2_dosage_...
,allele1_arm1_,allele2_arm1_,alleletype_arm1_...
,allele1_arm2_,allele2_arm2_,alleletype_arm2_...
,bu3_by_bu1_xref_,bu2_by_bu3_xref_...
,eu1_by_dud_xref_,bu1_by_eu1_xref_,bo1_by_bu1_xref_...
,eu2_by_dud_xref_,bu2_by_eu2_xref_,bo2_by_bu2_xref_...
,atfr_arm1_perm_and_,atfr_arm1_perm_and_alpha_,atfr_arm1_perm_and_delta_...
,atfr_arm1_perm_xor_,atfr_arm1_perm_xor_alpha_,atfr_arm1_perm_xor_delta_...
,atfr_arm1_perm_nor_,atfr_arm1_perm_nor_alpha_,atfr_arm1_perm_nor_delta_...
,A_p_arm1_perm_and_,A_p_arm1_perm_and_alpha_,A_p_arm1_perm_and_delta_...
,A_p_arm1_perm_xor_,A_p_arm1_perm_xor_alpha_,A_p_arm1_perm_xor_delta_...
,A_p_arm1_perm_nor_,A_p_arm1_perm_nor_alpha_,A_p_arm1_perm_nor_delta_...
,atfr_arm2_perm_and_,atfr_arm2_perm_and_alpha_,atfr_arm2_perm_and_delta_...
,atfr_arm2_perm_xor_,atfr_arm2_perm_xor_alpha_,atfr_arm2_perm_xor_delta_...
,atfr_arm2_perm_nor_,atfr_arm2_perm_nor_alpha_,atfr_arm2_perm_nor_delta_...
,A_p_arm2_perm_and_,A_p_arm2_perm_and_alpha_,A_p_arm2_perm_and_delta_...
,A_p_arm2_perm_xor_,A_p_arm2_perm_xor_alpha_,A_p_arm2_perm_xor_delta_...
,A_p_arm2_perm_nor_,A_p_arm2_perm_nor_alpha_,A_p_arm2_perm_nor_delta_...
,V_rx_trnx_perm_and_,V_rx_trnx_perm_xor_,V_rx_trnx_perm_nor_...
,V_rx_tstx_perm_and_,V_rx_tstx_perm_xor_,V_rx_tstx_perm_nor_...
);

error_code_null = 0;      % <-- No error; alleles match and minor-allele-frequency is commensurate. ;
error_code_synonym = 1;   % <-- Synonym; allele names must be interchanged with synonyms (i.e,. A<-->T and C<-->G). ;
error_code_dominance = 2; % <-- Dominance; the alleles must be switched (i.e., XY --> YX) so that the first allele listed corresponds to the minor (i.e., less frequent) allele. ;
error_code_both = 3;      % <-- Both; the allele names must be interchanged *and* the order must be switched. ;
error_code_ambiguous = 8; % <-- Ambiguity; I cannot determine the correct alignment. ;
error_code_missing = 9;   % <-- Missing; this snp is not in the overlap between the two sets. ;

%%%%%%%%%%%%%%%%;
% Now step through the dosage file line by line. ;
% updating the inner-products as we go. ;
%%%%%%%%%%%%%%%%;
ent_cutoff = 0.03; 
flag_V_trnx = 0; %<-- use V from training set (all snps from arm1). ;
flag_V_tstx = 1; %<-- use V from testing set (only snps in intersection of arm1 and arm2). ;
flag_normalize_atfr = 0; %<-- use atfr to normalize. ;
flag_normalize_A_p = 1; %<-- use A_p to normalize. ;
assert(~(flag_normalize_atfr & flag_normalize_A_p));
flag_normalize_arm1 = 0; %<-- use arm1 to normalize. ;
flag_normalize_arm2 = 1; %<-- use arm2 to normalize. ;
assert(~(flag_normalize_arm1 & flag_normalize_arm2));
flag_match_arm1 = 0; %<-- swap dosage data to match arm1. ;
flag_match_arm2 = 1; %<-- swap dosage data to match arm2. ;
assert(~(flag_match_arm1 & flag_match_arm2));
tmp_dosage_imputed_row_AnV_ = zeros(1,2);
tmp_dosage_rounded_row_AnV_ = zeros(1,2);
tmp_dosage_genotyp_row_AnV_ = zeros(1,2);
tmp_dosage_genotyp_row_and_ = zeros(n_dosage_snp,1);
tmp_dosage_genotyp_row_xor_ = zeros(n_dosage_snp,1);
tmp_dosage_genotyp_row_nor_ = zeros(n_dosage_snp,1);
tmp_dosage_genotyp_row_and_normalized_ = zeros(n_dosage_snp,1);
tmp_dosage_genotyp_row_xor_normalized_ = zeros(n_dosage_snp,1);
tmp_dosage_genotyp_row_nor_normalized_ = zeros(n_dosage_snp,1);
mismatch_d1_imputed_ = 9*ones(n_dosage_snp,1);
mismatch_d1_rounded_ = 9*ones(n_dosage_snp,1);
mismatch_d1_genotyp_ = 9*ones(n_dosage_snp,1);
mismatch_d2_imputed_ = 9*ones(n_dosage_snp,1);
mismatch_d2_rounded_ = 9*ones(n_dosage_snp,1);
mismatch_d2_genotyp_ = 9*ones(n_dosage_snp,1);
flag_bo1_all_sum = 0; flag_bo2_and_sum = 0; flag_bo2_xor_sum = 0; flag_bo2_nor_sum = 0; n_z_a = [0,0];
%%%%%%%%;
for ndosage_snp=1:n_dosage_snp;
if (mod(ndosage_snp,1000)==0); disp(sprintf(' %% ndosage_snp %d/%d',ndosage_snp,n_dosage_snp)); end;
ndud = ndosage_snp;
%%%%%%%%;
% check arm1 for snp;
%%%%%%%%;
neu1 = find(sum(eu1_by_dud_xref_(:,ndud),2)); assert(length(neu1)<=1);
if (length(neu1)==1);
flag_bo1_all_sum = flag_bo1_all_sum + 1;
nbu1 = find(sum(bu1_by_eu1_xref_(:,neu1),2)); assert(length(nbu1)==1);
tmp_bo1_ = find(sum(bo1_by_bu1_xref_(:,nbu1),2)); assert(length(tmp_bo1_)>=1);
nbo1 = tmp_bo1_(1);
assert(strcmp(bo1_{nbo1},dud_{ndud}));
assert(strcmp(bo1_{nbo1},eu1_{neu1}));
assert(strcmp(bo1_{nbo1},bu1_{nbu1}));
tmp_dosage_imputed_val_ = [tmp_dosage_row_1_(ndud),tmp_dosage_row_2_(ndud)];
tmp_dosage_rounded_val_ = max(0,min(2,round(tmp_dosage_imputed_val_)));
tmp_dosage_imputed_val_nor_ = tmp_dosage_imputed_val_(1:2:end-1);
tmp_dosage_imputed_val_xor_ = tmp_dosage_imputed_val_(2:2:end-0);
tmp_dosage_imputed_val_and_ = ones(1,1) - tmp_dosage_imputed_val_nor_ - tmp_dosage_imputed_val_xor_;
tmp_dosage_rounded_val_nor_ = tmp_dosage_rounded_val_(1:2:end-1);
tmp_dosage_rounded_val_xor_ = tmp_dosage_rounded_val_(2:2:end-0);
tmp_dosage_rounded_val_and_ = ones(1,1) - tmp_dosage_rounded_val_nor_ - tmp_dosage_rounded_val_xor_;
%%%%%%%%;
% check arm2 for snp;
%%%%%%%%;
tmp_error_code_d2 = error_code_missing;
flag_bo2_and = 0; flag_bo2_xor = 0; flag_bo2_nor = 0;
nbu3 = find(sum(bu3_by_bu1_xref_(:,nbu1),2)); assert(length(nbu3)<=1);
if (~isempty(nbu3));
nbu2 = find(sum(bu2_by_bu3_xref_(:,nbu3),2)); assert(length(nbu2)==1);
tmp_bo2_ = find(sum(bo2_by_bu2_xref_(:,nbu2),2)); assert(length(tmp_bo2_)>=1);
%%%%%%%%;
% check for mismatch in d2. ;
%%%%%%%%;
nbo2 = tmp_bo2_(1);
assert(strcmp(bo2_{nbo2},dud_{ndud}));
assert(strcmp(bo2_{nbo2},bu1_{nbu1}));
tmp_error_code_d2 = lisa_allele_mismatch_ver0(...
  allele1_dosage_(ndud),allele2_dosage_(ndud),allele1_arm2_(nbo2),allele2_arm2_(nbo2)...
  );
tmp_error_code_d2 = mismatch_d2_genotyp__{nstudy_arm2}(ndud);
mismatch_d2_genotyp_(ndud) = tmp_error_code_d2;
mismatch_d2_imputed_(ndud) = tmp_error_code_d2;
%%%%%%%%;
for nlo2=1:length(tmp_bo2_);
nbo2 = tmp_bo2_(nlo2);
if strcmp(alleletype_arm2_{nbo2},'and'); flag_bo2_and = 1; flag_bo2_and_sum = flag_bo2_and_sum + 1; end;
if strcmp(alleletype_arm2_{nbo2},'xor'); flag_bo2_xor = 1; flag_bo2_xor_sum = flag_bo2_xor_sum + 1; end;
if strcmp(alleletype_arm2_{nbo2},'nor'); flag_bo2_nor = 1; flag_bo2_nor_sum = flag_bo2_nor_sum + 1; end;
end;%for nlo2=1:length(tmp_bo2_);
end;%if (~isempty(nbu3));
tmp_dosage_genotyp_val_nor_ = flag_bo2_nor * tmp_dosage_rounded_val_nor_;
tmp_dosage_genotyp_val_xor_ = flag_bo2_xor * tmp_dosage_rounded_val_xor_;
tmp_dosage_genotyp_val_and_ = flag_bo2_and * tmp_dosage_rounded_val_and_;
%%%%%%%%;
% check for mismatch in d1. ;
%%%%%%%%;
tmp_error_code_d1 = lisa_allele_mismatch_ver0(...
  allele1_dosage_(ndud),allele2_dosage_(ndud),allele1_arm1_(nbo1),allele2_arm1_(nbo1)...
  );
tmp_error_code_d1 = mismatch_d1_genotyp__{nstudy_arm2}(ndud);
mismatch_d1_genotyp_(ndud) = tmp_error_code_d1;
if (flag_bo2_and || flag_bo2_xor || flag_bo2_nor); mismatch_d1_imputed_(ndud) = tmp_error_code_d1; end;
%%%%%%%%;
% match to either arm1 or arm2 if necessary. ;
%%%%%%%%;
if ( flag_match_arm1 & ( (tmp_error_code_d1==error_code_dominance) | (tmp_error_code_d1==error_code_both) ) )...
 | ( flag_match_arm2 & ( (tmp_error_code_d2==error_code_dominance) | (tmp_error_code_d2==error_code_both) ) );
tmp_ = tmp_dosage_imputed_val_nor_ ; tmp_dosage_imputed_val_nor_ = tmp_dosage_imputed_val_and_; tmp_dosage_imputed_val_and_ = tmp_ ; %<-- switching dominant allele ;
tmp_ = tmp_dosage_rounded_val_nor_ ; tmp_dosage_rounded_val_nor_ = tmp_dosage_rounded_val_and_; tmp_dosage_rounded_val_and_ = tmp_ ; %<-- switching dominant allele ;
tmp_ = tmp_dosage_genotyp_val_nor_ ; tmp_dosage_genotyp_val_nor_ = tmp_dosage_genotyp_val_and_; tmp_dosage_genotyp_val_and_ = tmp_ ; %<-- switching dominant allele ;
end;%if match. ;
%%%%%%%%;
tmp_dosage_genotyp_row_and_(ndud) = tmp_dosage_genotyp_val_and_;
tmp_dosage_genotyp_row_xor_(ndud) = tmp_dosage_genotyp_val_xor_;
tmp_dosage_genotyp_row_nor_(ndud) = tmp_dosage_genotyp_val_nor_;
%%%%%%%%;
% set to zero if necessary. ;
%%%%%%%%;
if ( flag_match_arm1 & ( 0*(tmp_error_code_d1==error_code_ambiguous) | 1*(tmp_error_code_d1==error_code_missing) ) ) ...
 | ( flag_match_arm2 & ( 0*(tmp_error_code_d2==error_code_ambiguous) | 1*(tmp_error_code_d2==error_code_missing) ) ) ;
tmp_dosage_imputed_val_and_normalized_ = zeros(1,1);
tmp_dosage_imputed_val_xor_normalized_ = zeros(1,1);
tmp_dosage_imputed_val_nor_normalized_ = zeros(1,1);
tmp_dosage_rounded_val_and_normalized_ = zeros(1,1);
tmp_dosage_rounded_val_xor_normalized_ = zeros(1,1);
tmp_dosage_rounded_val_nor_normalized_ = zeros(1,1);
tmp_dosage_genotyp_val_and_normalized_ = zeros(1,1);
tmp_dosage_genotyp_val_xor_normalized_ = zeros(1,1);
tmp_dosage_genotyp_val_nor_normalized_ = zeros(1,1);
%%%%%%%%;
% otherwise:
%%%%%%%%;
else;
if flag_normalize_atfr & flag_normalize_arm1
assert(isfinite(atfr_arm1_perm_and_alpha_(ndud)) & atfr_arm1_perm_and_delta_(ndud));
tmp_dosage_imputed_val_and_normalized_ = (2*tmp_dosage_imputed_val_and_ - 1 - atfr_arm1_perm_and_alpha_(ndud)) * sqrt(atfr_arm1_perm_and_delta_(ndud));
tmp_dosage_imputed_val_xor_normalized_ = (2*tmp_dosage_imputed_val_xor_ - 1 - atfr_arm1_perm_xor_alpha_(ndud)) * sqrt(atfr_arm1_perm_xor_delta_(ndud));
tmp_dosage_imputed_val_nor_normalized_ = (2*tmp_dosage_imputed_val_nor_ - 1 - atfr_arm1_perm_nor_alpha_(ndud)) * sqrt(atfr_arm1_perm_nor_delta_(ndud));
tmp_dosage_rounded_val_and_normalized_ = (2*tmp_dosage_rounded_val_and_ - 1 - atfr_arm1_perm_and_alpha_(ndud)) * sqrt(atfr_arm1_perm_and_delta_(ndud));
tmp_dosage_rounded_val_xor_normalized_ = (2*tmp_dosage_rounded_val_xor_ - 1 - atfr_arm1_perm_xor_alpha_(ndud)) * sqrt(atfr_arm1_perm_xor_delta_(ndud));
tmp_dosage_rounded_val_nor_normalized_ = (2*tmp_dosage_rounded_val_nor_ - 1 - atfr_arm1_perm_nor_alpha_(ndud)) * sqrt(atfr_arm1_perm_nor_delta_(ndud));
tmp_dosage_genotyp_val_and_normalized_ = (2*tmp_dosage_genotyp_val_and_ - 1 - atfr_arm1_perm_and_alpha_(ndud)) * sqrt(atfr_arm1_perm_and_delta_(ndud));
tmp_dosage_genotyp_val_xor_normalized_ = (2*tmp_dosage_genotyp_val_xor_ - 1 - atfr_arm1_perm_xor_alpha_(ndud)) * sqrt(atfr_arm1_perm_xor_delta_(ndud));
tmp_dosage_genotyp_val_nor_normalized_ = (2*tmp_dosage_genotyp_val_nor_ - 1 - atfr_arm1_perm_nor_alpha_(ndud)) * sqrt(atfr_arm1_perm_nor_delta_(ndud));
end;%if flag_normalize_atfr & flag_normalize_arm1
if flag_normalize_A_p & flag_normalize_arm1
assert(isfinite(A_p_arm1_perm_and_alpha_(ndud)) & A_p_arm1_perm_and_delta_(ndud));
tmp_dosage_imputed_val_and_normalized_ = (2*tmp_dosage_imputed_val_and_ - 1 - A_p_arm1_perm_and_alpha_(ndud)) * sqrt(A_p_arm1_perm_and_delta_(ndud));
tmp_dosage_imputed_val_xor_normalized_ = (2*tmp_dosage_imputed_val_xor_ - 1 - A_p_arm1_perm_xor_alpha_(ndud)) * sqrt(A_p_arm1_perm_xor_delta_(ndud));
tmp_dosage_imputed_val_nor_normalized_ = (2*tmp_dosage_imputed_val_nor_ - 1 - A_p_arm1_perm_nor_alpha_(ndud)) * sqrt(A_p_arm1_perm_nor_delta_(ndud));
tmp_dosage_rounded_val_and_normalized_ = (2*tmp_dosage_rounded_val_and_ - 1 - A_p_arm1_perm_and_alpha_(ndud)) * sqrt(A_p_arm1_perm_and_delta_(ndud));
tmp_dosage_rounded_val_xor_normalized_ = (2*tmp_dosage_rounded_val_xor_ - 1 - A_p_arm1_perm_xor_alpha_(ndud)) * sqrt(A_p_arm1_perm_xor_delta_(ndud));
tmp_dosage_rounded_val_nor_normalized_ = (2*tmp_dosage_rounded_val_nor_ - 1 - A_p_arm1_perm_nor_alpha_(ndud)) * sqrt(A_p_arm1_perm_nor_delta_(ndud));
tmp_dosage_genotyp_val_and_normalized_ = (2*tmp_dosage_genotyp_val_and_ - 1 - A_p_arm1_perm_and_alpha_(ndud)) * sqrt(A_p_arm1_perm_and_delta_(ndud));
tmp_dosage_genotyp_val_xor_normalized_ = (2*tmp_dosage_genotyp_val_xor_ - 1 - A_p_arm1_perm_xor_alpha_(ndud)) * sqrt(A_p_arm1_perm_xor_delta_(ndud));
tmp_dosage_genotyp_val_nor_normalized_ = (2*tmp_dosage_genotyp_val_nor_ - 1 - A_p_arm1_perm_nor_alpha_(ndud)) * sqrt(A_p_arm1_perm_nor_delta_(ndud));
end;%if flag_normalize_A_p & flag_normalize_arm1
%%%%%%%%;
if flag_normalize_atfr & flag_normalize_arm2;
assert(isfinite(atfr_arm2_perm_and_alpha_(ndud)) & atfr_arm2_perm_and_delta_(ndud));
tmp_dosage_imputed_val_and_normalized_ = (2*tmp_dosage_imputed_val_and_ - 1 - atfr_arm2_perm_and_alpha_(ndud)) * sqrt(atfr_arm2_perm_and_delta_(ndud));
tmp_dosage_imputed_val_xor_normalized_ = (2*tmp_dosage_imputed_val_xor_ - 1 - atfr_arm2_perm_xor_alpha_(ndud)) * sqrt(atfr_arm2_perm_xor_delta_(ndud));
tmp_dosage_imputed_val_nor_normalized_ = (2*tmp_dosage_imputed_val_nor_ - 1 - atfr_arm2_perm_nor_alpha_(ndud)) * sqrt(atfr_arm2_perm_nor_delta_(ndud));
tmp_dosage_rounded_val_and_normalized_ = (2*tmp_dosage_rounded_val_and_ - 1 - atfr_arm2_perm_and_alpha_(ndud)) * sqrt(atfr_arm2_perm_and_delta_(ndud));
tmp_dosage_rounded_val_xor_normalized_ = (2*tmp_dosage_rounded_val_xor_ - 1 - atfr_arm2_perm_xor_alpha_(ndud)) * sqrt(atfr_arm2_perm_xor_delta_(ndud));
tmp_dosage_rounded_val_nor_normalized_ = (2*tmp_dosage_rounded_val_nor_ - 1 - atfr_arm2_perm_nor_alpha_(ndud)) * sqrt(atfr_arm2_perm_nor_delta_(ndud));
tmp_dosage_genotyp_val_and_normalized_ = (2*tmp_dosage_genotyp_val_and_ - 1 - atfr_arm2_perm_and_alpha_(ndud)) * sqrt(atfr_arm2_perm_and_delta_(ndud));
tmp_dosage_genotyp_val_xor_normalized_ = (2*tmp_dosage_genotyp_val_xor_ - 1 - atfr_arm2_perm_xor_alpha_(ndud)) * sqrt(atfr_arm2_perm_xor_delta_(ndud));
tmp_dosage_genotyp_val_nor_normalized_ = (2*tmp_dosage_genotyp_val_nor_ - 1 - atfr_arm2_perm_nor_alpha_(ndud)) * sqrt(atfr_arm2_perm_nor_delta_(ndud));
end;%if flag_normalize_atfr & flag_normalize_arm2
if flag_normalize_A_p & flag_normalize_arm2;
assert(isfinite(A_p_arm2_perm_and_alpha_(ndud)) & A_p_arm2_perm_and_delta_(ndud));
tmp_dosage_imputed_val_and_normalized_ = (2*tmp_dosage_imputed_val_and_ - 1 - A_p_arm2_perm_and_alpha_(ndud)) * sqrt(A_p_arm2_perm_and_delta_(ndud));
tmp_dosage_imputed_val_xor_normalized_ = (2*tmp_dosage_imputed_val_xor_ - 1 - A_p_arm2_perm_xor_alpha_(ndud)) * sqrt(A_p_arm2_perm_xor_delta_(ndud));
tmp_dosage_imputed_val_nor_normalized_ = (2*tmp_dosage_imputed_val_nor_ - 1 - A_p_arm2_perm_nor_alpha_(ndud)) * sqrt(A_p_arm2_perm_nor_delta_(ndud));
tmp_dosage_rounded_val_and_normalized_ = (2*tmp_dosage_rounded_val_and_ - 1 - A_p_arm2_perm_and_alpha_(ndud)) * sqrt(A_p_arm2_perm_and_delta_(ndud));
tmp_dosage_rounded_val_xor_normalized_ = (2*tmp_dosage_rounded_val_xor_ - 1 - A_p_arm2_perm_xor_alpha_(ndud)) * sqrt(A_p_arm2_perm_xor_delta_(ndud));
tmp_dosage_rounded_val_nor_normalized_ = (2*tmp_dosage_rounded_val_nor_ - 1 - A_p_arm2_perm_nor_alpha_(ndud)) * sqrt(A_p_arm2_perm_nor_delta_(ndud));
tmp_dosage_genotyp_val_and_normalized_ = (2*tmp_dosage_genotyp_val_and_ - 1 - A_p_arm2_perm_and_alpha_(ndud)) * sqrt(A_p_arm2_perm_and_delta_(ndud));
tmp_dosage_genotyp_val_xor_normalized_ = (2*tmp_dosage_genotyp_val_xor_ - 1 - A_p_arm2_perm_xor_alpha_(ndud)) * sqrt(A_p_arm2_perm_xor_delta_(ndud));
tmp_dosage_genotyp_val_nor_normalized_ = (2*tmp_dosage_genotyp_val_nor_ - 1 - A_p_arm2_perm_nor_alpha_(ndud)) * sqrt(A_p_arm2_perm_nor_delta_(ndud));
end;%if flag_normalize_A_p & flag_normalize_arm2
end;%if set to zero. ;
%%%%%%%%;
tmp_dosage_genotyp_row_and_normalized_(ndud) = tmp_dosage_genotyp_val_and_normalized_;
tmp_dosage_genotyp_row_xor_normalized_(ndud) = tmp_dosage_genotyp_val_xor_normalized_;
tmp_dosage_genotyp_row_nor_normalized_(ndud) = tmp_dosage_genotyp_val_nor_normalized_;
%%%%%%%%;
if flag_V_trnx;
tmp_dosage_imputed_row_AnV_ = tmp_dosage_imputed_row_AnV_ ...
   + tmp_dosage_imputed_val_and_normalized_*V_rx_trnx_perm_and_(ndud,:) ...
   + tmp_dosage_imputed_val_xor_normalized_*V_rx_trnx_perm_xor_(ndud,:) ...
   + tmp_dosage_imputed_val_nor_normalized_*V_rx_trnx_perm_nor_(ndud,:) ...
  ;
tmp_dosage_rounded_row_AnV_ = tmp_dosage_rounded_row_AnV_ ...
   + tmp_dosage_rounded_val_and_normalized_*V_rx_trnx_perm_and_(ndud,:) ...
   + tmp_dosage_rounded_val_xor_normalized_*V_rx_trnx_perm_xor_(ndud,:) ...
   + tmp_dosage_rounded_val_nor_normalized_*V_rx_trnx_perm_nor_(ndud,:) ...
  ;
tmp_dosage_genotyp_row_AnV_ = tmp_dosage_genotyp_row_AnV_ ...
   + tmp_dosage_genotyp_val_and_normalized_*V_rx_trnx_perm_and_(ndud,:) ...
   + tmp_dosage_genotyp_val_xor_normalized_*V_rx_trnx_perm_xor_(ndud,:) ...
   + tmp_dosage_genotyp_val_nor_normalized_*V_rx_trnx_perm_nor_(ndud,:) ...
  ;
n_z_a = n_z_a + (V_rx_trnx_perm_and_(ndud,:)~=0) + (V_rx_trnx_perm_xor_(ndud,:)~=0) + (V_rx_trnx_perm_nor_(ndud,:)~=0);
end;%if flag_V_trnx;
if flag_V_tstx;
tmp_dosage_imputed_row_AnV_ = tmp_dosage_imputed_row_AnV_ ...
   + tmp_dosage_imputed_val_and_normalized_*V_rx_tstx_perm_and_(ndud,:) ...
   + tmp_dosage_imputed_val_xor_normalized_*V_rx_tstx_perm_xor_(ndud,:) ...
   + tmp_dosage_imputed_val_nor_normalized_*V_rx_tstx_perm_nor_(ndud,:) ...
  ;
tmp_dosage_rounded_row_AnV_ = tmp_dosage_rounded_row_AnV_ ...
   + tmp_dosage_rounded_val_and_normalized_*V_rx_tstx_perm_and_(ndud,:) ...
   + tmp_dosage_rounded_val_xor_normalized_*V_rx_tstx_perm_xor_(ndud,:) ...
   + tmp_dosage_rounded_val_nor_normalized_*V_rx_tstx_perm_nor_(ndud,:) ...
  ;
tmp_dosage_genotyp_row_AnV_ = tmp_dosage_genotyp_row_AnV_ ...
   + tmp_dosage_genotyp_val_and_normalized_*V_rx_tstx_perm_and_(ndud,:) ...
   + tmp_dosage_genotyp_val_xor_normalized_*V_rx_tstx_perm_xor_(ndud,:) ...
   + tmp_dosage_genotyp_val_nor_normalized_*V_rx_tstx_perm_nor_(ndud,:) ...
  ;
n_z_a = n_z_a + (V_rx_tstx_perm_and_(ndud,:)~=0) + (V_rx_tstx_perm_xor_(ndud,:)~=0) + (V_rx_tstx_perm_nor_(ndud,:)~=0);
end;%if flag_V_tstx;
end;%if (length(neu1)==1);
end;%for ndosage_snp=1:n_dosage_snp;

%%%%%%%%;
disp(sprintf(' %% flag_bo1_all_sum = %d; flag_bo2_and_sum = %d; flag_bo2_xor_sum = %d; flag_bo2_nor_sum = %d',flag_bo1_all_sum,flag_bo2_and_sum,flag_bo2_xor_sum,flag_bo2_nor_sum));
tmp_h_mismatch_d1_imputed_ = hist(mismatch_d1_imputed_,[0,1,2,3,8,9]);
tmp_h_mismatch_d1_genotyp_ = hist(mismatch_d1_genotyp_,[0,1,2,3,8,9]);
disp(sprintf(' %% tmp_h_d1_mismatched_imputed_: %0.10d %0.10d %0.10d %0.10d %0.10d %0.10d',tmp_h_mismatch_d1_imputed_));
disp(sprintf(' %% tmp_h_d1_mismatched_genotyp_: %0.10d %0.10d %0.10d %0.10d %0.10d %0.10d',tmp_h_mismatch_d1_genotyp_));
tmp_h_mismatch_d2_imputed_ = hist(mismatch_d2_imputed_,[0,1,2,3,8,9]);
tmp_h_mismatch_d2_genotyp_ = hist(mismatch_d2_genotyp_,[0,1,2,3,8,9]);
disp(sprintf(' %% tmp_h_d2_mismatched_imputed_: %0.10d %0.10d %0.10d %0.10d %0.10d %0.10d',tmp_h_mismatch_d2_imputed_));
disp(sprintf(' %% tmp_h_d2_mismatched_genotyp_: %0.10d %0.10d %0.10d %0.10d %0.10d %0.10d',tmp_h_mismatch_d2_genotyp_));
tmp_mismatch_d1_imputed_swap_ = ( (mismatch_d1_imputed_==error_code_dominance) | (mismatch_d1_imputed_==error_code_both) );
tmp_mismatch_d1_genotyp_swap_ = ( (mismatch_d1_genotyp_==error_code_dominance) | (mismatch_d1_genotyp_==error_code_both) );
tmp_mismatch_d2_imputed_swap_ = ( (mismatch_d2_imputed_==error_code_dominance) | (mismatch_d2_imputed_==error_code_both) );
tmp_mismatch_d2_genotyp_swap_ = ( (mismatch_d2_genotyp_==error_code_dominance) | (mismatch_d2_genotyp_==error_code_both) );
tmp_mismatch_d1_imputed_keep_ = ( (mismatch_d1_imputed_==error_code_null) | (mismatch_d1_imputed_==error_code_synonym) );
tmp_mismatch_d1_genotyp_keep_ = ( (mismatch_d1_genotyp_==error_code_null) | (mismatch_d1_genotyp_==error_code_synonym) );
tmp_mismatch_d2_imputed_keep_ = ( (mismatch_d2_imputed_==error_code_null) | (mismatch_d2_imputed_==error_code_synonym) );
tmp_mismatch_d2_genotyp_keep_ = ( (mismatch_d2_genotyp_==error_code_null) | (mismatch_d2_genotyp_==error_code_synonym) );
tmp_mismatch_dx_imputed_keep_ = tmp_mismatch_d1_imputed_keep_ & tmp_mismatch_d2_imputed_keep_ ;
tmp_mismatch_dx_genotyp_keep_ = tmp_mismatch_d1_genotyp_keep_ & tmp_mismatch_d2_genotyp_keep_ ;
tmp_mismatch_dx_imputed_swap_ = tmp_mismatch_d1_imputed_swap_ & tmp_mismatch_d2_imputed_swap_ ;
tmp_mismatch_dx_genotyp_swap_ = tmp_mismatch_d1_genotyp_swap_ & tmp_mismatch_d2_genotyp_swap_ ;
tmp_mismatch_dx_imputed_miss_ = tmp_mismatch_d1_imputed_swap_ & tmp_mismatch_d2_imputed_keep_ ;
tmp_mismatch_dx_genotyp_miss_ = tmp_mismatch_d1_genotyp_swap_ & tmp_mismatch_d2_genotyp_keep_ ;
assert(sum( tmp_mismatch_dx_imputed_keep_ & tmp_mismatch_dx_imputed_swap_ )==0);
assert(sum( tmp_mismatch_dx_imputed_keep_ & tmp_mismatch_dx_imputed_miss_ )==0);
assert(sum( tmp_mismatch_dx_imputed_swap_ & tmp_mismatch_dx_imputed_miss_ )==0);
assert(sum( tmp_mismatch_dx_genotyp_keep_ & tmp_mismatch_dx_genotyp_swap_ )==0);
assert(sum( tmp_mismatch_dx_genotyp_keep_ & tmp_mismatch_dx_genotyp_miss_ )==0);
assert(sum( tmp_mismatch_dx_genotyp_swap_ & tmp_mismatch_dx_genotyp_miss_ )==0);
disp(sprintf(' %% imputed: keep %d swap %d miss %d; genotyp: keep %d swap %d miss %d'...
,sum(tmp_mismatch_dx_imputed_keep_),sum(tmp_mismatch_dx_imputed_swap_),sum(tmp_mismatch_dx_imputed_miss_)...
,sum(tmp_mismatch_dx_genotyp_keep_),sum(tmp_mismatch_dx_genotyp_swap_),sum(tmp_mismatch_dx_genotyp_miss_)...
	     ));
%%%%%%%%;
