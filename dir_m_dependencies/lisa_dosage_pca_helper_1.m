function error_code = lisa_dosage_pca_helper_1(a1,a2,b1,b2,f_a_and,f_a_xor,f_a_nor,f_b_and,f_b_xor,f_b_nor,ent_cutoff);
% returns error type associated with alleles a1,a2 (from one list) and b1,b2 (from another list). ;
% If frequencies f_a and f_b are given, we use them. ;
%%%%%%%%%%%%%%%%;
% We define different kinds of allelic mismatch errors: ;
%%%%%%%%%%%%%%%%;
error_code_null = 0;      % <-- No error; alleles match and minor-allele-frequency is commensurate. ;
error_code_synonym = 1;   % <-- Synonym; allele names must be interchanged with synonyms (i.e,. A<-->T and C<-->G). ;
error_code_dominance = 2; % <-- Dominance; the alleles must be switched (i.e., XY --> YX) so that the first allele listed corresponds to the minor (i.e., less frequent) allele. ;
error_code_both = 3;      % <-- Both; the allele names must be interchanged *and* the order must be switched. ;
error_code_ambiguous = 8; % <-- Ambiguity; I cannot determine the correct alignment. ;
error_code_missing = 9;   % <-- Missing; this snp is not in the overlap between the two sets. ;

error_code = error_code_ambiguous;

if nargin<5; flag_f=0; end;
if nargin>5; flag_f=1; end;
if nargin<11; ent_cutoff = 0.03; end;

if flag_f;
tmp_a = f_a_and + f_a_xor + f_a_nor;
f_a_and = f_a_and/tmp_a;
f_a_xor = f_a_xor/tmp_a;
f_a_nor = f_a_nor/tmp_a;
p_a_opt = f_a_and + 0.5*f_a_xor; q_a_opt = f_a_nor + 0.5*f_a_xor;
I_a_opt = f_a_and.*log(f_a_and./(p_a_opt.^2)) + ...
             f_a_xor.*log(f_a_xor./(2*p_a_opt.*q_a_opt)) + ...
             f_a_nor.*log(f_a_nor./(q_a_opt.^2)) ;
tmp_b = f_b_and + f_b_xor + f_b_nor;
f_b_and = f_b_and/tmp_b;
f_b_xor = f_b_xor/tmp_b;
f_b_nor = f_b_nor/tmp_b;
p_b_opt = f_b_and + 0.5*f_b_xor; q_b_opt = f_b_nor + 0.5*f_b_xor;
I_b_opt = f_b_and.*log(f_b_and./(p_b_opt.^2)) + ...
             f_b_xor.*log(f_b_xor./(2*p_b_opt.*q_b_opt)) + ...
             f_b_nor.*log(f_b_nor./(q_b_opt.^2)) ;
flag_p_hweq = (I_a_opt < ent_cutoff) & (I_b_opt < ent_cutoff) ;
flag_p_match = ((p_a_opt-0.5)*(p_b_opt-0.5)>0);
flag_q_match = ((p_a_opt-0.5)*(q_b_opt-0.5)>0);
flag_p_close = (abs(p_a_opt-p_b_opt)<0.05);
end;%if flag_f;

allele_original__ = zeros(84,1); 
allele_original__(65) = 1; %<-- A ;
allele_original__(67) = 2; %<-- C ;
allele_original__(71) = 3; %<-- G ;
allele_original__(84) = 4; %<-- T ;
allele_complement__ = zeros(84,1); 
allele_complement__(65) = 4; %<-- A<-T ;
allele_complement__(67) = 3; %<-- C<-G ;
allele_complement__(71) = 2; %<-- G<-C ;
allele_complement__(84) = 1; %<-- T<-A ;
%%%%%%%%;
flag_a_synonym = ( (allele_original_(a1)==allele_original_(a2)) | (allele_original_(a1)==allele_complement_(a2)) );
flag_b_synonym = ( (allele_original_(b1)==allele_original_(b2)) | (allele_original_(b1)==allele_complement_(b2)) );
flag_a_distinct = ( (allele_original_(a1)~=allele_original_(a2)) & (allele_original_(a1)~=allele_complement_(a2)) );
flag_b_distinct = ( (allele_original_(b1)~=allele_original_(b2)) & (allele_original_(b1)~=allele_complement_(b2)) );
flag_ab_synonym = (flag_a_synonym | flag_b_synonym);
flag_ab_distinct = (flag_a_distinct & flag_b_distinct);
%%%%%%%%;
flag_11_match = (allele_original_(a1)==allele_original_(b1));
flag_11_synonym = (allele_original_(a1)==allele_complement_(b1));
flag_22_match = (allele_original_(a2)==allele_original_(b2));
flag_22_synonym = (allele_original_(a2)==allele_complement_(b2));
flag_1122_match = (flag_11_match & flag_22_match);
flag_1122_synonym = (flag_11_match & flag_22_synonym) | (flag_11_synonym & flag_22_match) | (flag_11_synonym & flag_22_synonym) ;
%%%%%%%%;
flag_12_match = (allele_original_(a1)==allele_original_(b2));
flag_12_synonym = (allele_original_(a1)==allele_complement_(b2));
flag_21_match = (allele_original_(a2)==allele_original_(b1));
flag_21_synonym = (allele_original_(a2)==allele_complement_(b1));
flag_1221_match = (flag_12_match & flag_21_match);
flag_1221_synonym = (flag_12_match & flag_21_synonym) | (flag_12_synonym & flag_21_match) | (flag_12_synonym & flag_21_synonym) ;
%%%%%%%%;

if flag_ab_synonym;
if flag_f==0; error_code = error_code_ambiguous; end;
if flag_f==1; if flag_p_hweq;
if (  flag_p_close ) ; error_code = error_code_ambiguous; end;
if (  flag_p_match & ~flag_p_close ) ; error_code = error_code_null; end;
if ( ~flag_p_match & ~flag_p_close ) ; error_code = error_code_dominance; end;
end;end;%if flag_f==1; if flag_p_hweq;
end;%if flag_ab_synonym;

if flag_ab_distinct;
%%%%%%%%;
if (  flag_1122_match );
if flag_f==0; error_code = error_code_null; end;
if flag_f==1; if flag_p_hweq;
if (  flag_p_close ); error_code = error_code_null; end; 
if (  flag_p_match & ~flag_p_close ) ; error_code = error_code_null ; end;
if ( ~flag_p_match & ~flag_p_close ) ; error_code = error_code_ambiguous ; end;
end;end;%if flag_f==1; if flag_p_hweq;
end;%if (  flag_1122_match );
%%%%%%%%;
if (  flag_1221_match );
if flag_f==0; error_code = error_code_dominance; end;
if flag_f==1; if flag_p_hweq;
if (  flag_p_close ); error_code = error_code_dominance; end; 
if (  flag_q_match & ~flag_p_close ) ; error_code = error_code_dominance ; end;
if ( ~flag_q_match & ~flag_p_close ) ; error_code = error_code_ambiguous ; end;
end;end;%if flag_f==1; if flag_p_hweq;
end;%if (  flag_1221_match );
%%%%%%%%;
if (  flag_1122_synonym );
if flag_f==0; error_code = error_code_synonym; end;
if flag_f==1; if flag_p_hweq;
if (  flag_p_close ); error_code = error_code_synonym; end; 
if (  flag_p_match & ~flag_p_close ) ; error_code = error_code_synonym ; end;
if ( ~flag_p_match & ~flag_p_close ) ; error_code = error_code_ambiguous ; end;
end;end;%if flag_f==1; if flag_p_hweq;
end;%if (  flag_1122_synonym );
%%%%%%%%;
if (  flag_1221_synonym );
if flag_f==0; error_code = error_code_both; end;
if flag_f==1; if flag_p_hweq;
if (  flag_p_close ); error_code = error_code_both; end; 
if (  flag_q_match & ~flag_p_close ) ; error_code = error_code_both ; end;
if ( ~flag_q_match & ~flag_p_close ) ; error_code = error_code_ambiguous ; end;
end;end;%if flag_f==1; if flag_p_hweq;
end;%if (  flag_1221_synonym );
%%%%%%%%;
end;%if flag_ab_distinct;


