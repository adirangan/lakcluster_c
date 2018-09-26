%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Compile with:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make -f lakcluster_ver18.make lakcluster_ver18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Run with:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
./lakcluster_ver18 < some_input_file.in ;
(example input_files below):

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Test for errors with:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
./lakcluster_ver18 < dir_in/An_ajdk_v_error.in
./lakcluster_ver18 < dir_in/AnAt_vv_error.in
./lakcluster_ver18 < dir_in/An_v_error.in
./lakcluster_ver18 < dir_in/An_ZtSn_ww_error.in
./lakcluster_ver18 < dir_in/An_ZtSWn_ww_error.in
./lakcluster_ver18 < dir_in/An_ZtSWn_Yt_vv_error.in
./lakcluster_ver18 < dir_in/AnZt_S_WnYt_vv_error.in
./lakcluster_ver18 < dir_in/An_ZtSWn_Yt_ww_error.in
./lakcluster_ver18 < dir_in/AnZt_vv_error.in
./lakcluster_ver18 < dir_in/AnZt_xx_error.in
./lakcluster_ver18 < dir_in/AtTAn_vv_error.in
./lakcluster_ver18 < dir_in/AtTYn_vv_error.in
./lakcluster_ver18 < dir_in/AttYn____WtsZn_vv_error.in
./lakcluster_ver18 < dir_in/At_T_YnWt_S_Zn_vv_error.in
./lakcluster_ver18 < dir_in/AtTYn____WtSZn_vv_error.in
./lakcluster_ver18 < dir_in/At_T_YnWt_S_Zn_ww_error.in
./lakcluster_ver18 < dir_in/At_T_YnWt_ww_error.in
./lakcluster_ver18 < dir_in/bcc_flattenloop_error.in
./lakcluster_ver18 < dir_in/bcc_init_error.in
./lakcluster_ver18 < dir_in/bcc_lf_AnZtSWnYt_error.in
./lakcluster_ver18 < dir_in/bcc_lf_AtTYnWtSZn_error.in
./lakcluster_ver18 < dir_in/bcc_lrup_error.in
./lakcluster_ver18 < dir_in/bcc_lrup_flattenloop_error.in
./lakcluster_ver18 < dir_in/bcc_QX_error.in
./lakcluster_ver18 < dir_in/bcc_sumscores_error.in
./lakcluster_ver18 < dir_in/binary_M_setup_test.in
./lakcluster_ver18 < dir_in/D_AtTn_ZtSn_vv_error.in
./lakcluster_ver18 < dir_in/dcc_init_error.in
./lakcluster_ver18 < dir_in/dcc_lf_D_AtTn_ZtSn_error.in
./lakcluster_ver18 < dir_in/dcc_lf_TAnZtS_error.in
./lakcluster_ver18 < dir_in/dcc_sumscores_error.in

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Test timing with:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
./lakcluster_ver18 < dir_in/An_ajdk_v_speed.in
./lakcluster_ver18 < dir_in/AnAt_vv_speed.in
./lakcluster_ver18 < dir_in/An_v_speed.in
./lakcluster_ver18 < dir_in/An_ZtSn_ww_speed.in
./lakcluster_ver18 < dir_in/An_ZtSWn_Yt_vv_speed.in
./lakcluster_ver18 < dir_in/An_ZtSWn_Yt_ww_speed.in
./lakcluster_ver18 < dir_in/AnZt_vv_speed2.in
./lakcluster_ver18 < dir_in/AnZt_vv_speed.in
./lakcluster_ver18 < dir_in/AtTYn_vv_speed.in
./lakcluster_ver18 < dir_in/At_T_YnWt_S_Zn_ww_speed.in
./lakcluster_ver18 < dir_in/At_T_YnWt_ww_speed.in
./lakcluster_ver18 < dir_in/bcc_flattenloop_speed.in
./lakcluster_ver18 < dir_in/bcc_lrup_sumscores_speed.in
./lakcluster_ver18 < dir_in/bcc_sumscores_speed2.in
./lakcluster_ver18 < dir_in/bcc_sumscores_speed.in
./lakcluster_ver18 < dir_in/D_AtTn_ZtSn_vv_speed.in
./lakcluster_ver18 < dir_in/dcc_time_sumscores_speed.in
./lakcluster_ver18 < dir_in/ynwt_vv_speed.in
./lakcluster_ver18 < dir_in/ztsvn_vv_speed.in

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Matlab drivers included within ./dir_m/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Each of these drivers can be run (in Matlab).
Depending on your Matlab path, it may be necessary to
change directory to the subdirectory "./dir_m".
Each driver takes the following 8 input arguments:

1. dir_trunk (string): The name of the directory from which the
   command is called (typically "pwd").

2. N (integer): Parameter governing the size of the data-array (some
   drivers will increase or decrease the number of columns or rows as
   necessary).

3. X_factor (double): Parameter governing the size of the 'important'
   bicluster embedded in the data-array (usually the one we are trying
   to find). The size "m" of the important bicluster is typically N
   raised to the power of X_factor. Other spurious biclusters are
   often larger.

4. X_esm (double): Parameter governing the noise of the 'important'
   bicluster embedded in the data-array. This parameter is equal to
   the spectral noise 'epsilon' multiplied by sqrt(n).

5. gamma (double): The fraction of rows/columns thrown out each
   step. If gamma=0.5, then half the remaining rows/columns are
   ejected each step. Smaller values of gamma are better. We usually
   pick gamma - 0.05 or smaller.

6. B_MLT (integer): The number of bits retained within internal
   matrix-matrix multiplication. B_MLT=24 or 32 is usually sufficient
   for an accurate calculation, however in certain situations B_MLT
   can be as low as 8 with good results.

7. rng_num (integer): The random-number-generator seed set before
   generating random matrices. With this input argument we ensure that
   subsequent calls to the same driver (with the same inputs)
   construct the same example.

8. pt_num (integer): This "part-number" specifies the suffix that will
   be added to the output file. This is only used when passing
   multiple calls to the same driver within a parallel computing
   environment. If not used, either omit or set to "-1".

Example drivers include:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lakcluster_test_AAAA_ver5.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests the loop-counting algorithm using a simple case-only example
with no controls or covariates. For example, you may try:
>> N=1024;X_factor=0.55;X_esm=0.1;gamma=0.05;B_MLT=32;rng_num=1;
>> lakcluster_test_AAAA_ver5(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lakcluster_test_AAYY_ver2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests example with genetic-controls (no covariates). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lakcluster_test_ADAADA_ver3.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests case-only example with varying sparsity (no controls or
covariates).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lakcluster_test_AZZA_ver2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests example with patient-controls as well as cases (no covariates).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lakcluster_test_uAAAA_ver1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests case-only example with categorical covariates.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lakcluster_test_AATAA_ver4.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests case-only example with continuous-covariates. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lakcluster_test_AZWY_ver2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests example with patient- and genetic-controls (no covariates).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lakcluster_test_uAZZA_ver4.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests example with patient-controls and categorical covariates.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lakcluster_test_uADZSZDA_Ireq1_ver7.m
lakcluster_test_uADZSZDA_Ireq2_ver7.m
lakcluster_test_uADZSZDA_Ireq3_ver7.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
These test examples with patient-controls, 3 categorical-covariates
and 2 continuous-covariates, as well as varying sparsity
coefficients. The parameter "Ireq" refers to the number of
categorical-covariates required for an embedded bicluster to be
considered by the algorithm.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dexcluster_test_AAAA_ver5.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests half-loop algorithm with case-only example (no controls or
covariates).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dexcluster_test_uAZSZA_ver5.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tests half-loop algorithm with patient-controls and categorical
covariates.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dexcluster_test_uADZSZDA_Ireq1_ver7.m
dexcluster_test_uADZSZDA_Ireq2_ver7.m
dexcluster_test_uADZSZDA_Ireq3_ver7.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
These are analogous to lakcluster_test_uADZSZDA_Ireq?_ver7.m, except
that they test the half-loop algorithm rather than the full
loop-counting algorithm.



