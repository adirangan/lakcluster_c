void GLOBAL_ops_reset_all(){ int nt; for (nt=0;nt<MAX_THREADS;nt++){ GLOBAL_ops_f_[nt]=0; GLOBAL_ops_b_[nt]=0;}}

void GLOBAL_ops_reset_one(int nt){ GLOBAL_ops_f_[nt]=0; GLOBAL_ops_b_[nt]=0;}

void GLOBAL_ops_count_one(int nt,unsigned long long int f,unsigned long long int b){ GLOBAL_ops_f_[nt] += f; GLOBAL_ops_b_[nt] += b;}

void GLOBAL_ops_addup_all(){ int nt; for (nt=0;nt<MAX_THREADS;nt++){ GLOBAL_ops_f_sum += GLOBAL_ops_f_[nt]; GLOBAL_ops_b_sum += GLOBAL_ops_b_[nt];}}

void GLOBAL_ops_addup_one(int nt){ GLOBAL_ops_f_sum += GLOBAL_ops_f_[nt]; GLOBAL_ops_b_sum += GLOBAL_ops_b_[nt];}

void GLOBAL_ops_printf_all(int verbose,char *prefix)
{ 
  int nt=0;
  if (verbose>1){
    for (nt=0;nt<MAX_THREADS;nt++){
      if (GLOBAL_ops_f_[nt]>0 || GLOBAL_ops_b_[nt]>0){ printf("%snt %d: \tflops 10^%0.2f = %14ld \tblops 10^%0.2f = %14ld\n",prefix,nt,log10(GLOBAL_ops_f_[nt]),GLOBAL_ops_f_[nt],log10(GLOBAL_ops_b_[nt]),GLOBAL_ops_b_[nt]);}
      /* for (nt=0;nt<MAX_THREADS;nt++){ } */}
    /* if (verbose>1){ } */}
  if (verbose>0){
    printf("%sall : \tflops 10^%0.2f = %14ld \tblops 10^%0.2f = %14ld\n",prefix,log10(GLOBAL_ops_f_sum),GLOBAL_ops_f_sum,log10(GLOBAL_ops_b_sum),GLOBAL_ops_b_sum);
    /* if (verbose>0){ } */}
}

void GLOBAL_tic(int nx){ nx = maximum(0,minimum(GLOBAL_NTICKS-1,nx)); GLOBAL_t_start[nx] = clock(); gettimeofday(&GLOBAL_d_start[nx],NULL);}

void GLOBAL_toc(int nx,int verbose,char *prefix)
{ 
  double r=0;
  nx = maximum(0,minimum(GLOBAL_NTICKS-1,nx)); 
  GLOBAL_t_final[nx] = clock(); gettimeofday(&GLOBAL_d_final[nx],NULL); 
  GLOBAL_l_ssec[nx] =  GLOBAL_d_final[nx].tv_sec -  GLOBAL_d_start[nx].tv_sec; 
  GLOBAL_l_usec[nx] = GLOBAL_d_final[nx].tv_usec - GLOBAL_d_start[nx].tv_usec; 
  GLOBAL_l_msec[nx] = ((GLOBAL_l_ssec[nx]*1000) + GLOBAL_l_usec[nx]/1000.0) + 0.5; 
  GLOBAL_elct[nx] = (double)(1000*(GLOBAL_t_final[nx]-GLOBAL_t_start[nx])/CLOCKS_PER_SEC)/(double)1000; 
  GLOBAL_elrt[nx] = (double)GLOBAL_l_msec[nx]/(double)1000; 
  r = GLOBAL_elct[nx]/GLOBAL_elrt[nx];
  if (verbose){ 
    if (finite(r)){ printf("%sct/rt %0.3f/%0.3f = %.1f,",prefix,GLOBAL_elct[nx],GLOBAL_elrt[nx],r); wkspace_printf();}
    else{ printf("%sct/rt %0.3f/%0.3f = 0,",prefix,GLOBAL_elct[nx],GLOBAL_elrt[nx]); wkspace_printf();}
    /* if (verbose){ } */}
}

void GLOBAL_ops_toc(int nt,int nx,int verbose,char *prefix)
{ 
  double r=0;
  unsigned long long int f=0;
  unsigned long long int b=0;
  nx = maximum(0,minimum(GLOBAL_NTICKS-1,nx)); 
  GLOBAL_t_final[nx] = clock(); gettimeofday(&GLOBAL_d_final[nx],NULL); 
  GLOBAL_l_ssec[nx] =  GLOBAL_d_final[nx].tv_sec -  GLOBAL_d_start[nx].tv_sec; 
  GLOBAL_l_usec[nx] = GLOBAL_d_final[nx].tv_usec - GLOBAL_d_start[nx].tv_usec; 
  GLOBAL_l_msec[nx] = ((GLOBAL_l_ssec[nx]*1000) + GLOBAL_l_usec[nx]/1000.0) + 0.5; 
  GLOBAL_elct[nx] = (double)(1000*(GLOBAL_t_final[nx]-GLOBAL_t_start[nx])/CLOCKS_PER_SEC)/(double)1000; 
  GLOBAL_elrt[nx] = (double)GLOBAL_l_msec[nx]/(double)1000; 
  r = GLOBAL_elct[nx]/GLOBAL_elrt[nx];
  if (nt<0){ f = GLOBAL_ops_f_sum; b = GLOBAL_ops_b_sum;} else{ f = GLOBAL_ops_f_[nt]; b = GLOBAL_ops_b_[nt];}
  if (verbose){ 
    if (finite(r)){ printf("%sct/rt %0.3f/%0.3f = %.1f, %0.2f Gf %0.2f Gb,",prefix,GLOBAL_elct[nx],GLOBAL_elrt[nx],r,f/GLOBAL_elrt[nx]/1e9,b/GLOBAL_elrt[nx]/1e9); wkspace_printf();}
    else{ printf("%sct/rt %0.3f/%0.3f = 0, ?? ",prefix,GLOBAL_elct[nx],GLOBAL_elrt[nx]); wkspace_printf();}
    /* if (verbose){ } */}
}

void GLOBAL_pthread_tic(){ GLOBAL_nf_cur = periodize(GLOBAL_nf_cur+1,0,GLOBAL_thread_count); if (GLOBAL_verbose>1){ printf(" nf_cur: %d[%d/%d];\n",GLOBAL_nf_cur,GLOBAL_nf_opn,GLOBAL_thread_count);}}

void GLOBAL_pthread_toc(){ GLOBAL_nf_opn += GLOBAL_nf_cur>0; if (GLOBAL_nf_cur==0){ for (GLOBAL_nf_cur=1;GLOBAL_nf_cur<GLOBAL_thread_count;GLOBAL_nf_cur++){ pthread_join(GLOBAL_threads[GLOBAL_nf_cur],NULL);} GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;}}

void GLOBAL_pthread_tuc(){ if (GLOBAL_nf_opn && GLOBAL_nf_cur==0){ for (GLOBAL_nf_cur=1;GLOBAL_nf_cur<GLOBAL_thread_count;GLOBAL_nf_cur++){ pthread_join(GLOBAL_threads[GLOBAL_nf_cur],NULL);} GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;} else if (GLOBAL_nf_opn && GLOBAL_nf_cur>0){ for (GLOBAL_nf_ind=1;GLOBAL_nf_ind<=GLOBAL_nf_cur;GLOBAL_nf_ind++){ pthread_join(GLOBAL_threads[GLOBAL_nf_ind],NULL);} GLOBAL_nf_ind=0; GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;}}
