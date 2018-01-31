void updateglobals(char *vname)
{
  /* given a variable name, scan the appropriate format into the appropriate global variable */
  int verbose=GLOBAL_verbose;
  char comma_vs_semicolon[1],tmpchar[FNAMESIZE]; int length=0,nv=0;
  if (0){ /* do nothing */ }
  else if (strcmp(vname,"GLOBAL_verbose")==0){ scanf("%d",&GLOBAL_verbose); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_verbose);}}
  else if (strcmp(vname,"GLOBAL_thread_count")==0){ 
    scanf("%d",&GLOBAL_thread_count); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_thread_count);}
    length = GLOBAL_thread_count; 
    /* else if (strcmp(vname,"GLOBAL_thread_count")==0){ } */}
  else if (strcmp(vname,"GLOBAL_omp_type")==0){ scanf("%d",&GLOBAL_omp_type); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_omp_type);}}
  else if (strcmp(vname,"GLOBAL_LBITS")==0){ scanf("%d",&GLOBAL_LBITS); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_LBITS);}}
  else if (strcmp(vname,"GLOBAL_TEST_TYPE")==0){ scanf("%[^,;]",GLOBAL_TEST_TYPE);  if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_TEST_TYPE);}}
  else if (strcmp(vname,"GLOBAL_TEST_TYP2")==0){ scanf("%[^,;]",GLOBAL_TEST_TYP2);  if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_TEST_TYP2);}}
  else if (strcmp(vname,"GLOBAL_TEST_mrand")==0){ scanf("%lf",&GLOBAL_TEST_mrand); if (verbose>0){ printf("%s read to be %0.2f\n",vname,GLOBAL_TEST_mrand);}}
  else if (strcmp(vname,"GLOBAL_TEST_sparse")==0){ scanf("%d",&GLOBAL_TEST_sparse); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_TEST_sparse);}}
  else if (strcmp(vname,"GLOBAL_TEST_A_n_rows")==0){ scanf("%d",&GLOBAL_TEST_A_n_rows); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_TEST_A_n_rows);}}
  else if (strcmp(vname,"GLOBAL_TEST_A_n_cols")==0){ scanf("%d",&GLOBAL_TEST_A_n_cols); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_TEST_A_n_cols);}}
  else if (strcmp(vname,"GLOBAL_TEST_Z_n_rows")==0){ scanf("%d",&GLOBAL_TEST_Z_n_rows); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_TEST_Z_n_rows);}}
  else if (strcmp(vname,"GLOBAL_TEST_Y_n_cols")==0){ scanf("%d",&GLOBAL_TEST_Y_n_cols); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_TEST_Y_n_cols);}}
  else if (strcmp(vname,"GLOBAL_TEST_T_n_cols")==0){ scanf("%d",&GLOBAL_TEST_T_n_cols); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_TEST_T_n_cols);}}
  else if (strcmp(vname,"GLOBAL_TEST_niter")==0){ scanf("%d",&GLOBAL_TEST_niter); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_TEST_niter);}}
  else if (strcmp(vname,"GLOBAL_QR_strategy")==0){ scanf("%[^,;]",GLOBAL_QR_strategy);  if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_QR_strategy);}}
  else if (strcmp(vname,"GLOBAL_QC_strategy")==0){ scanf("%[^,;]",GLOBAL_QC_strategy);  if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_QC_strategy);}}
  else if (strcmp(vname,"GLOBAL_D_MLT")==0){ scanf("%lld",&GLOBAL_D_MLT); if (verbose>0){ printf("%s read to be %lld\n",vname,GLOBAL_D_MLT);} GLOBAL_B_MLT = (int)llround(log2((double)GLOBAL_D_MLT));}
  else if (strcmp(vname,"GLOBAL_B_MLT")==0){ scanf("%d",&GLOBAL_B_MLT); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_B_MLT);} GLOBAL_D_MLT = (long long int)pow(2.0,(double)GLOBAL_B_MLT);}
  else if (strcmp(vname,"GLOBAL_gamma_type")==0){ scanf("%d",&GLOBAL_gamma_type); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_gamma_type);}}
  else if (strcmp(vname,"GLOBAL_gamma")==0){ scanf("%lf",&GLOBAL_gamma); if (verbose>0){ printf("%s read to be %f\n",vname,GLOBAL_gamma);}}
  else if (strcmp(vname,"GLOBAL_Ireq")==0){ scanf("%d",&GLOBAL_Ireq); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_Ireq);}}
  else if (strcmp(vname,"GLOBAL_NBINS")==0){ 
    scanf("%d",&GLOBAL_NBINS); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_NBINS);} 
    if (GLOBAL_NBINS>0){ 
      length=GLOBAL_NBINS; 
      GLOBAL_A_n_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_A_t_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_A_n_rows_ = (int *) wkspace_alloc(GLOBAL_NBINS*sizeof(int)); fill_uchar_zero((unsigned char *)GLOBAL_A_n_rows_,GLOBAL_NBINS*sizeof(int)); for (nv=0;nv<length;nv++){ GLOBAL_A_n_rows_[nv]=0;}
      GLOBAL_A_n_rind_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_Z_n_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_Z_t_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_Z_n_rows_ = (int *) wkspace_alloc(GLOBAL_NBINS*sizeof(int)); fill_uchar_zero((unsigned char *)GLOBAL_Z_n_rows_,GLOBAL_NBINS*sizeof(int)); for (nv=0;nv<length;nv++){ GLOBAL_Z_n_rows_[nv]=0;}
      GLOBAL_Z_n_rind_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_Y_n_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_Y_t_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_W_n_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_W_t_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_T_n_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_T_t_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_S_n_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      GLOBAL_S_t_name_ = (char **) wkspace_alloc(GLOBAL_NBINS*sizeof(char *)); 
      for (nv=0;nv<length;nv++){ 
	GLOBAL_A_n_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_A_n_name_[nv][0]='\0';
	GLOBAL_A_t_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_A_t_name_[nv][0]='\0';
	GLOBAL_A_n_rind_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_A_n_rind_[nv][0]='\0';
	GLOBAL_Z_n_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_Z_n_name_[nv][0]='\0';
	GLOBAL_Z_t_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_Z_t_name_[nv][0]='\0';
	GLOBAL_Z_n_rind_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_Z_n_rind_[nv][0]='\0';
	GLOBAL_Y_n_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_Y_n_name_[nv][0]='\0';
	GLOBAL_Y_t_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_Y_t_name_[nv][0]='\0';
	GLOBAL_W_n_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_W_n_name_[nv][0]='\0';
	GLOBAL_W_t_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_W_t_name_[nv][0]='\0';
	GLOBAL_T_n_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_T_n_name_[nv][0]='\0';
	GLOBAL_T_t_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_T_t_name_[nv][0]='\0';
	GLOBAL_S_n_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_S_n_name_[nv][0]='\0';
	GLOBAL_S_t_name_[nv] = (char *) wkspace_alloc(FNAMESIZE*sizeof(char)); GLOBAL_S_t_name_[nv][0]='\0';
	/* for (nv=0;nv<length;nv++){ } */}
      /* if (GLOBAL_NBINS>0){ } */}
    /* else if (strcmp(vname,"GLOBAL_NBINS")==0){ } */}
  else if (strcmp(vname,"GLOBAL_A_n_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_A_n_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_A_n_name_[nv]); sprintf(GLOBAL_A_n_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_A_n_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_A_n_name")==0){ scanf("%[^,;]",GLOBAL_A_n_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_A_n_name);}} */
  else if (strcmp(vname,"GLOBAL_A_t_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_A_t_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_A_t_name_[nv]); sprintf(GLOBAL_A_t_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_A_t_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_A_t_name")==0){ scanf("%[^,;]",GLOBAL_A_t_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_A_t_name);}} */
  else if (strcmp(vname,"GLOBAL_A_n_rows_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_A_n_rows_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_A_n_rows_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_A_n_rows")==0){ scanf("%d",&GLOBAL_A_n_rows); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_A_n_rows);}} */
  else if (strcmp(vname,"GLOBAL_A_n_cols")==0){ scanf("%d",&GLOBAL_A_n_cols); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_A_n_cols);}}
  else if (strcmp(vname,"GLOBAL_A_n_rind_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_A_n_rind_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_A_n_rind_[nv]); sprintf(GLOBAL_A_n_rind_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_A_n_rind_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_A_n_rind")==0){ scanf("%[^,;]",GLOBAL_A_n_rind); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_A_n_rind);}} */
  else if (strcmp(vname,"GLOBAL_A_n_cind")==0){ scanf("%[^,;]",GLOBAL_A_n_cind); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_A_n_cind); sprintf(GLOBAL_A_n_cind,"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_A_n_cind);}}
  else if (strcmp(vname,"GLOBAL_Z_n_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_Z_n_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_Z_n_name_[nv]); sprintf(GLOBAL_Z_n_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_Z_n_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_Z_n_name")==0){ scanf("%[^,;]",GLOBAL_Z_n_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_Z_n_name);}} */
  else if (strcmp(vname,"GLOBAL_Z_t_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_Z_t_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_Z_t_name_[nv]); sprintf(GLOBAL_Z_t_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_Z_t_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_Z_t_name")==0){ scanf("%[^,;]",GLOBAL_Z_t_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_Z_t_name);}} */
  else if (strcmp(vname,"GLOBAL_Z_n_rows_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%d",&(GLOBAL_Z_n_rows_[nv])); if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %d\n",vname,nv,GLOBAL_Z_n_rows_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_Z_n_rows")==0){ scanf("%d",&GLOBAL_Z_n_rows); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_Z_n_rows);}} */
  else if (strcmp(vname,"GLOBAL_Z_n_rind_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_Z_n_rind_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_Z_n_rind_[nv]); sprintf(GLOBAL_Z_n_rind_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_Z_n_rind_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_Z_n_rind")==0){ scanf("%[^,;]",GLOBAL_Z_n_rind); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_Z_n_rind);}} */
  else if (strcmp(vname,"GLOBAL_Y_n_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_Y_n_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_Y_n_name_[nv]); sprintf(GLOBAL_Y_n_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_Y_n_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_Y_n_name")==0){ scanf("%[^,;]",GLOBAL_Y_n_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_Y_n_name);}} */
  else if (strcmp(vname,"GLOBAL_Y_t_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_Y_t_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_Y_t_name_[nv]); sprintf(GLOBAL_Y_t_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_Y_t_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_Y_t_name")==0){ scanf("%[^,;]",GLOBAL_Y_t_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_Y_t_name);}} */
  else if (strcmp(vname,"GLOBAL_Y_n_cols")==0){ scanf("%d",&GLOBAL_Y_n_cols); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_Y_n_cols);}}
  else if (strcmp(vname,"GLOBAL_Y_n_cind")==0){ scanf("%[^,;]",GLOBAL_Y_n_cind); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_Y_n_cind); sprintf(GLOBAL_Y_n_cind,"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_Y_n_cind);}}
  else if (strcmp(vname,"GLOBAL_W_n_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_W_n_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_W_n_name_[nv]); sprintf(GLOBAL_W_n_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_W_n_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_W_n_name")==0){ scanf("%[^,;]",GLOBAL_W_n_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_W_n_name);}} */
  else if (strcmp(vname,"GLOBAL_W_t_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_W_t_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_W_t_name_[nv]); sprintf(GLOBAL_W_t_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_W_t_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_W_t_name")==0){ scanf("%[^,;]",GLOBAL_W_t_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_W_t_name);}} */
  else if (strcmp(vname,"GLOBAL_T_n_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_T_n_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_T_n_name_[nv]); sprintf(GLOBAL_T_n_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_T_n_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_T_n_name")==0){ scanf("%[^,;]",GLOBAL_T_n_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_T_n_name);}} */
  else if (strcmp(vname,"GLOBAL_T_t_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_T_t_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_T_t_name_[nv]); sprintf(GLOBAL_T_t_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_T_t_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_T_t_name")==0){ scanf("%[^,;]",GLOBAL_T_t_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_T_t_name);}} */
  else if (strcmp(vname,"GLOBAL_T_n_cols")==0){ scanf("%d",&GLOBAL_T_n_cols); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_T_n_cols);}}
  else if (strcmp(vname,"GLOBAL_T_n_cind")==0){ scanf("%[^,;]",GLOBAL_T_n_cind); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_T_n_cind); sprintf(GLOBAL_T_n_cind,"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_T_n_cind);}}
  else if (strcmp(vname,"GLOBAL_S_n_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_S_n_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_S_n_name_[nv]); sprintf(GLOBAL_S_n_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_S_n_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_S_n_name")==0){ scanf("%[^,;]",GLOBAL_S_n_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_S_n_name);}} */
  else if (strcmp(vname,"GLOBAL_S_t_name_")==0){ length=GLOBAL_NBINS; for (nv=0;nv<length;nv++){ scanf("%[^,;]",GLOBAL_S_t_name_[nv]); if (strcmp(GLOBAL_DIR_XPRE,"\0")){ sprintf(tmpchar,GLOBAL_S_t_name_[nv]); sprintf(GLOBAL_S_t_name_[nv],"%s/%s",GLOBAL_DIR_XPRE,tmpchar);} if (nv<length-1){ scanf("%c",comma_vs_semicolon);} if (verbose>0){ printf("%s[%d] read to be %s\n",vname,nv,GLOBAL_S_t_name_[nv]);}}}
  /* else if (strcmp(vname,"GLOBAL_S_t_name")==0){ scanf("%[^,;]",GLOBAL_S_t_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_S_t_name);}} */
  else if (strcmp(vname,"GLOBAL_DIR_XPRE")==0){ scanf("%[^,;]",GLOBAL_DIR_XPRE); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_DIR_XPRE);} /* else if (strcmp(vname,"GLOBAL_DIR_XPRE")==0){ } */}
  else if (strcmp(vname,"GLOBAL_DIR_BASE")==0){ scanf("%[^,;]",GLOBAL_DIR_BASE); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_DIR_BASE);} /* else if (strcmp(vname,"GLOBAL_DIR_BASE")==0){ } */}
  else if (strcmp(vname,"GLOBAL_out_name")==0){ scanf("%[^,;]",GLOBAL_out_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_out_name);} /* else if (strcmp(vname,"GLOBAL_out_name")==0){ } */}
  else if (strcmp(vname,"GLOBAL_scorebox_out_xdrop")==0){ scanf("%[^,;]",GLOBAL_scorebox_out_xdrop); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_scorebox_out_xdrop);} /* else if (strcmp(vname,"GLOBAL_scorebox_out_xdrop")==0){ } */}
  else if (strcmp(vname,"GLOBAL_scorebox_row_max")==0){ scanf("%d",&GLOBAL_scorebox_row_max); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_scorebox_row_max);}}
  else if (strcmp(vname,"GLOBAL_scorebox_row_num")==0){ scanf("%d",&GLOBAL_scorebox_row_num); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_scorebox_row_num);}}
  else if (strcmp(vname,"GLOBAL_scorebox_col_max")==0){ scanf("%d",&GLOBAL_scorebox_col_max); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_scorebox_col_max);}}
  else if (strcmp(vname,"GLOBAL_scorebox_col_num")==0){ scanf("%d",&GLOBAL_scorebox_col_num); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_scorebox_col_num);}}
  else if (strcmp(vname,"GLOBAL_DIR_NAME")==0){ scanf("%[^,;]",GLOBAL_DIR_NAME); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_DIR_NAME);} /* else if (strcmp(vname,"GLOBAL_DIR_NAME")==0){ } */}
  else if (strcmp(vname,"END")==0){ /* do nothing */ if (verbose>0){ printf("end of input reached\n");}}
/*   else if (strcmp(vname,"yy")==0){ scanf("%zz",&yy); if (verbose>0){ printf("%s read to be %zz\n",vname,yy);}} */
  else /* if anything else */{ printf(" %% Error! vname %s in updateglobals\n",vname); exit(RET_READ_FAIL); }
}


void read_input()
{
  /* This reads the piped input file, or standard input.
     the variable names should not be longer than 128 characters */
  int verbose=GLOBAL_verbose;
  char vname[128],equals[128],space[128],semicolon[128];
  do{
    scanf("%[^=]",vname);scanf("%s",equals);scanf("%c",space);updateglobals(vname);scanf("%c",semicolon);
    if (verbose>1){ printf("At this point variable name is (%s), equals is (%s), semicolon is (%s)\n",vname,equals,semicolon);} 
    scanf("%c",semicolon);}
  while (strcmp(vname,"END")!=0);
}
