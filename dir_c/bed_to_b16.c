#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

unsigned long long int fsize_0(char *fname)
{
  unsigned long long s=0;
  FILE *fp=NULL;
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning, could not open %s.\n",fname); exit(RET_READ_FAIL);}
  fseek(fp, 0l, SEEK_END); s = (unsigned long long int) ftell(fp); rewind(fp);
  return s;
}

unsigned long int wc_0(char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  unsigned long int n_line=0;
  char *line=NULL;
  size_t n_tmp=0;
  ssize_t n_read=0;
  unsigned long long int n_char=0;
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning, could not open %s.\n",fname); exit(RET_READ_FAIL);}
  while ((n_read = getline(&line,&n_tmp,fp)) != -1){ n_char += n_read; n_line += 1;} free(line);
  fclose(fp); fp=NULL;
  if (verbose){ printf(" %% n_line %d n_char %lld\n",n_line,n_char);}
  return n_line;
}

void bed_to_b16_test()
{
  int verbose=1;
  int MDA_d_[2];
  char fname_bed[1024];
  char fname_bim[1024];
  char fname_fam[1024];
  FILE *fp_bim=NULL;
  FILE *fp_bed=NULL;
  FILE *fp_fam=NULL;
  double snp_mss_threshold = 0.01;
  double snp_maf_threshold = 0.10;
  double pat_mss_threshold = 0.04;
  double snp_I_threshold = 0.01;
  unsigned long int n_bim = 999;//unsigned long int n_bim = 11342;
  //unsigned long int n_bim = 265;
  unsigned long int n_fam = 488377;
  unsigned long int n_bam = 0;
  unsigned long long int n_bed = 0;
  unsigned long long int size_bed = 0;
  int *flip_allele_=NULL; /* denotes whether a snp reads the first allele as minor and the second as major (+1) or vice-versa (-1) */
  unsigned char flip_value_[4];
  unsigned long long int *snp__=NULL,*pat__=NULL,*tot_=NULL;
  unsigned char *bed_line_=NULL;
  unsigned char bed_byte=0,bed_value=0;
  unsigned char key0=0;
  unsigned char key1=0;
  unsigned char key2=0;
  unsigned char ub=0;
  int nbim=0,nfam=0,nbam=0,n_read=0;
  int snp_and_tag=3;
  int snp_xor_tag=2;
  int snp_nor_tag=0;
  int snp_mss_tag=1;
  int n_val = 4;
  int snp_tab=0,pat_tab=0;
  int *snp_mode_=NULL;
  int *pat_mode_=NULL;
  double *snp_frq_and_=NULL,*snp_frq_xor_=NULL,*snp_frq_nor_=NULL,*snp_frq_mss_=NULL,*snp_p_opt_=NULL,*snp_q_opt_=NULL,*snp_I_opt_=NULL;
  unsigned long long int *snp_tot_dat_=NULL;
  double *pat_frq_and_=NULL,*pat_frq_xor_=NULL,*pat_frq_nor_=NULL,*pat_frq_mss_=NULL;
  unsigned long long int *pat_tot_dat_=NULL;
  double tmp_min=0,tmp_and=0,tmp_xor=0,tmp_nor=0;
  int *index_snp_retain_from_bim_=NULL;
  int *index_snp_bim_from_retain_=NULL;
  int nsnp_retain=0,n_snp_retain=0;
  int *index_pat_retain_from_fam_=NULL;
  int *index_pat_fam_from_retain_=NULL;
  int npat_retain=0,n_pat_retain=0;
  int n_allele=0;
  int *index_allele_from_allele_=NULL;
  int *index_snp_retain_from_allele_=NULL;
  int *index_allele_from_snp_retain__=NULL;
  double *allele_frequency_=NULL;
  double *allele_maf_=NULL;
  double *allele_mss_=NULL;
  double *allele_I_opt_=NULL;
  double tmp_maf=0;
  int *index_allele_orig_from_sort_=NULL;
  int *index_allele_sort_from_orig_=NULL;
  int *i_workspace_=NULL;
  double *d_workspace_=NULL;
  if (verbose){ printf(" %% [entering bed_to_b16_test]\n");}
  sprintf(fname_bed,"/home/rangan/dir_bcc/dir_ukb/calls/ukb_cal_chr21_v2.bed");
  sprintf(fname_bim,"/home/rangan/dir_bcc/dir_ukb/bim/ukb_snp_chr21_v2.bim");
  sprintf(fname_fam,"/home/rangan/dir_bcc/dir_ukb/fam_calls/ukb43036_cal_chr21_v2_s488288.fam");
  //sprintf(fname_bed,"/home/rangan/dir_bcc/dir_ukb/calls/ukb_cal_chrMT_v2.bed");
  //sprintf(fname_bim,"/home/rangan/dir_bcc/dir_ukb/bim/ukb_snp_chrMT_v2.bim");
  //sprintf(fname_fam,"/home/rangan/dir_bcc/dir_ukb/fam_calls/ukb43036_cal_chrMT_v2_s488288.fam");
  if (verbose){ printf(" %% bim: %s -> wc %d\n",fname_bim,wc_0(fname_bim));}
  if (verbose){ printf(" %% fam: %s -> wc %d\n",fname_fam,wc_0(fname_fam));}
  if ((fp_bed=fopen(fname_bed,"r"))==NULL){ printf(" %% Warning, could not open %s.\n",fname_bed); exit(RET_READ_FAIL);}
  if ((fp_bim=fopen(fname_bim,"r"))==NULL){ printf(" %% Warning, could not open %s.\n",fname_bim); exit(RET_READ_FAIL);}
  if ((fp_fam=fopen(fname_fam,"r"))==NULL){ printf(" %% Warning, could not open %s.\n",fname_fam); exit(RET_READ_FAIL);}
  fread(&key0,sizeof(unsigned char),1,fp_bed); if (verbose){ printf(" %% bed: key0 read %.2x\n",key0);}
  if (key0!=108){ printf(" %% Warning, key0 %d!=108\n",key0);}
  fread(&key1,sizeof(unsigned char),1,fp_bed); if (verbose){ printf(" %% bed: key1 read %.2x\n",key1);}
  if (key1!= 27){ printf(" %% Warning, key1 %d!= 27\n",key1);}
  fread(&key2,sizeof(unsigned char),1,fp_bed); if (verbose){ printf(" %% bed: key2 read %.2x\n",key2);}
  if (key2!=  1){ printf(" %% Warning, key2 %d!=  1\n",key2);}
  n_bam = n_fam/4 + ( n_fam%4==0 ? 0 : 1);
  n_bed = 3 + n_bim*n_bam;
  size_bed = fsize_0(fname_bed);
  if (verbose){ printf(" %% n_bed %lld vs size_bed %lld --> difference %lld\n",n_bed,size_bed,n_bed-size_bed);}
  if (size_bed<n_bed){ printf(" %% Warning, size_bed %lld < n_bed %lld\n",size_bed,n_bed);}
  flip_allele_ = (int *)wkspace_all0c(sizeof(int)*n_bim);
  for (nbim=0;nbim<n_bim;nbim++){ flip_allele_[nbim] = +1; /* for (nbim=0;nbim<n_bim;nbim++){ } */}
  snp__ = (unsigned long long int *)wkspace_all0c(sizeof(unsigned long long int)*n_val*n_bim);
  pat__ = (unsigned long long int *)wkspace_all0c(sizeof(unsigned long long int)*n_val*n_val*n_bam);
  tot_ = (unsigned long long int *)wkspace_all0c(sizeof(unsigned long long int)*n_val);
  bed_line_ = (unsigned char *)wkspace_all0c(sizeof(unsigned char)*n_bam);
  /* 
     Here the tags are: 
     00 --> 0 = 'nor' = homozygous (minor)
     01 --> 1 = 'mss' = missing
     10 --> 2 = 'xor' = heterozygous
     11 --> 3 = 'and' = homozygous (major)
     Note that if the snp is to be 'flipped', then nor should count as and, and vice versa.
     However, xor and mss should remain unchanged.
     To accomplish this we simply use a table lookup.
   */
  flip_value_[0] = (unsigned char)3; flip_value_[1] = (unsigned char)1; flip_value_[2] = (unsigned char)2; flip_value_[3] = (unsigned char)0;
  snp_tab=0;
  for (nbim=0;nbim<n_bim;nbim++){
    if ((nbim%1000)==0){ printf(" %% nbim %d/%d\n",nbim,n_bim);}
    n_read = fread(&(bed_line_[0]),sizeof(unsigned char),n_bam,fp_bed); if (n_read!=n_bam){ printf(" %% Warning, n_read %d < n_bam %d\n",n_read,n_bam);}
    if (flip_allele_[nbim]==+1){
    nfam=0; pat_tab=0;
    for (nbam=0;nbam<n_bam;nbam++){
      bed_byte= (unsigned char)(bed_line_[nbam]);
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 0) << 6) >> 6);
      snp__[snp_tab+bed_value]++; pat__[pat_tab+bed_value]++; tot_[bed_value]++; nfam+=1; pat_tab+=n_val;
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 2) << 6) >> 6);
      snp__[snp_tab+bed_value]++; pat__[pat_tab+bed_value]++; tot_[bed_value]++; nfam+=1; pat_tab+=n_val;
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 4) << 6) >> 6);
      snp__[snp_tab+bed_value]++; pat__[pat_tab+bed_value]++; tot_[bed_value]++; nfam+=1; pat_tab+=n_val;
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 6) << 0) >> 0);
      snp__[snp_tab+bed_value]++; pat__[pat_tab+bed_value]++; tot_[bed_value]++; nfam+=1; pat_tab+=n_val;
      /* for (nbam=0;nbam<n_bam;nbam++){ } */}
    /* if (flip_allele_[nbim]==+1){ } */}
    if (flip_allele_[nbim]==0){ /* do nothing */ }
    if (flip_allele_[nbim]==-1){
    nfam=0; pat_tab=0;
    for (nbam=0;nbam<n_bam;nbam++){
      bed_byte= (unsigned char)(bed_line_[nbam]);
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 0) << 6) >> 6)];
      snp__[snp_tab+bed_value]++; pat__[pat_tab+bed_value]++; tot_[bed_value]++; nfam+=1; pat_tab+=n_val;
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 2) << 6) >> 6)];
      snp__[snp_tab+bed_value]++; pat__[pat_tab+bed_value]++; tot_[bed_value]++; nfam+=1; pat_tab+=n_val;
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 4) << 6) >> 6)];
      snp__[snp_tab+bed_value]++; pat__[pat_tab+bed_value]++; tot_[bed_value]++; nfam+=1; pat_tab+=n_val;
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 6) << 0) >> 0)];
      snp__[snp_tab+bed_value]++; pat__[pat_tab+bed_value]++; tot_[bed_value]++; nfam+=1; pat_tab+=n_val;
      /* for (nbam=0;nbam<n_bam;nbam++){ } */}
    /* if (flip_allele_[nbim]==-1){ } */}
    snp_tab+=n_val;
    /* for (nbim=0;nbim<n_bim;nbim++){ } */}
  raprintf(tot_,"unsigned long long int",1,n_val," %% tot_: ");
  fclose(fp_bed); fp_bed = NULL;
  fclose(fp_bim); fp_bim = NULL;
  fclose(fp_fam); fp_fam = NULL;
  if (verbose>0){
    MDA_d_[0] = n_val; MDA_d_[1] = n_bim; MDA_write_ulli(2,MDA_d_,snp__,"/home/rangan/dir_bcc/dir_ukb/dir_mda/snp__.mda");
    MDA_d_[0] = n_val; MDA_d_[1] = n_val*n_bam; MDA_write_ulli(2,MDA_d_,pat__,"/home/rangan/dir_bcc/dir_ukb/dir_mda/pat__.mda");
    MDA_d_[0] = n_val; MDA_d_[1] = 1; MDA_write_ulli(2,MDA_d_,tot_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/tot_.mda");
    /* if (verbose>0){ } */}
  snp_mode_ = (int *)wkspace_all0c(sizeof(int)*n_bim);
  snp_tab=0;
  for (nbim=0;nbim<n_bim;nbim++){
    bed_value=1; if (snp__[snp_tab+bed_value]>snp__[snp_tab+snp_mode_[nbim]]){ snp_mode_[nbim]=bed_value;};
    bed_value=2; if (snp__[snp_tab+bed_value]>snp__[snp_tab+snp_mode_[nbim]]){ snp_mode_[nbim]=bed_value;};
    bed_value=3; if (snp__[snp_tab+bed_value]>snp__[snp_tab+snp_mode_[nbim]]){ snp_mode_[nbim]=bed_value;};
    snp_tab+=n_val;
    /* for (nbim=0;nbim<n_bim;nbim++){ } */}
  pat_mode_ = (int *)wkspace_all0c(sizeof(int)*n_fam);
  pat_tab=0;
  for (nfam=0;nfam<n_fam;nfam++){
    bed_value=1; if (pat__[pat_tab+bed_value]>pat__[pat_tab+pat_mode_[nfam]]){ pat_mode_[nfam]=bed_value;};
    bed_value=2; if (pat__[pat_tab+bed_value]>pat__[pat_tab+pat_mode_[nfam]]){ pat_mode_[nfam]=bed_value;};
    bed_value=3; if (pat__[pat_tab+bed_value]>pat__[pat_tab+pat_mode_[nfam]]){ pat_mode_[nfam]=bed_value;};
    pat_tab+=n_val;
    /* for (nfam=0;nfam<n_fam;nfam++){ } */}
  if (verbose>0){
    MDA_d_[0] = n_bim; MDA_d_[1] = 1; MDA_write_i4(2,MDA_d_,snp_mode_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/snp_mode_.mda");
    MDA_d_[0] = n_fam; MDA_d_[1] = 1; MDA_write_i4(2,MDA_d_,pat_mode_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/pat_mode_.mda");
    /* if (verbose>0){ } */}
  /*
    Calculating relative-entropy (i.e., KL-divergence) for each snp within each study.
    If we assume that a snp has minor allele frequency q and major allele frequency p,
    then the distribution of (n_and,n_xor,n_nor) should be proportional to:
    (p*p, 2*p*q, q*q), respectively.
    If the true distribution is given by (frq_and,frq_xor,frq_nor) = (n_and,n_xor,n_nor)/n_patient_keep,
    then the relative entropy between the true distribution and the expected distribution is:
    I = frq_and*log(frq_and/p^2) + frq_xor*log(frq_xor/(2*p*q)) + frq_nor*log(frq_nor/q^2).
    The derivative of I with respect to p is:
    dI/dp = -2*frq_and/p - (q-p)*frq_xor/(p*q) + 2*frq_nor/q.
    Setting this to 0, we see that the optimal p and q are given by:
    2*p*frq_nor - 2*q*frq_and = (q-p)*frq_xor,
    p*(2*frq_nor+2*frq_and+2*frq_xor) = 2*frq_and + frq_xor,
    or simply:
    p_opt = frq_and + 0.5*frq_xor, '));
    q_opt = frq_nor + 0.5*frq_xor. ,'));
    I_opt = frq_and*log(frq_and/p_opt^2) + frq_xor*log(frq_xor/(2*p_opt*q_opt)) + frq_nor*log(frq_nor/q_opt^2).
  */
  snp_frq_and_=(double *)wkspace_all0c(sizeof(double)*n_bim);
  snp_frq_xor_=(double *)wkspace_all0c(sizeof(double)*n_bim);
  snp_frq_nor_=(double *)wkspace_all0c(sizeof(double)*n_bim);
  snp_frq_mss_=(double *)wkspace_all0c(sizeof(double)*n_bim);
  snp_tot_dat_=(unsigned long long int *)wkspace_all0c(sizeof(unsigned long long int)*n_bim);
  snp_p_opt_=(double *)wkspace_all0c(sizeof(double)*n_bim);
  snp_q_opt_=(double *)wkspace_all0c(sizeof(double)*n_bim);
  snp_I_opt_=(double *)wkspace_all0c(sizeof(double)*n_bim);
  snp_tab=0;
  for (nbim=0;nbim<n_bim;nbim++){
    snp_tot_dat_[nbim] = snp__[snp_tab+snp_and_tag]+snp__[snp_tab+snp_xor_tag]+snp__[snp_tab+snp_nor_tag];
    tmp_min = (double)1/(double)snp_tot_dat_[nbim];
    snp_frq_mss_[nbim] = (double)snp__[snp_tab+snp_mss_tag]/(double)n_fam;
    snp_frq_and_[nbim] = (double)snp__[snp_tab+snp_and_tag]/(double)snp_tot_dat_[nbim];
    snp_frq_xor_[nbim] = (double)snp__[snp_tab+snp_xor_tag]/(double)snp_tot_dat_[nbim];
    snp_frq_nor_[nbim] = (double)snp__[snp_tab+snp_nor_tag]/(double)snp_tot_dat_[nbim];
    snp_p_opt_[nbim] = snp_frq_and_[nbim] + 0.5*snp_frq_xor_[nbim];
    snp_q_opt_[nbim] = snp_frq_nor_[nbim] + 0.5*snp_frq_xor_[nbim];
    tmp_and = snp_frq_and_[nbim]*log(snp_frq_and_[nbim]/maximum(tmp_min,pow(snp_p_opt_[nbim],2)));
    tmp_xor = snp_frq_xor_[nbim]*log(snp_frq_xor_[nbim]/maximum(tmp_min,(2*snp_p_opt_[nbim]*snp_q_opt_[nbim])));
    tmp_nor = snp_frq_nor_[nbim]*log(snp_frq_nor_[nbim]/maximum(tmp_min,pow(snp_q_opt_[nbim],2)));
    snp_I_opt_[nbim] = tmp_and + tmp_xor + tmp_nor;
    snp_tab+=n_val;
    /* for (nbim=0;nbim<n_bim;nbim++){ } */}
    /* find mode */
  unsigned char *snp_tag_mode_=NULL;
  snp_tag_mode_ = (unsigned char *)wkspace_all0c(sizeof(unsigned char)*n_bim);
  snp_tab=0;
  for (nbim=0;nbim<n_bim;nbim++){
    snp_tag_mode_[nbim] = snp_xor_tag;
    if ( (snp__[snp_tab+snp_xor_tag]>=snp__[snp_tab+snp_and_tag]) && (snp__[snp_tab+snp_xor_tag]>=snp__[snp_tab+snp_nor_tag]) ){ snp_tag_mode_[nbim] = snp_xor_tag;}
    if ( (snp__[snp_tab+snp_nor_tag]>=snp__[snp_tab+snp_and_tag]) && (snp__[snp_tab+snp_nor_tag]>=snp__[snp_tab+snp_xor_tag]) ){ snp_tag_mode_[nbim] = snp_nor_tag;}
    if ( (snp__[snp_tab+snp_and_tag]>=snp__[snp_tab+snp_xor_tag]) && (snp__[snp_tab+snp_and_tag]>=snp__[snp_tab+snp_nor_tag]) ){ snp_tag_mode_[nbim] = snp_and_tag;}
    snp_tab+=n_val;
    /* for (nbim=0;nbim<n_bim;nbim++){ } */}
  if (verbose>0){
    MDA_d_[0] = n_bim; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,snp_frq_mss_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/snp_frq_mss_.mda");
    MDA_d_[0] = n_bim; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,snp_frq_and_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/snp_frq_and_.mda");
    MDA_d_[0] = n_bim; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,snp_frq_xor_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/snp_frq_xor_.mda");
    MDA_d_[0] = n_bim; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,snp_frq_nor_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/snp_frq_nor_.mda");
    MDA_d_[0] = n_bim; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,snp_p_opt_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/snp_p_opt_.mda");
    MDA_d_[0] = n_bim; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,snp_q_opt_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/snp_q_opt_.mda");
    MDA_d_[0] = n_bim; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,snp_I_opt_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/snp_I_opt_.mda");
    /* if (verbose>0){ } */}
  pat_frq_and_=(double *)wkspace_all0c(sizeof(double)*n_fam);
  pat_frq_xor_=(double *)wkspace_all0c(sizeof(double)*n_fam);
  pat_frq_nor_=(double *)wkspace_all0c(sizeof(double)*n_fam);
  pat_frq_mss_=(double *)wkspace_all0c(sizeof(double)*n_fam);
  pat_tot_dat_=(unsigned long long int *)wkspace_all0c(sizeof(unsigned long long int)*n_fam);
  pat_tab=0;
  for (nfam=0;nfam<n_fam;nfam++){
    pat_tot_dat_[nfam] = pat__[pat_tab+snp_and_tag]+pat__[pat_tab+snp_xor_tag]+pat__[pat_tab+snp_nor_tag];
    tmp_min = (double)1/(double)pat_tot_dat_[nfam];
    pat_frq_mss_[nfam] = (double)pat__[pat_tab+snp_mss_tag]/(double)n_bim;
    pat_frq_and_[nfam] = (double)pat__[pat_tab+snp_and_tag]/(double)pat_tot_dat_[nfam];
    pat_frq_xor_[nfam] = (double)pat__[pat_tab+snp_xor_tag]/(double)pat_tot_dat_[nfam];
    pat_frq_nor_[nfam] = (double)pat__[pat_tab+snp_nor_tag]/(double)pat_tot_dat_[nfam];
    pat_tab+=n_val;
    /* for (nfam=0;nfam<n_fam;nfam++){ } */}
  if (verbose>0){
    MDA_d_[0] = n_fam; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,pat_frq_mss_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/pat_frq_mss_.mda");
    MDA_d_[0] = n_fam; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,pat_frq_and_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/pat_frq_and_.mda");
    MDA_d_[0] = n_fam; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,pat_frq_xor_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/pat_frq_xor_.mda");
    MDA_d_[0] = n_fam; MDA_d_[1] = 1; MDA_write_r8(2,MDA_d_,pat_frq_nor_,"/home/rangan/dir_bcc/dir_ukb/dir_mda/pat_frq_nor_.mda");
    /* if (verbose>0){ } */}
  /* find snp indices that satisfy thresholds */
  index_snp_retain_from_bim_ = (int *)wkspace_all0c(sizeof(int)*n_bim);
  index_snp_bim_from_retain_ = (int *)wkspace_all0c(sizeof(int)*n_bim);
  nsnp_retain = 0;
  for (nbim=0;nbim<n_bim;nbim++){
    if ( (snp_tot_dat_[nbim]>0) && (snp_frq_and_[nbim]>0) && (snp_frq_xor_[nbim]>0) && (snp_frq_nor_[nbim]>0) && (minimum(snp_p_opt_[nbim],snp_q_opt_[nbim])>=snp_maf_threshold) && (snp_frq_mss_[nbim]<=snp_mss_threshold) && (snp_I_opt_[nbim]<=snp_I_threshold) ){ index_snp_bim_from_retain_[nsnp_retain] = nbim; index_snp_retain_from_bim_[nbim] = nsnp_retain; nsnp_retain++;}
    /* for (nbim=0;nbim<n_bim;nbim++){ } */}
  n_snp_retain = nsnp_retain;
  /* find pat indices that satisfy thresholds */
  index_pat_retain_from_fam_ = (int *)wkspace_all0c(sizeof(int)*n_fam);
  index_pat_fam_from_retain_ = (int *)wkspace_all0c(sizeof(int)*n_fam);
  npat_retain = 0;
  for (nfam=0;nfam<n_fam;nfam++){
    if ( (pat_tot_dat_[nfam]>0) && (pat_frq_mss_[nfam]<=pat_mss_threshold) ){ index_pat_fam_from_retain_[npat_retain] = nfam; index_pat_retain_from_fam_[nfam] = npat_retain; npat_retain++;}
    /* for (nfam=0;nfam<n_fam;nfam++){ } */}
  n_pat_retain = npat_retain;
  /* extract alleles (nor,xor,and) from retained snps */
  /* alleles are stored so that allele-type varies quickly, and snp-number varies slowly */
  int allele_nor_tag = 0; int allele_xor_tag = 1; int allele_and_tag = 2;
  n_allele = 3*n_snp_retain;
  index_allele_from_allele_ = (int *)wkspace_all0c(sizeof(int)*n_allele);
  index_snp_retain_from_allele_ = (int *)wkspace_all0c(sizeof(int)*n_allele);
  index_allele_from_snp_retain__ = (int *)wkspace_all0c(sizeof(int)*n_allele); /* Contains a stack of 3 arrays, each of length n_snp_retain. The arrays correspond (respectively) to nor, xor and and. */
  allele_frequency_ = (double *)wkspace_all0c(sizeof(double)*n_allele);
  allele_maf_ = (double *)wkspace_all0c(sizeof(double)*n_allele);
  allele_mss_ = (double *)wkspace_all0c(sizeof(double)*n_allele);
  allele_I_opt_ = (double *)wkspace_all0c(sizeof(double)*n_allele);
  for (nsnp_retain=0;nsnp_retain<n_snp_retain;nsnp_retain++){
    nbim = index_snp_bim_from_retain_[nsnp_retain];
    index_allele_from_allele_[allele_nor_tag + nsnp_retain*3] = allele_nor_tag + nsnp_retain*3;
    index_allele_from_allele_[allele_xor_tag + nsnp_retain*3] = allele_xor_tag + nsnp_retain*3;
    index_allele_from_allele_[allele_and_tag + nsnp_retain*3] = allele_and_tag + nsnp_retain*3;
    index_snp_retain_from_allele_[allele_nor_tag + nsnp_retain*3] = nsnp_retain;
    index_snp_retain_from_allele_[allele_xor_tag + nsnp_retain*3] = nsnp_retain;
    index_snp_retain_from_allele_[allele_and_tag + nsnp_retain*3] = nsnp_retain;
    index_allele_from_snp_retain__[nsnp_retain + allele_nor_tag*n_snp_retain] = allele_nor_tag;
    index_allele_from_snp_retain__[nsnp_retain + allele_xor_tag*n_snp_retain] = allele_xor_tag;
    index_allele_from_snp_retain__[nsnp_retain + allele_and_tag*n_snp_retain] = allele_and_tag;
    allele_frequency_[allele_nor_tag + nsnp_retain*3] = snp_frq_nor_[nbim];
    allele_frequency_[allele_xor_tag + nsnp_retain*3] = snp_frq_xor_[nbim];
    allele_frequency_[allele_and_tag + nsnp_retain*3] = snp_frq_and_[nbim];
    tmp_maf = minimum(snp_p_opt_[nbim],snp_q_opt_[nbim]);
    allele_maf_[allele_nor_tag + nsnp_retain*3] = tmp_maf;
    allele_maf_[allele_xor_tag + nsnp_retain*3] = tmp_maf;
    allele_maf_[allele_and_tag + nsnp_retain*3] = tmp_maf;
    allele_mss_[allele_nor_tag + nsnp_retain*3] = snp_frq_mss_[nbim];
    allele_mss_[allele_xor_tag + nsnp_retain*3] = snp_frq_mss_[nbim];
    allele_mss_[allele_and_tag + nsnp_retain*3] = snp_frq_mss_[nbim];
    allele_I_opt_[allele_nor_tag + nsnp_retain*3] = snp_I_opt_[nbim];
    allele_I_opt_[allele_xor_tag + nsnp_retain*3] = snp_I_opt_[nbim];
    allele_I_opt_[allele_and_tag + nsnp_retain*3] = snp_I_opt_[nbim];
    /* for (nsnp_retain=0;nsnp_retain<n_snp_retain;nsnp_retain++){ } */}
  /* sort the alleles by frequency */
  index_allele_orig_from_sort_ = (int *)wkspace_all0c(sizeof(int)*n_allele);
  index_allele_sort_from_orig_ = (int *)wkspace_all0c(sizeof(int)*n_allele);
  i_workspace_ = (int *)wkspace_all0c(sizeof(int)*n_allele);
  d_workspace_ = (double *)wkspace_all0c(sizeof(double)*n_allele);
  dQuickSort_index_index_driver(n_allele,allele_frequency_,1,d_workspace_,index_allele_orig_from_sort_,i_workspace_,index_allele_sort_from_orig_);
  /* build b16 array so that patients vary quickly and alleles vary slowly */
  if (verbose){ printf(" %% n_snp_retain %d; n_pat_retain %d; n_allele %d;\n",n_snp_retain,n_pat_retain,n_allele);}
  int bitj=16; int bit8=8;
  int n_allele_extend=0,l_allele=0;
  int n_pat_retain_extend=0,l_pat_retain=0;
  bitj = 16; bit8 = 8;
  n_allele_extend = (bitj - (n_allele%bitj))%bitj;
  l_allele = (n_allele + n_allele_extend)/bit8;
  n_pat_retain_extend = (bitj - (n_pat_retain%bitj))%bitj;
  l_pat_retain = (n_pat_retain + n_pat_retain_extend)/bit8;
  if (verbose){ printf(" %% size_A_n: %.2fGB\n",(double)l_pat_retain*(double)n_allele/1e9);}
  if (verbose){ printf(" %% size_A_t: %.2fGB\n",(double)l_allele*(double)n_pat_retain/1e9);}
  unsigned char *A_n_=NULL;
  A_n_ = (unsigned char *)wkspace_all0c(sizeof(unsigned char)*(size_t)l_pat_retain*(size_t)n_allele);
  unsigned char *snp_tag_line_=NULL;
  snp_tag_line_ = (unsigned char *)wkspace_all0c(sizeof(unsigned char)*n_val*n_bam);
  unsigned char *allele_nor_line_=NULL;
  unsigned char *allele_xor_line_=NULL;
  unsigned char *allele_and_line_=NULL;
  allele_nor_line_ = (unsigned char *)wkspace_all0c(sizeof(unsigned char)*l_pat_retain);
  allele_xor_line_ = (unsigned char *)wkspace_all0c(sizeof(unsigned char)*l_pat_retain);
  allele_and_line_ = (unsigned char *)wkspace_all0c(sizeof(unsigned char)*l_pat_retain);
  int nallele_nor_orig=0,nallele_xor_orig=0,nallele_and_orig=0;
  int nallele_nor_sort=0,nallele_xor_sort=0,nallele_and_sort=0;
  int nbim_target=0,nfam_target=0;
  int lpat_retain=0;
  unsigned long long int ltab=0;
  unsigned char tmp_pat_add=0;
  unsigned char tmp_a=0,tmp_b=0;
  unsigned char *b_=NULL;
  if ((fp_bed=fopen(fname_bed,"r"))==NULL){ printf(" %% Warning, could not open %s.\n",fname_bed); exit(RET_READ_FAIL);}
  fread(&key0,sizeof(unsigned char),1,fp_bed); if (verbose){ printf(" %% bed: key0 read %.2x\n",key0);}
  if (key0!=108){ printf(" %% Warning, key0 %d!=108\n",key0);}
  fread(&key1,sizeof(unsigned char),1,fp_bed); if (verbose){ printf(" %% bed: key1 read %.2x\n",key1);}
  if (key1!= 27){ printf(" %% Warning, key1 %d!= 27\n",key1);}
  fread(&key2,sizeof(unsigned char),1,fp_bed); if (verbose){ printf(" %% bed: key2 read %.2x\n",key2);}
  if (key2!=  1){ printf(" %% Warning, key2 %d!=  1\n",key2);}
  snp_tab=0; nsnp_retain=0; nbim_target = index_snp_bim_from_retain_[nsnp_retain];
  for (nbim=0;nbim<n_bim;nbim++){
    if ((nbim%1000)==0){ printf(" %% nbim %d/%d\n",nbim,n_bim);}
    n_read = fread(&(bed_line_[0]),sizeof(unsigned char),n_bam,fp_bed); if (n_read!=n_bam){ printf(" %% Warning, n_read %d < n_bam %d\n",n_read,n_bam);}
    if (flip_allele_[nbim]==+1){
    nfam=0;
    for (nbam=0;nbam<n_bam;nbam++){
      bed_byte= (unsigned char)(bed_line_[nbam]);
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 0) << 6) >> 6);
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 2) << 6) >> 6);
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 4) << 6) >> 6);
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 6) << 0) >> 0);
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      /* for (nbam=0;nbam<n_bam;nbam++){ } */}
    /* if (flip_allele_[nbim]==+1){ } */}
    if (flip_allele_[nbim]==0){ /* do nothing */ }
    if (flip_allele_[nbim]==-1){
    nfam=0;
    for (nbam=0;nbam<n_bam;nbam++){
      bed_byte= (unsigned char)(bed_line_[nbam]);
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 0) << 6) >> 6)];
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 2) << 6) >> 6)];
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 4) << 6) >> 6)];
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 6) << 0) >> 0)];
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      /* for (nbam=0;nbam<n_bam;nbam++){ } */}
    /* if (flip_allele_[nbim]==-1){ } */}
    if (nbim==nbim_target){
      if (index_snp_retain_from_bim_[nbim]!=nsnp_retain){ printf(" %% Warning, nbim %d index_snp_retain_from_bim_[nbim] %d nsnp_retain %d\n",(int)nbim,(int)index_snp_retain_from_bim_[nbim],(int)nsnp_retain);}
      nallele_nor_orig = index_allele_from_snp_retain__[nsnp_retain + allele_nor_tag*n_snp_retain] + nsnp_retain*3;
      nallele_xor_orig = index_allele_from_snp_retain__[nsnp_retain + allele_xor_tag*n_snp_retain] + nsnp_retain*3;
      nallele_and_orig = index_allele_from_snp_retain__[nsnp_retain + allele_and_tag*n_snp_retain] + nsnp_retain*3;
      nallele_nor_sort = index_allele_sort_from_orig_[nallele_nor_orig];
      nallele_xor_sort = index_allele_sort_from_orig_[nallele_xor_orig];
      nallele_and_sort = index_allele_sort_from_orig_[nallele_and_orig];
      for (lpat_retain=0;lpat_retain<l_pat_retain;lpat_retain++){ allele_nor_line_[lpat_retain]=0;}
      for (lpat_retain=0;lpat_retain<l_pat_retain;lpat_retain++){ allele_xor_line_[lpat_retain]=0;}
      for (lpat_retain=0;lpat_retain<l_pat_retain;lpat_retain++){ allele_and_line_[lpat_retain]=0;}
      npat_retain=0; nfam_target=index_pat_fam_from_retain_[npat_retain];
      for (nfam=0;nfam<n_fam;nfam++){
	if (nfam==nfam_target){
	  lpat_retain = npat_retain/bit8;
	  tmp_pat_add = ((unsigned char)1) << (7-(npat_retain%bit8));
	  allele_nor_line_[lpat_retain] = ( snp_tag_line_[nfam]==snp_nor_tag ? tmp_pat_add : 0 );
	  ltab=(unsigned long long int)lpat_retain+(unsigned long long int)nallele_nor_sort*(unsigned long long int)l_pat_retain; /* if (allele_nor_line_[lpat_retain]){ b_ = &(A_n_[nallele_nor_sort*l_pat_retain]); bset__on(b_,npat_retain);} */ A_n_[ltab]+=allele_nor_line_[lpat_retain];
	  //if (ltab==20195500 && allele_nor_line_[lpat_retain]){ printf(" %% nor: flipping 20195500 bit %d: nbim %.6d nsnp_retain %.6d --> orig (%.6d,%.6d,%.6d) --> sort (%.6d,%.6d,%.6d) nfam %.6d npat_retain %.6d\n",npat_retain%bit8,nbim,nsnp_retain,nallele_nor_orig,nallele_xor_orig,nallele_and_orig,nallele_nor_sort,nallele_xor_sort,nallele_and_sort,nfam,npat_retain);}
	  allele_xor_line_[lpat_retain] = ( snp_tag_line_[nfam]==snp_xor_tag ? tmp_pat_add : 0 );
	  ltab=(unsigned long long int)lpat_retain+(unsigned long long int)nallele_xor_sort*(unsigned long long int)l_pat_retain; /* if (allele_xor_line_[lpat_retain]){ b_ = &(A_n_[nallele_xor_sort*l_pat_retain]); bset__on(b_,npat_retain);} */ A_n_[ltab]+=allele_xor_line_[lpat_retain];
	  //if (ltab==20195500 && allele_xor_line_[lpat_retain]){ printf(" %% xor: flipping 20195500 bit %d\n",npat_retain%bit8);}
	  allele_and_line_[lpat_retain] = ( snp_tag_line_[nfam]==snp_and_tag ? tmp_pat_add : 0 );
	  ltab=(unsigned long long int)lpat_retain+(unsigned long long int)nallele_and_sort*(unsigned long long int)l_pat_retain; /* if (allele_and_line_[lpat_retain]){ b_ = &(A_n_[nallele_and_sort*l_pat_retain]); bset__on(b_,npat_retain);} */ A_n_[ltab]+=allele_and_line_[lpat_retain];
	  //if (ltab==20195500 && allele_and_line_[lpat_retain]){ printf(" %% and: flipping 20195500 bit %d\n",npat_retain%bit8);}
	  if (verbose>1){
	    if ( ( (nsnp_retain==0) && (npat_retain==119498) ) || ( (nsnp_retain==0) && (npat_retain==0) )){
	      printf(" %% nsnp_retain %.6d npat_retain %.6d\n",nsnp_retain,npat_retain);
	      printf(" %% nbim %.6d nbim_target %.6d\n",nbim,nbim_target);
	      printf(" %% nfam %.6d nfam_target %.6d\n",nfam,nfam_target);
	      printf(" %% nallele_nor_orig %.6d nallele_xor_orig %.6d nallele_and_orig %.6d\n",nallele_nor_orig,nallele_xor_orig,nallele_and_orig);
	      printf(" %% nallele_nor_sort %.6d nallele_xor_sort %.6d nallele_and_sort %.6d\n",nallele_nor_sort,nallele_xor_sort,nallele_and_sort);
	      printf(" %% lpat_retain %d tmp_pat_add %d\n",lpat_retain,tmp_pat_add);
	      printf(" %% allele_nor_line_[lpat_retain] %.3d allele_xor_line_[lpat_retain] %.3d allele_and_line_[lpat_retain] %.3d\n",(int)allele_nor_line_[lpat_retain],(int)allele_xor_line_[lpat_retain],(int)allele_and_line_[lpat_retain]);
	      ltab=(unsigned long long int)lpat_retain+(unsigned long long int)nallele_nor_sort*(unsigned long long int)l_pat_retain; printf(" %% nor: ltab %.6lld A_n_[ltab] %.3d --> %d\n",ltab,A_n_[ltab], ( (unsigned char)(A_n_[ltab]) >> (7-(npat_retain%bit8)) ) & (unsigned char)1 );
	      ltab=(unsigned long long int)lpat_retain+(unsigned long long int)nallele_xor_sort*(unsigned long long int)l_pat_retain; printf(" %% xor: ltab %.6lld A_n_[ltab] %.3d --> %d\n",ltab,A_n_[ltab], ( (unsigned char)(A_n_[ltab]) >> (7-(npat_retain%bit8)) ) & (unsigned char)1 );
	      ltab=(unsigned long long int)lpat_retain+(unsigned long long int)nallele_and_sort*(unsigned long long int)l_pat_retain; printf(" %% and: ltab %.6lld A_n_[ltab] %.3d --> %d\n",ltab,A_n_[ltab], ( (unsigned char)(A_n_[ltab]) >> (7-(npat_retain%bit8)) ) & (unsigned char)1 );
	      /* if ( (nsnp_retain==0) && (npat_retain==0) ){ } */}
	    /* if (verbose>1){ } */}
	  npat_retain += 1; nfam_target=index_pat_fam_from_retain_[npat_retain];
	  /* if (nfam==nfam_target){ } */}
	/* for (nfam=0;nfam<n_fam;nfam++){ } */}
      if (npat_retain!=n_pat_retain){ printf(" %% Warning, npat_retain %d n_pat_retain %d\n",npat_retain,n_pat_retain);}
      //ltab=(unsigned long long int)l_pat_(unsigned long long int)retain*nallele_nor_(unsigned long long int)sort; for (lpat_retain=0;lpat_retain<l_pat_retain;lpat_retain++){ A_n_[ltab] += allele_nor_line_[lpat_retain]; ltab+=1; /* for (lpat_retain=0;lpat_retain<l_pat_retain;lpat_retain++){ } */}
      //ltab=(unsigned long long int)l_pat_(unsigned long long int)retain*nallele_xor_(unsigned long long int)sort; for (lpat_retain=0;lpat_retain<l_pat_retain;lpat_retain++){ A_n_[ltab] += allele_xor_line_[lpat_retain]; ltab+=1; /* for (lpat_retain=0;lpat_retain<l_pat_retain;lpat_retain++){ } */}
      //ltab=(unsigned long long int)l_pat_(unsigned long long int)retain*nallele_and_(unsigned long long int)sort; for (lpat_retain=0;lpat_retain<l_pat_retain;lpat_retain++){ A_n_[ltab] += allele_and_line_[lpat_retain]; ltab+=1; /* for (lpat_retain=0;lpat_retain<l_pat_retain;lpat_retain++){ } */}
      nsnp_retain += 1; nbim_target = index_snp_bim_from_retain_[nsnp_retain];
      /* if (nbim==nbim_target){ } */}
    snp_tab+=n_val;
    /* for (nbim=0;nbim<n_bim;nbim++){ } */}
  if (nsnp_retain!=n_snp_retain){ printf(" %% Warning, nsnp_retain %d n_snp_retain %d\n",nsnp_retain,n_snp_retain);}
  fclose(fp_bed); fp_bed = NULL;
  /* test a few entries */
  int n_snp_test=0,nsnp_test=0;
  int n_pat_test=0,npat_test=0;
  n_snp_test = 4; n_pat_test = 5;
  for (nsnp_test=0;nsnp_test<n_snp_test;nsnp_test++){
    nsnp_retain = floor((double)(n_snp_retain-1)*(double)nsnp_test/(double)(n_snp_test-1));
    if (nsnp_test==n_snp_test-1){ nsnp_retain=n_snp_retain-1;}
    npat_retain = floor((double)(n_pat_retain-1)*(double)nsnp_test/(double)(n_snp_test-1));
    if (nsnp_test==n_snp_test-1){ npat_retain=n_pat_retain-1;}
    nbim_target = index_snp_bim_from_retain_[nsnp_retain];
    nfam_target = index_pat_fam_from_retain_[npat_retain];
    /* for (nsnp_test=0;nsnp_test<n_snp_test;nsnp_test++){ } */}
  int nsnp_retain_target=0,npat_retain_target=0;
  if ((fp_bed=fopen(fname_bed,"r"))==NULL){ printf(" %% Warning, could not open %s.\n",fname_bed); exit(RET_READ_FAIL);}
  fread(&key0,sizeof(unsigned char),1,fp_bed); if (verbose){ printf(" %% bed: key0 read %.2x\n",key0);}
  if (key0!=108){ printf(" %% Warning, key0 %d!=108\n",key0);}
  fread(&key1,sizeof(unsigned char),1,fp_bed); if (verbose){ printf(" %% bed: key1 read %.2x\n",key1);}
  if (key1!= 27){ printf(" %% Warning, key1 %d!= 27\n",key1);}
  fread(&key2,sizeof(unsigned char),1,fp_bed); if (verbose){ printf(" %% bed: key2 read %.2x\n",key2);}
  if (key2!=  1){ printf(" %% Warning, key2 %d!=  1\n",key2);}
  snp_tab=0; nsnp_test=0;
  nsnp_retain_target = floor((double)(n_snp_retain-1)*(double)nsnp_test/(double)(n_snp_test-1));
  if (nsnp_test==n_snp_test-1){ nsnp_retain_target=n_snp_retain-1;}
  nbim_target = index_snp_bim_from_retain_[nsnp_retain_target];
  for (nbim=0;nbim<n_bim;nbim++){
    if ((nbim%1000)==0){ printf(" %% nbim %d/%d\n",nbim,n_bim);}
    n_read = fread(&(bed_line_[0]),sizeof(unsigned char),n_bam,fp_bed); if (n_read!=n_bam){ printf(" %% Warning, n_read %d < n_bam %d\n",n_read,n_bam);}
    if (flip_allele_[nbim]==+1){
    nfam=0;
    for (nbam=0;nbam<n_bam;nbam++){
      bed_byte= (unsigned char)(bed_line_[nbam]);
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 0) << 6) >> 6);
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 2) << 6) >> 6);
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 4) << 6) >> 6);
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = (unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 6) << 0) >> 0);
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      /* for (nbam=0;nbam<n_bam;nbam++){ } */}
    /* if (flip_allele_[nbim]==+1){ } */}
    if (flip_allele_[nbim]==0){ /* do nothing */ }
    if (flip_allele_[nbim]==-1){
    nfam=0;
    for (nbam=0;nbam<n_bam;nbam++){
      bed_byte= (unsigned char)(bed_line_[nbam]);
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 0) << 6) >> 6)];
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 2) << 6) >> 6)];
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 4) << 6) >> 6)];
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      bed_value = flip_value_[(unsigned char)((unsigned char)((unsigned char)((unsigned char)bed_byte >> 6) << 0) >> 0)];
      snp_tag_line_[nfam] = ( bed_value==snp_mss_tag ? snp_tag_mode_[nbim] : bed_value ); nfam+=1;
      /* for (nbam=0;nbam<n_bam;nbam++){ } */}
    /* if (flip_allele_[nbim]==-1){ } */}
    if (nbim==nbim_target){
      if (index_snp_retain_from_bim_[nbim]!=nsnp_retain_target){ printf(" %% Warning, nbim %d index_snp_retain_from_bim_[nbim] %d nsnp_retain_target %d\n",(int)nbim,(int)index_snp_retain_from_bim_[nbim],(int)nsnp_retain_target);}
      nallele_nor_orig = index_allele_from_snp_retain__[nsnp_retain_target + allele_nor_tag*n_snp_retain] + nsnp_retain_target*3;
      nallele_xor_orig = index_allele_from_snp_retain__[nsnp_retain_target + allele_xor_tag*n_snp_retain] + nsnp_retain_target*3;
      nallele_and_orig = index_allele_from_snp_retain__[nsnp_retain_target + allele_and_tag*n_snp_retain] + nsnp_retain_target*3;
      nallele_nor_sort = index_allele_sort_from_orig_[nallele_nor_orig];
      nallele_xor_sort = index_allele_sort_from_orig_[nallele_xor_orig];
      nallele_and_sort = index_allele_sort_from_orig_[nallele_and_orig];
      printf(" %% nsnp_test %.6d nsnp_retain_target %.6d nbim %0.6d nbim_target %0.6d\n",nsnp_test,nsnp_retain_target,nbim,nbim_target);
      printf(" %% nallele_nor_orig %.6d nallele_xor_orig %.6d nallele_and_orig %.6d\n",nallele_nor_orig,nallele_xor_orig,nallele_and_orig);
      printf(" %% nallele_nor_sort %.6d nallele_xor_sort %.6d nallele_and_sort %.6d\n",nallele_nor_sort,nallele_xor_sort,nallele_and_sort);
      npat_test=0;
      npat_retain_target = floor((double)(n_pat_retain-1)*(double)npat_test/(double)(n_pat_test-1));
      if (npat_test==n_pat_test-1){ npat_retain_target=n_pat_retain-1;}
      nfam_target = index_pat_fam_from_retain_[npat_retain_target];
      for (nfam=0;nfam<n_fam;nfam++){
	if (nfam==nfam_target){
	  lpat_retain = npat_retain_target/bit8;
	  printf(" %% %% npat_retain_target %.6d lpat_retain %.6d\n",npat_retain_target,lpat_retain);
	  allele_nor_line_[lpat_retain] = ( snp_tag_line_[nfam_target]==snp_nor_tag ? (unsigned char)1 : (unsigned char)0 );
	  allele_xor_line_[lpat_retain] = ( snp_tag_line_[nfam_target]==snp_xor_tag ? (unsigned char)1 : (unsigned char)0 );
	  allele_and_line_[lpat_retain] = ( snp_tag_line_[nfam_target]==snp_and_tag ? (unsigned char)1 : (unsigned char)0 );
	  ltab=(unsigned long long int)lpat_retain+(unsigned long long int)nallele_nor_sort*(unsigned long long int)l_pat_retain; tmp_a = allele_nor_line_[lpat_retain]; tmp_b = ( (unsigned char)(A_n_[ltab]) >> (7-(npat_retain_target%bit8)) ) & (unsigned char)1;
	  printf(" %% %% nsnp_test %.6d npat_test %.6d nbim_target %.6d nfam_target %.6d ; ltab %0.6lld A_n_[ltab] %.3d ; nor: tmp_a %.3d tmp_b %.3d\n",nsnp_test,npat_test,nbim_target,nfam_target,ltab,(int)A_n_[ltab],(int)tmp_a,(int)tmp_b);
	  ltab=(unsigned long long int)lpat_retain+(unsigned long long int)nallele_xor_sort*(unsigned long long int)l_pat_retain; tmp_a = allele_xor_line_[lpat_retain]; tmp_b = ( (unsigned char)(A_n_[ltab]) >> (7-(npat_retain_target%bit8)) ) & (unsigned char)1;
	  printf(" %% %% nsnp_test %.6d npat_test %.6d nbim_target %.6d nfam_target %.6d ; ltab %0.6lld A_n_[ltab] %.3d ; xor: tmp_a %.3d tmp_b %.3d\n",nsnp_test,npat_test,nbim_target,nfam_target,ltab,(int)A_n_[ltab],(int)tmp_a,(int)tmp_b);
	  ltab=(unsigned long long int)lpat_retain+(unsigned long long int)nallele_and_sort*(unsigned long long int)l_pat_retain; tmp_a = allele_and_line_[lpat_retain]; tmp_b = ( (unsigned char)(A_n_[ltab]) >> (7-(npat_retain_target%bit8)) ) & (unsigned char)1;
	  printf(" %% %% nsnp_test %.6d npat_test %.6d nbim_target %.6d nfam_target %.6d ; ltab %0.6lld A_n_[ltab] %.3d ; and: tmp_a %.3d tmp_b %.3d\n",nsnp_test,npat_test,nbim_target,nfam_target,ltab,(int)A_n_[ltab],(int)tmp_a,(int)tmp_b);
	  printf(" %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% \n");
	  npat_test += 1;
	  npat_retain_target = floor((double)(n_pat_retain-1)*(double)npat_test/(double)(n_pat_test-1));
	  if (npat_test==n_pat_test-1){ npat_retain_target=n_pat_retain-1;}
	  nfam_target = index_pat_fam_from_retain_[npat_retain_target];
	  /* if (nfam==nfam_target){ } */}
	/* for (nfam=0;nfam<n_fam;nfam++){ } */}
      if (npat_test!=n_pat_test){ printf(" %% Warning, npat_test %d n_pat_test %d\n",npat_test,n_pat_test);}
      nsnp_test += 1;
      nsnp_retain_target = floor((double)(n_snp_retain-1)*(double)nsnp_test/(double)(n_snp_test-1));
      if (nsnp_test==n_snp_test-1){ nsnp_retain_target=n_snp_retain-1;}
      nbim_target = index_snp_bim_from_retain_[nsnp_retain_target];
      /* if (nbim==nbim_target){ } */}
    snp_tab+=n_val;
    /* for (nbim=0;nbim<n_bim;nbim++){ } */}
  if (nsnp_test!=n_snp_test){ printf(" %% Warning, nsnp_test %d n_snp_test %d\n",nsnp_test,n_snp_test);}
  fclose(fp_bed); fp_bed = NULL;
  wkspace_printf();
  if (verbose){ printf(" %% [finished bed_to_b16_test]\n");}
}

