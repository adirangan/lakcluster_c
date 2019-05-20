	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->Y_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_TYnWtS_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb2]->M_Zn->rpop_j - 0*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_Wn->cpop_j);
	  dra_plusdivequals_s___m_m(&(E->QR_TYnWtS_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_TYnWtS->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_TYnWtS_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_TYnWtS_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->Y_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->Y_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_TYnYtT_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb2]->M_An->rpop_j - 1*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_Yn->cpop_j);
	  dra_plusdivequals_s___m_m(&(E->QR_TYnYtT_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_TYnYtT->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_TYnYtT_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_TYnYtT_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->Y_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->A_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_TAnZtS_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb2]->M_Zn->rpop_j - 0*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_Zn->cpop_j);
	  dra_plusdivequals_s___m_m(&(E->QR_TAnZtS_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_TAnZtS->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_TAnZtS_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_TAnZtS_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->A_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->A_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb2]->M_An->rpop_j - 1*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_An->cpop_j);
	  dra_plusdivequals_s___m_m(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_TAnAtT->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_TAnAtT_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->A_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->Z_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_TYnWtS_nrm[nbx*D->Y_ncols*D->T_ncols])),D->Y_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb1]->M_An->rpop_j) * (long long int)(E_[nb2]->M_Wn->rpop_j);
	  dra_plusdivequals_s___m_m(&(D->QC_TYnWtS_nrm[nbx*D->Y_ncols*D->T_ncols]),D->Y_ncols,D->T_ncols,F->QC_TYnWtS->lf,(double)lld,D->Y_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_TYnWtS_nrm[nbx*D->Y_ncols*D->T_ncols]),"double_trn",D->Y_ncols,D->T_ncols," %%%% D->QC_TYnWtS_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->Z_rbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->A_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_TYnYtT_nrm[nbx*D->Y_ncols*D->T_ncols])),D->Y_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb1]->M_An->rpop_j) * (long long int)(E_[nb2]->M_An->rpop_j - 1*(nb1==nb2?1:0));
	  dra_plusdivequals_s___m_m(&(D->QC_TYnYtT_nrm[nbx*D->Y_ncols*D->T_ncols]),D->Y_ncols,D->T_ncols,F->QC_TYnYtT->lf,(double)lld,D->Y_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_TYnYtT_nrm[nbx*D->Y_ncols*D->T_ncols]),"double_trn",D->Y_ncols,D->T_ncols," %%%% D->QC_TYnYtT_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->A_rbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->Z_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_TAnZtS_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb1]->M_An->rpop_j) * (long long int)(E_[nb2]->M_Zn->rpop_j);
	  dra_plusdivequals_s___m_m(&(D->QC_TAnZtS_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_TAnZtS->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_TAnZtS_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_TAnZtS_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->Z_rbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->A_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb1]->M_An->rpop_j) * (long long int)(E_[nb2]->M_An->rpop_j - 1*(nb1==nb2?1:0));
	  dra_plusdivequals_s___m_m(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_TAnAtT->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_TAnAtT_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->A_rbother){ } */}
