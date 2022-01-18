function lisa = lisa_struct_make_ver0(mr_string,mc_string,cl_num,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;

lisa = struct(...
	      'verbose',0 ...
	      ,'mr_string',mr_string...
	      ,'mc_string',mc_string...
	      ,'cl_num',cl_num...
	      ,'flag_dex_vs_lak',flag_dex_vs_lak...
	      ,'gamma',gamma...
	      ,'B_MLT',B_MLT...
	      ,'n_mds',n_mds...
	      ,'flag_reverse',flag_reverse...
	      ,'n_maf',n_maf...
	      ,'n_cov',n_cov...
	      ,'n_scramble',n_scramble...
	      ,'n_shuffle',n_shuffle...
	      );
