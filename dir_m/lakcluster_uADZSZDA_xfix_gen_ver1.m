function [output_string,gamma_d] = lakcluster_uADZSZDA_xfix_gen_ver1(rev_flag,A_n_rind_,Z_n_rind_,T_n_cind,GLOBAL_TEST_sparse,gamma,B_MLT,Ireq,shuffle_num);

nbins=length(A_n_rind_);
Z_bother = 0;
for nb1=0:nbins-1;
if (length(Z_n_rind_{1+nb1})>0); Z_bother = 1; end;
end;%for nb1=0:nbins-1;
if ~Z_bother; rev_flag = 0; end;

if Z_bother & rev_flag==1; rev_str = '_X'; end;
if Z_bother & rev_flag==0; rev_str = '_D'; end;
if ~Z_bother; rev_str = ''; end;
gamma_d = floor(gamma*1000); gamma_str = sprintf('_g%.3d',gamma_d);
if nbins>1; u_str = 'u'; else; u_str = ''; end;
if GLOBAL_TEST_sparse==1; D_str = 'D'; else; D_str = ''; end;
if length(T_n_cind)>1; S_str = 'S'; else; S_str = ''; end;
if Z_bother; Z_str = 'Z'; else; Z_str = 'A'; end;
type_str = sprintf('_%sA%s%s%s%s%sA',u_str,D_str,Z_str,S_str,Z_str,D_str);
if strcmp(type_str,'_AAAA'); type_str = ''; end;
B_str = sprintf('_B%.2d',B_MLT);
if nbins>1; I_str = sprintf('_I%d',Ireq); else; I_str = ''; end;
if shuffle_num>0; shuffle_str = sprintf('_s%.4d',shuffle_num); else; shuffle_str = ''; end;
output_string = sprintf('lak%s%s%s%s%s%s',type_str,rev_str,gamma_str,B_str,I_str,shuffle_str);
