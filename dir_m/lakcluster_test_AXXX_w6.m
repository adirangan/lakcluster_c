function lakcluster_test_AXXX_w6(test_string);
% designed to work with lakcluster_ver18 ;
%{
  % look for number of files with: ;
  find . -maxdepth 1 -type d -print0 | xargs -0 -I {} sh -c 'echo -e $(find {} | wc -l) {}' | sort -n

  % copies m-files from local to hpc ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/Makefile ./Makefile ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/lakcluster_ver18.c ./lakcluster_ver18.c ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/?*.[m] ./dir_m/ ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_c/?*.[ch] ./dir_c/ ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_in/?*.in ./dir_in/ ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/lakcluster_test_AXXX_w6.m ./dir_m/lakcluster_test_AXXX_w6.m ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/hpc_launcher.m /scratch/avr209/dir_m/hpc_launcher.m ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/slurm_launcher.m /scratch/avr209/dir_m/slurm_launcher.m ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/lakcluster_test_AXXX_w6.m /scratch/avr209/dir_m/lakcluster_test_AXXX_w6.m ;
  %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rerunning the more complex simulations studies for paper: 081817 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq3'); 
nbins = 3; N_ra = 2.^[11]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:8; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=7; 
end;%if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq2'); 
nbins = 3; N_ra = 2.^[11]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:8; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=7; 
end;%if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq1'); 
nbins = 3; N_ra = 2.^[11]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:8; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=7; 
end;%if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test_string,'dexcluster_test_uADZSZDA_Ireq3'); 
nbins = 3; N_ra = 2.^[11]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:32; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=7; 
end;%if strcmp(test_string,'dexcluster_test_uADZSZDA_Ireq3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test_string,'dexcluster_test_uADZSZDA_Ireq2'); 
nbins = 3; N_ra = 2.^[11]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:32; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=7; 
end;%if strcmp(test_string,'dexcluster_test_uADZSZDA_Ireq2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test_string,'dexcluster_test_uADZSZDA_Ireq1'); 
nbins = 3; N_ra = 2.^[11]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:32; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=7; 
end;%if strcmp(test_string,'dexcluster_test_uADZSZDA_Ireq1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_trunk = pwd;
dir_base = sprintf('%s/dir_%s',dir_trunk,test_string);
aucname = sprintf('%s/xx_out.txt',dir_base);
nnodes = 1; tpn = 8; memdecl = 18; walltime = 48;

lN = length(N_ra);
lf = length(X_factor_d_ra);
le = length(X_esm_d_ra);
lg = length(gamma_d_ra);
lB = length(B_ra);
lr = length(rng_num_ra);
tot_n_length = lN*lf*le*lg*lB*lr;
tot_n_dat = zeros(lN,lf,le,lg,lB,lr);
if exist(aucname,'file');
disp(sprintf(' %% found %s',aucname));
x = csvread(aucname);
for nx=1:size(x,1);
N = x(nx,1);
X_factor_d = x(nx,2);
X_esm_d = x(nx,3);
gamma_d = x(nx,4);
B_MLT = x(nx,5);
rng_num = x(nx,6);
bN = find(abs(N_ra-N)<=1);
bf = find(abs(X_factor_d_ra-X_factor_d)<=1);
be = find(abs(X_esm_d_ra-X_esm_d)<=1);
bg = find(abs(gamma_d_ra-gamma_d)<=1);
bB = find(abs(B_ra-B_MLT)<=1);
br = find(rng_num_ra==rng_num);
tot_n_dat(bN,bf,be,bg,bB,br) = 1;
end;%for nx=1:size(x,1);
end;% if exist(aucname,'file');
disp(sprintf(' %% currently finished %d/%d (%.2d%%)',sum(tot_n_dat(:)),tot_n_length,floor(100*sum(tot_n_dat(:))/(tot_n_length))));

pt_num=1; cmd_lim = 17*8; command_string = ''; ncommands=0; 

for ng=1:length(gamma_ra); gamma = gamma_ra(ng);
for nB=1:length(B_ra); B_MLT = B_ra(nB);
for nN=1:length(N_ra); N = N_ra(nN);
for nr=1:length(rng_num_ra); rng_num = rng_num_ra(nr);
for nf=1:length(X_factor_ra); X_factor = X_factor_ra(nf);
for ne=1:length(X_esm_ra); X_esm = X_esm_ra(ne);
disp(sprintf(' %% checking: nN %d nf %d ne %d ng %d nB %d nr %d',[nN,nf,ne,ng,nB,nr]));
if tot_n_dat(nN,nf,ne,ng,nB,nr)==1;
disp(sprintf(' %% found; deleting... ')); 
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
 else; 
disp(sprintf(' %% not yet found; running... ')); 
command_string = sprintf('%s\n%s_ver%d(''%s'',%d,%0.16f,%0.16f,%0.16f,%d,%d,%d);',command_string,test_string,ver_num,dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num); 
command_string = sprintf('%s xfix_del(''%s'',''%s'',%d,%0.16f,%0.16f,%0.16f,%d,%d,%d);',command_string,dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins); 
ncommands = ncommands + 1;
if ncommands>cmd_lim; slurm_launcher(command_string,nnodes,tpn,memdecl,walltime); pt_num = pt_num + 1; command_string=''; ncommands=0; end;
end;% if tot_n_dat(nN,nf,ne,ng,nB,nr)==1;
end;%for ne=1:length(X_esm_ra); X_esm = X_esm_ra(ne);
end;%for nf=1:length(X_factor_ra); X_factor = X_factor_ra(nf);
end;%for nr=1:length(rng_num_ra); rng_num = rng_num_ra(nr);
end;%for nN=1:length(N_ra); N = N_ra(nN);
end;%for nB=1:length(B_ra); B_MLT = B_ra(nB);
end;%for ng=1:length(gamma_ra); gamma = gamma_ra(ng);
if ncommands>0; slurm_launcher(command_string,nnodes,tpn,memdecl,walltime); pt_num = pt_num + 1; command_string=''; ncommands=0; end;



