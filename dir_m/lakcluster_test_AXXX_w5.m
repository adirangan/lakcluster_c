function lakcluster_test_AXXX_w5(test_string);
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
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/lakcluster_test_AXXX_w5.m ./dir_m/lakcluster_test_AXXX_w5.m ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/hpc_launcher.m /scratch/avr209/dir_m/hpc_launcher.m ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/slurm_launcher.m /scratch/avr209/dir_m/slurm_launcher.m ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/lakcluster_test_AXXX_w5.m /scratch/avr209/dir_m/lakcluster_test_AXXX_w5.m ;
  %}

nbins = 1;

N_ra = 2.^[6:13];
X_factor_ra = [0.45:0.025:0.65];
X_esm_ra = 10.^[-3:0.5:1];
gamma_ra = [0.01,0.05,0.15,0.50];%gamma_ra = [0.00,0.05,0.15,0.50];
B_ra = [0];%B_ra = [0,4,8,16];
rng_num_ra = 1:16;
X_factor_d_ra = floor(100*X_factor_ra);
X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra));
gamma_d_ra = floor(10*gamma_ra);
ver_num=1;

if strcmp(test_string,'lakcluster_test_AAAA_sym'); 
nbins = 1; N_ra = 2.^[6:13]; X_factor_ra = [0.40:0.0125:0.60]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=0; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_AAAA_sym'); 
if strcmp(test_string,'lakcluster_test_AAAA'); 
nbins = 1; N_ra = 2.^[6:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=4; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_AAAA'); 
if strcmp(test_string,'lakcluster_test_AAAA_gamma'); 
test_string = 'lakcluster_test_AAAA';
nbins = 1; N_ra = 2.^[10]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = 0.5.^[0:7]; B_ra = [0];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=4; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_AAAA_gamma'); 
if strcmp(test_string,'lakcluster_test_AZZA'); 
nbins = 1; N_ra = 2.^[6:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=1; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_AZZA'); 
if strcmp(test_string,'lakcluster_test_AZWY'); 
nbins = 1; N_ra = 2.^[6:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=1; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_AZWY'); 
if strcmp(test_string,'lakcluster_test_AAYY'); nbins = 1; end;
if strcmp(test_string,'lakcluster_test_ADAADA'); 
nbins = 1; N_ra = 2.^[6:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=0; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_ADAADA'); 
if strcmp(test_string,'lakcluster_test_AATAA'); 
nbins = 1; N_ra = 2.^[6:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=2; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_AATAA');
if strcmp(test_string,'lakcluster_test_uAAAA'); 
nbins = 2; N_ra = 2.^[6:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=0; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_uAAAA'); 
if strcmp(test_string,'lakcluster_test_uAZZA'); 
nbins = 2; N_ra = 2.^[6:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=2; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_uAZZA'); 
if strcmp(test_string,'lakcluster_test_uADZSZDA'); 
nbins = 3; N_ra = 2.^[6:2:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:64; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=6; 
end;%if strcmp(test_string,'lakcluster_test_uADZSZDA');
if strcmp(test_string,'lakcluster_test_uAZSZA'); 
nbins = 3; N_ra = 2.^[9:12]; X_factor_ra = [0.55:0.025:0.75]; X_esm_ra = 10.^[-3:0.5:+1]; gamma_ra = [0.05,0.50];%gamma_ra = [0.00,0.05,0.15,0.50]; 
B_ra = [1,4,8,16]; rng_num_ra = 1:16; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=1; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'lakcluster_test_uAZSZA'); 
if strcmp(test_string,'aucward_test_uAZSZA');
nbins = 1; N_ra = 2.^[6:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:128; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=1; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'aucward_test_uAZSZA');
if strcmp(test_string,'dexcluster_test_uAZSZA');
nbins = 1; N_ra = 2.^[6:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [0];
rng_num_ra = 1:128; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=4; disp(sprintf('%% Warning! check version number for %s',test_string)); return;
end;%if strcmp(test_string,'dexcluster_test_uAZSZA');
if strcmp(test_string,'dexcluster_test_uADZSZDA');
%nbins = 3; N_ra = 2.^[6:2:13]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
nbins = 3; N_ra = 2.^[10]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:128; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=1;
end;%if strcmp(test_string,'dexcluster_test_uADZSZDA');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rerunning the more complex simulations studies for paper: 081817 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq3'); 
nbins = 3; N_ra = 2.^[11]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:32; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=7; 
end;%if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq2'); 
nbins = 3; N_ra = 2.^[11]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:32; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=7; 
end;%if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(test_string,'lakcluster_test_uADZSZDA_Ireq1'); 
nbins = 3; N_ra = 2.^[11]; X_factor_ra = [0.45:0.0125:0.65]; X_esm_ra = 10.^[-3:0.25:1]; gamma_ra = [0.01]; B_ra = [32];
rng_num_ra = 1:32; X_factor_d_ra = floor(100*X_factor_ra); X_esm_d_ra = 10 + floor(-10*log10(X_esm_ra)); gamma_d_ra = floor(10*gamma_ra); ver_num=7; 
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
disp(sprintf(' %% currently finished %d/%d (%.2d%%)',sum(tot_n_dat(:)),lN*lf*le*lg*lB*lr,floor(100*sum(tot_n_dat(:))/(lN*lf*le*lg*lB*lr))));

pt_num=1; cmd_lim = 128; command_string = ''; ncommands=0; 

for ng=1:length(gamma_ra); gamma = gamma_ra(ng);
for nB=1:length(B_ra); B_MLT = B_ra(nB);
for nN=1:length(N_ra); N = N_ra(nN);
rng_len = length(rng_num_ra);
if 0;
elseif log2(N)>=13; rng_len = ceil(length(rng_num_ra)/4) ; 
elseif log2(N)>=12; rng_len = ceil(length(rng_num_ra)/2) ; 
elseif log2(N)>=11; rng_len = ceil(length(rng_num_ra)/2) ; 
elseif log2(N)>=10; rng_len = ceil(length(rng_num_ra)/1) ; 
elseif log2(N)>= 9; rng_len = ceil(length(rng_num_ra)/1) ; 
elseif log2(N)>= 8; rng_len = ceil(length(rng_num_ra)/1) ; 
end;%if;
rng_num_ra_use = 1:rng_len;
for nr=1:length(rng_num_ra_use); rng_num = rng_num_ra_use(nr);
for nf=1:length(X_factor_ra); X_factor = X_factor_ra(nf);
for ne=1:length(X_esm_ra); X_esm = X_esm_ra(ne);
disp(sprintf(' checking: nN %d nf %d ne %d ng %d nB %d nr %d',[nN,nf,ne,ng,nB,nr]));
if tot_n_dat(nN,nf,ne,ng,nB,nr)==1;
disp(sprintf(' deleting... ')); 
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
 else; 
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



