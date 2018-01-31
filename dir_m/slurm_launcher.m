function slurm_launcher(command_string,nnodes,tpn,memdecl,walltime);
% launches an instance of matlab set to execute 'command_string';
% nnodes (number of nodes -- integer);
% tpn (tasks-per-node -- integer);
% walltime (runtim in hours -- integer);

disp(command_string);

dir_slurm_launcher = sprintf('%s/dir_slurm_launcher',pwd);
slurm_numfile = sprintf('%s/slurm_launcher_numfile.txt',dir_slurm_launcher);

nd=0;
if (~exist(dir_slurm_launcher,'dir')); system(sprintf('mkdir %s',dir_slurm_launcher)); end;
if (~exist(slurm_numfile,'file'));
while (exist(sprintf('%s/slurm_launcher_%.4d.sh',dir_slurm_launcher,nd),'file')); nd = nd+1; end;%while;
disp(sprintf('searching... name %s/slurm_launcher_%.4d.sh',dir_slurm_launcher,nd));
elseif exist(slurm_numfile,'file');
fid = fopen(slurm_numfile,'r'); nd_find = textscan(fid,'%d',1); fclose(fid);
nd = nd_find{1}+1;
disp(sprintf('choosing: name %s/slurm_launcher_%.4d.sh',dir_slurm_launcher,nd));
end; % if exist slurm_numfile ;
fid = fopen(slurm_numfile,'w'); fprintf(fid,'%d',nd); fclose(fid);

filename_infix = sprintf('slurm_launcher_%.4d',nd);
filename_m = sprintf('%s/slurm_launcher_%.4d.m',dir_slurm_launcher,nd);
disp(sprintf('trying to open %s',filename_m));
fid_m = fopen(filename_m,'w'); fprintf(fid_m,'%s\n',command_string); fclose(fid_m);
filename_sh = sprintf('%s/slurm_launcher_%.4d.sh',dir_slurm_launcher,nd);
fid_sh = fopen(filename_sh,'w'); 
fprintf(fid_sh,'#!/bin/sh \n');
fprintf(fid_sh,'#\n');
fprintf(fid_sh,'#SBATCH --verbose\n');
fprintf(fid_sh,'#SBATCH --job-name=%s\n',filename_infix);
filename_out = sprintf('%s/slurm_launcher_%.4d.out',dir_slurm_launcher,nd);
fprintf(fid_sh,'#SBATCH --output=%s\n',filename_out);
filename_err = sprintf('%s/slurm_launcher_%.4d.err',dir_slurm_launcher,nd);
fprintf(fid_sh,'#SBATCH --error=%s\n',filename_err);
fprintf(fid_sh,'#SBATCH --time=%d:00:00\n',walltime);
fprintf(fid_sh,'#SBATCH --nodes=%d --ntasks-per-node=%d\n',nnodes,tpn);
fprintf(fid_sh,'#SBATCH --mem=%dGB\n',memdecl);
fprintf(fid_sh,'\n');
fprintf(fid_sh,'/bin/hostname\n');
fprintf(fid_sh,'/bin/pwd\n');
fprintf(fid_sh,'cd /scratch/avr209/dir_m\n');
fprintf(fid_sh,'module load matlab/2016b\n');
fprintf(fid_sh,'matlab -nodisplay < %s > %s\n',filename_m,filename_out);
fclose(fid_sh);
filename_log = sprintf('%s/slurm_launcher_log.txt',dir_slurm_launcher);
fid_log = fopen(filename_log,'a'); fprintf(fid_log,'%s: %s; nn %d tpn %d walltime %d \n',filename_m,command_string,nnodes,tpn,walltime); fclose(fid_log);
system(sprintf('sbatch %s',filename_sh));

