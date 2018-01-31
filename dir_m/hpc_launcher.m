function hpc_launcher(command_string,nnodes,ppn,memdecl,walltime);
% launches an instance of matlab set to execute 'command_string';
% nnodes (number of nodes -- integer);
% ppn (processors per node -- integer);
% walltime (runtim in hours -- integer);

disp(command_string);

dir_hpc_launcher = sprintf('%s/dir_hpc_launcher',pwd);
hpc_numfile = sprintf('%s/hpc_launcher_numfile.txt',dir_hpc_launcher);

nd=0;
if (~exist(dir_hpc_launcher,'dir')); system(sprintf('mkdir %s',dir_hpc_launcher)); end;
if (~exist(hpc_numfile,'file'));
while (exist(sprintf('%s/hpc_launcher_%.4d.pbs',dir_hpc_launcher,nd),'file')); nd = nd+1; end;%while;
disp(sprintf('searching... name %s/hpc_launcher_%.4d.pbs',dir_hpc_launcher,nd));
elseif exist(hpc_numfile,'file');
fid = fopen(hpc_numfile,'r'); nd_find = textscan(fid,'%d',1); fclose(fid);
nd = nd_find{1}+1;
disp(sprintf('choosing: name %s/hpc_launcher_%.4d.pbs',dir_hpc_launcher,nd));
end; % if exist hpc_numfile ;
fid = fopen(hpc_numfile,'w'); fprintf(fid,'%d',nd); fclose(fid);

filename_m = sprintf('%s/hpc_launcher_%.4d.m',dir_hpc_launcher,nd);
disp(sprintf('trying to open %s',filename_m));
fid_m = fopen(filename_m,'w'); fprintf(fid_m,'%s;\n',command_string); fclose(fid_m);
filename_pbs = sprintf('%s/hpc_launcher_%.4d.pbs',dir_hpc_launcher,nd);
fid_pbs = fopen(filename_pbs,'w'); 
fprintf(fid_pbs,'#!/bin/bash \n');
fprintf(fid_pbs,'#PBS -N hpc_launcher_%.4d \n',nd);
fprintf(fid_pbs,'#PBS -l nodes=%d:ppn=%d \n',nnodes,ppn);
fprintf(fid_pbs,'#PBS -l mem=%dGb \n',memdecl);
fprintf(fid_pbs,'#PBS -l walltime=%d:59:59 \n',walltime);
filename_out = sprintf('%s/hpc_launcher_%.4d.out',dir_hpc_launcher,nd);
fprintf(fid_pbs,'#PBS -o %s \n',filename_out);
filename_err = sprintf('%s/hpc_launcher_%.4d.err',dir_hpc_launcher,nd);
fprintf(fid_pbs,'#PBS -e %s \n',filename_err);
fprintf(fid_pbs,'#PBS -j oe \n');
fprintf(fid_pbs,'cd %s \n',pwd);
%fprintf(fid_pbs,'module load matlab/2014a \n',filename_m,filename_out);
fprintf(fid_pbs,'module load matlab \n',filename_m,filename_out);
fprintf(fid_pbs,'matlab -nodisplay < %s > %s  \n',filename_m,filename_out);
fclose(fid_pbs);
filename_log = sprintf('%s/hpc_launcher_log.txt',dir_hpc_launcher);
fid_log = fopen(filename_log,'a'); fprintf(fid_log,'%s: %s; nn %d ppn %d walltime %d \n',filename_m,command_string,nnodes,ppn,walltime); fclose(fid_log);
system(sprintf('qsub %s',filename_pbs));

