%%%%%%%%;
% data-set at: ;
% https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc33k ;
%%%%%%%%;

clear;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

flag_recalc = 0;
flag_replot = 0;
tolerance_master = 1e-2;
nf=0;

dir_data = sprintf('/%s/rangan/dir_bcc/dir_10xgenomics/dir_pbmc33k',string_root);
if (~exist(sprintf('%s_mat',dir_data),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_data)); mkdir(sprintf('%s_mat',dir_data)); end;
if (~exist(sprintf('%s_jpg',dir_data),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_data)); mkdir(sprintf('%s_jpg',dir_data)); end;

fname_tsne = sprintf('%s/analysis/tsne/projection.csv',dir_data);
fp = fopen(fname_tsne,'r'); tmp_ = textscan(fp,'%s%s%s','headerlines',1,'Delimiter',','); fclose(fp);
loading_tsne_ = [ cellfun(@str2num,tmp_{2}) , cellfun(@str2num,tmp_{3}) ];
label_tsne_ = tmp_{1};
