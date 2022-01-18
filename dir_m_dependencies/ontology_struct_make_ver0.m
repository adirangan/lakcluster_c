function ontology = ontology_struct_make_ver0(fname_ontology);

verbose=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% First we load an ontology;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%fname_ontology = '/data/rangan/dir_bcc/dir_seek_072916/data/human/go_bp_iea.txt';
fcheck(fname_ontology);
fid=fopen(fname_ontology);
n_line = 0; flag_continue=1;
flag_continue=1;
while (flag_continue); tmp_ = fgetl(fid); if tmp_~=-1; n_line = n_line+1; else flag_continue=0; end; end;
fclose(fid);
if (verbose>0); disp(sprintf(' %% reading %s: found %d lines',fname_ontology,n_line)); end;
n_pathway = floor(n_line/2);
pathway_name_ = cell(n_pathway,1);
pathway_EZid__ = cell(n_pathway,1);
pathway_size_ = zeros(n_pathway,1);
fid=fopen(fname_ontology);
for npathway=1:n_pathway;
tmp_ = fgetl(fid);
pathway_name_{npathway} = tmp_;
tmp_ = fgetl(fid);
pathway_EZid__{npathway} = str2num(tmp_);
pathway_size_(npathway) = length(pathway_EZid__{npathway});
end;%for npathway=1:n_pathway;
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% If we like we can calculate the pairwise intersection (i.e., gene-overlap) between pathways in the ontology. ;
% This takes a while, and is not all that informative; the summary is that ;
% there are several pathways that are large and strongly overlap with one another, ;
% but also many pathways which are small and rather disparate. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
flag_cap=0;
if flag_cap;
disp(sprintf(' %% calculating pairwise intersection between pathways'));
pathway_cap_ = zeros(n_pathway,n_pathway);
for np1=1:n_pathway; 
pathway_cap_(np1,np1) = 1.0;
if (mod(np1,100)==0); disp(sprintf(' %% np1 %d/%d',np1,n_pathway)); end;
for np2=np1+1:n_pathway;
tmp = length(intersect(pathway_EZid__{np1},pathway_EZid__{np2}))/max(1,min(pathway_size_(np1),pathway_size_(np2)));
pathway_cap_(np1,np2) = tmp; pathway_cap_(np2,np1) = tmp;
end;%for np2=1:n_pathway;
end;%for np1=1:n_pathway; 
end;%if flag_cap;

ontology = struct(...
		  'verbose',0 ...
		  ,'fname_ontology',fname_ontology...
		  ,'n_line',n_line...
		  ,'n_pathway',n_pathway...
		  );
ontology.pathway_name_ = pathway_name_;
ontology.pathway_EZid__ = pathway_EZid__;
ontology.pathway_size_ = pathway_size_;

