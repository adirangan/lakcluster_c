function [ij_] = localmax(Z_);
% find all local maxima of Z_;
n_Z = numel(Z_);
[~,ij_] = max(Z_);
for nZ=2:n_Z-1;
if is_imax(Z_,nZ); ij_ = [ij_;nZ]; end;
end;%for nZ=2:n_Z-1;
ij_ = unique(ij_);
[~,p_] = sort(Z_(ij_),'descend');
ij_ = ij_(p_);

function flag_decrease_and = is_imax(Z_,ij);
% determines whether or not index ij is an interior maximum of Z_. ;
verbose=0;
n_Z = numel(Z_);
%%%%%%%%;
flag_decrease_neg = 0;
tab = ij-1; continue_flag=(tab>0);
while (continue_flag);
if (continue_flag==1 & Z_(tab)<Z_(ij)); flag_decrease_neg = 1; continue_flag=0; end;
if (continue_flag==1 & Z_(tab)>Z_(ij)); flag_decrease_neg = 0; continue_flag=0; end;
if (continue_flag==1 & Z_(tab)==Z_(ij)); tab=tab-1; continue_flag=1; end;
if (tab<=0); continue_flag=0; end;
end;%while;
if (verbose); disp(sprintf(' %% ij %d tab %d neg %d',ij,tab,flag_decrease_neg)); end;
%%%%%%%%;
flag_decrease_pos = 0;
tab = ij+1; continue_flag=(tab<n_Z);
while (continue_flag);
if (continue_flag==1 & Z_(tab)<Z_(ij)); flag_decrease_pos = 1; continue_flag=0; end;
if (continue_flag==1 & Z_(tab)>Z_(ij)); flag_decrease_pos = 0; continue_flag=0; end;
if (continue_flag==1 & Z_(tab)==Z_(ij)); tab=tab+1; continue_flag=1; end;
if (tab>=n_Z); continue_flag=0; end;
end;%while;
if (verbose); disp(sprintf(' %% ij %d tab %d pos %d',ij,tab,flag_decrease_pos)); end;
%%%%%%%%;
flag_decrease_and = flag_decrease_neg & flag_decrease_pos;