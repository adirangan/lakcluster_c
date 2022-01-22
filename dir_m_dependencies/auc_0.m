function [auc] = auc_0(v_ctrl,v_case);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [auc] = auc_0(v_ctrl,v_case);
%
% This function gets the auc (area under the curve) of the
% roc curve (receiver-operating-characteristic curve) for data v_ctrl and v_case,
% where v_case is expected to typically take on higher values than v_ctrl.
% The inputs are v_ctrl and v_case, which will be treated as arrays.
% The output is the auc (a double between 0 and 1). 
%
% test by running with no arguments:
% i.e., >> auc_0();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_flag=0;

if (nargin<2);
b=0.05;n=4*1024; n=12340;
disp(sprintf(' testing auc_0: '));
disp(sprintf(' generating test data with %d points',n));
%v_ctrl = max(0,min(1,randn(n+2,1)/10+0.5-b));v_case = max(0,min(1,randn(n+3,1)/10+0.5+b));
%K=12; v_ctrl = round(K*max(0,min(1,randn(n-4,1)/10+0.5-b)))/K; v_case = round(K*max(0,min(1,randn(n+13,1)/10+0.5+b)))/K;
K=64; v_ctrl = round(K*max(0,min(1,randn(n-3,1)/10+0.5-b))); v_case = round(K*max(0,min(1,randn(n+2,1)/10+0.5+b)));
plot_flag=1;
end;%if (isempty(v_ctrl) | isempty(v_case));

if (isempty(v_ctrl) | isempty(v_case));
auc = 0.5;
else%if (isempty(v_ctrl) | isempty(v_case));
v_ctrl = reshape(v_ctrl,1,length(v_ctrl));
v_case = reshape(v_case,1,length(v_case));
buc = [v_ctrl , v_case ; zeros(1,length(v_ctrl)) , ones(1,length(v_case))];
[buc_y,buc_i] = sort(buc(1,:),'ascend');
buc = buc(:,buc_i);
vv = buc(1,:); xy = buc(2,:);
vgvm = [0 (vv(2:end)>vv(1:end-1))]; vvp = zeros(1,length(vv)); vvp(1)=1; ptmp=vvp(1); for (nl=2:length(vvp)); if (vgvm(nl)); vvp(nl)=nl; ptmp=nl; else; vvp(nl)=ptmp; end; end; vvp = vvp-1; vvp(find(vvp<1))=1;
vlvp = [(vv(1:end-1)<vv(2:end)) 0]; vvq = zeros(1,length(vv)); vvq(end)=length(vv); qtmp=vvq(end); for (nl=length(vvq)-1:-1:1); if (vlvp(nl)); vvq(nl)=nl; qtmp=nl; else; vvq(nl)=qtmp; end;end; vvq = vvq+1; vvq(find(vvq>length(vvq)))=length(vvq)+1;
yy = cumsum(xy(end:-1:1)); yy = [yy(end:-1:1) 0];
xij = find(xy==0);
auc = sum( yy(vvq(xij)) + yy(vvp(xij)+1) ) / 2 / length(v_ctrl) / length(v_case);

if plot_flag;
disp(sprintf('plotting fpr and tpr...'));
v_ctrl = reshape(v_ctrl,1,length(v_ctrl));
v_case = reshape(v_case,1,length(v_case));
hbins=linspace(0,K,32);hv_ctrl=hist(v_ctrl,hbins);hv_case=hist(v_case,hbins);
[D,I1,I2] = union(v_ctrl,v_case);
tpr = sum(repmat(v_case,length(D),1)>=repmat(transpose(D),1,length(v_case)),2)/length(v_case);
fpr = sum(repmat(v_ctrl,length(D),1)>=repmat(transpose(D),1,length(v_ctrl)),2)/length(v_ctrl);
figpos();
subplot(1,2,1);
hold on;
b=bar(hbins-0.25,hv_ctrl);set(b,'FaceColor','r');
b=bar(hbins+0.25,hv_case);set(b,'FaceColor','b');
hold off;
xlabel('variable value'); ylabel('# in bin'); title('histogram of v_ctrl (red) vs v_case (blue)');
subplot(1,2,2);
hold on;
plot(fpr,tpr,'-');
plot([0 1],[0,1],':');
hold off;
axis square; xlabel('false positive rate'); ylabel('true positive rate'); title('ROC-curve');
Auc = trapz(fpr(end:-1:1),tpr(end:-1:1)); disp(sprintf(' Auc (area under curve) = Auc = %0.4f, should be the same as the calculated auc = %0.4f',Auc,auc));
end;%if plot_flag;

end;%if (isempty(v_ctrl) | isempty(v_case));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fig1] = figpos()
set(0,'Units','pixels') ;
scnsize = get(0,'ScreenSize');
fig1 = figure;
position = get(fig1,'Position');
outerpos = get(fig1,'OuterPosition');
borders = outerpos - position;

edge = -borders(1)/2;
pos1 = [edge - 1800,  scnsize(4) * (2/3) - 500,  scnsize(3) - edge + 300,  scnsize(4) + 300];
set(fig1,'OuterPosition',pos1);
