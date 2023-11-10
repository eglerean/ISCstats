% a script to test iscstats_ttest2_np.m
close all
clear all
rng(0)
% let's generate some connectivity matrices for our subjects
T=1000;
NG1=20;
NG2=20;
R = 246; % number of nodes
mkdir('networks')
mkdir('networks/group1')
mkdir('networks/group2')
for s=1:NG1
    adj=corr(randn(T,R));
    save(['networks/group1/net_' num2str(s) '.mat'],'adj');
end

for s=1:NG2
    adj=corr(randn(T,R));
    save(['networks/group2/net_' num2str(s) '.mat'],'adj');
end


% let's pick one node and test if group one is more similar in the pattern
% of links coming out of the node, versus group two

r=100;
% load the pattern of links coming out of that node for each subject

data_G1=[];
for s=1:NG1
    load(['networks/group1/net_' num2str(s) '.mat']); % we have a variable called adj
    vec = adj(:,r);
    % let's remove the r_th row because it will have all ones
    vec(r)=[];
    data_G1=[data_G1 vec];
end

data_G2=[];
for s=1:NG2
    load(['networks/group2/net_' num2str(s) '.mat']); % we have a variable called adj
    vec = adj(:,r);
    % let's remove the r_th row because it will have all ones
    vec(r)=[];
    data_G2=[data_G2 vec];
end



iscdata=corr([data_G1 data_G2]); % isc matrix with within group and across groups FC_ISC for the chosen node
out=iscstats_ttest2_np(iscdata,[ones(1,NG1) 2*ones(1,NG2)],5000);
disp(['   T-value = ' num2str(out.tval,2) ', p-value = ' num2str(out.pval,2)])

figure(1)
imagesc(iscdata,.1*[-1 1]);
axis square
xlabel('Subjects')
ylabel('Subjects')
title(['   T-value = ' num2str(out.tval,2) ', p-value = ' num2str(out.pval,2)])
colorbar

