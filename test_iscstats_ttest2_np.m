% a script to test iscstats_ttest2_np.m
close all
clear all
rng(0)
% let's generate some time series
T=1000;
NG1=20;
NG2=20;
G1=randn(T,NG1); % data for group one
G2=randn(T,NG2); % data for group two

% case 1, there is nothing different between the two groups
iscdata=corr([G1 G2]); % isc matrix with within group and across groups ISC
out=iscstats_ttest2_np(iscdata,[ones(1,NG1) 2*ones(1,NG2)],5000);
disp('When the two groups are similar:')
disp(['   T-value = ' num2str(out.tval,2) ', p-value = ' num2str(out.pval,2)])

figure(1)
subplot(1,2,1)
imagesc(iscdata,.1*[-1 1]);
axis square
xlabel('Subjects')
ylabel('Subjects')
title(['   T-value = ' num2str(out.tval,2) ', p-value = ' num2str(out.pval,2)])
colorbar

% case 2; group one is more synchronised

G1=G1+0.2*repmat(randn(T,1),1,NG1); % all subjects in group 1 have a common signal added
iscdata=corr([G1 G2]); % isc matrix with within group and across groups ISC
out=iscstats_ttest2_np(iscdata,[ones(1,NG1) 2*ones(1,NG2)],5000);
disp('When the two groups are different:')
disp(['   T-value = ' num2str(out.tval,2) ', p-value = ' num2str(out.pval,2)])
subplot(1,2,2)
imagesc(iscdata,.1*[-1 1]);
axis square
xlabel('Subjects')
ylabel('Subjects')
title(['   T-value = ' num2str(out.tval,2) ', p-value = ' num2str(out.pval,2)])
colorbar
