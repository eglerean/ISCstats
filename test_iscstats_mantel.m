% a script to test iscstats_mantel.m
close all
clear all
rng(0)
% let's generate some time series for NG subjects
T=1000;
N=20;


%% case 1 no similarity structure between subjects 
data=randn(T,N); % time series data
% let's create a mantel matrix from some subjects NB behavioural scores
NB = 10;
behav_data=round(20*rand(NB,N));
iscdata=corr(data);
behav_sim=corr(behav_data);
[r p]=iscstats_mantel(iscdata,behav_sim,5000,'spearman')

%% case 2, some subjects contain the same signal and the same subjects are more similar behaviourly
data=data+1/(N)*repmat([1:N],T,1).*repmat(randn(T,1),1,N);
% let's create a mantel matrix from some subjects NB behavioural scores
NB = 10;
behav_data=round(20*(rand(NB,N)+1/(N)*repmat([1:N],NB,1).*repmat(rand(NB,1),1,N)));
iscdata=corr(data);
behav_sim=corr(behav_data);
[r p]=iscstats_mantel(iscdata,behav_sim,5000,'spearman')


