function stats = iscstats_ttest2_np(data,design,niter)
% ISCSTATS_TTEST2_NP - A "non parametric" two-sample T-test for comparing 
% inter-subject correlation valuesthat. Instead of
% relying on the t-distribtion, uses permutations of group labels to
% estimate the null distribution. 
%
%   - Usage:
%   stats = iscstats_ttest2_np(data,design,niter)
%   - Input:
%   data = an ISC matrix where each column and each row is a subject
%   design = a row vector containing the numbers 1 and 2 for the two groups to compares
%			
%   niter = number of permutations (recommended at least 5000)
%
%   - Output:
%   stats = a struct with the subfields
%       pvals = p-values for each datapoint; it returns in order the p-values
%       for the right tail and for the left tail
%       tvals = T-values for datapoint, positive tvals mean group 1 > group 2
%
% Notes: the null distribution is estimated using the matlab function
% ksdensity by interpolating the permuted data. The distribution is
% estimated over 200 points if niter<=5000, otherwise it is estimated over
% round(200*niter/5000) points, for greater precision.
%
% The script implements the approach explained by Gang Chen et al. in 10.1016/j.neuroimage.2016.05.023 and 10.1016/j.neuroimage.2016.08.029

% EG 2018-02-10 - enrico.glerean@aalto.fi

%% let's do some input validation first
Nsubj = size(data,1);   % total number of subjects
if(size(data,1) ~= size(data,2)) 
	error('The input data matrix should be squared')
end
if(sum(diag(data))~=Nsubj) % diagonals should have all ones as it is an ISC matrix
	error('The input data matrix should be an ISC matrix with ones on the main diagonal')
end

if(size(design,2) ~= Nsubj)
    error('Mismatched number of subjects: the number of columns of data variable  should match the number of columns of the design variable')
end
if(size(design,1) ~= 1)
    error('The design variable should only contain 1 row')
end
% let's get the group IDs
g1 = find(design==1);
g2 = find(design==2);
if((length(g1)+length(g2))~=Nsubj)
    error('The design variable should only contain numbers 1 and 2')
end
% a check on niter
if(niter<0)
    disp('The variable niter should be a positive integer, function will continue assuming niter=5000')
    niter=5000;
end
if(niter==0)
    disp('I will not do permutations');
end
stats.tval=isc_tt_np(data,g1,g2);
% computing pvalues
if(niter>0)
    pval = isc_tt_np_pval(data,g1,g2,niter,stats.tval);
end

stats.pval=pval;
end

function [tval iscmodel]=isc_tt_np(iscdata,g1,g2)
    % helper function similar to matlab function ttest2.m for the case of
    % groups with difference variance
	N = size(iscdata,1);
	iscmodel=zeros(N);
	iscmodel(g1,g1)=ones(length(g1));
	iscmodel(g2,g2)=2*ones(length(g2));
	iscmodel=triu(iscmodel,1); % only top triangle
	ids1=find(iscmodel==1);
	ids2=find(iscmodel==2);
    xnans = isnan(iscdata(ids1));
    if any(xnans(:))
        nx = sum(~xnans,2);
    else
        nx = length(ids1); 
    end
    ynans = isnan(iscdata(ids2));
    if any(ynans(:))
        ny = sum(~ynans,2);
    else
        ny = length(ids2); % a scalar, => a scalar call to tinv
    end

    difference = nanmean(iscdata(ids1)) - nanmean(iscdata(ids2));
    
    s2x = nanvar(iscdata(ids1));
    s2y = nanvar(iscdata(ids2));
    s2xbar = s2x ./ nx;
    s2ybar = s2y ./ ny;
    se = sqrt(s2xbar + s2ybar);
    if(any(se == 0) || any(isnan(se)))
        se(find(se==0))=Inf;
		se(find(isnan(se)))=Inf;
		%error('Group variance seems to be null or NaN, please check your data')
    end
    tval = difference ./ se;

end

function pval = isc_tt_np_pval(iscdata,g1,g2,niter,tval)
    outiter=zeros(niter,1);
    ND=length(g1)+length(g2);
    for iter=1:niter
        perm=randperm(ND);
        % one could add a test to see that they are indeed permuted
        temp=iscdata(perm,perm);
        outiter(iter)=isc_tt_np(temp,g1,g2);
    end
    NCDF=200;
    if(niter>5000)
        NCDF=round(200*niter/5000);
    end
    [fi xi]=ksdensity(outiter,'function','cdf','npoints',NCDF);
    
    % trick to avoid NaNs, we approximate the domain of the CDF between
    % -Inf and Inf using the atanh function and the eps matlab precision
    % variable
    
    pval_left=interp1([atanh(-1+eps) xi atanh(1-eps)],[0 fi 1],tval); 
    pval_right=1-pval_left;
    pval=min([pval_right pval_left]);
end


