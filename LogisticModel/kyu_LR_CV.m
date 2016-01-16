function [w_opt,lambda_opt,bestperf] = kyu_LR_CV(X_train,Y_train,ID,args_logistic)

% regularized logistic regression
% inputs:
% L: order of regularization
% lambda: a list of different lambda values to test

% optional arguments:
% 'metric': which metric to use for choosing the best L2 coefficient
%           'MSE' (mean square error) or 'AUC' (area under the ROC curve)
% 'cutoff': decision making boundary on logistic probability
%           y=2 when P(X)>cutoff, 1 otherwise
%           default = 0.5
%           doesn't affect AUC (b/c all the cutoff values are averaged)
% outputs:
% w_opt: logistic model coefficients (the last element is an intercept)
% lambda_opt: the optimized lambda value 
% bestperf: 

CVrep = args_logistic.CVrep; % number of repetition that shuffles CV partitioning
CVsize = args_logistic.CVsize; % size of cross validation set
L = args_logistic.L; % order of a regularization term
lambda = args_logistic.lambda; % choices of lambda to try
alpha = args_logistic.alpha; % elastic regression params. used only when L=1
metric = args_logistic.metric; % which metric to use for choosing the best L2 coefficient
cutoff = args_logistic.cutoff;
maxiter = args_logistic.maxiter;

% initialize parameters
l = size(X_train,1);
npar = size(X_train,2);
err_trntrn = zeros(numel(lambda),CVrep);
err_trnval = zeros(numel(lambda),CVrep);
auc_trnval = zeros(numel(lambda),CVrep);
ROCcurves = cell(numel(lambda),1);
lst = zeros(numel(lambda),4);

dim = zeros(numel(lambda),CVrep);
w_temp = zeros(numel(lambda),npar);
perf = zeros(numel(lambda),7);
perf_stdev = zeros(numel(lambda),7);

% this tunes a coefficient for regularization 
% skipped if no regul. order is provided
if L ~= 0 

    for i = 1:numel(lambda)

        p_RP = [];
        perf_temp = zeros(CVrep,7);
        for k = 1:CVrep
        % N-fold cross validation for each lambda
            [trnindx,valindx] = randdraw(numel(Y_train),CVsize);
            X_trntrn = X_train(trnindx,:); 
            Y_trntrn = Y_train(trnindx,:);
            X_trnval = X_train(valindx,:);
            Y_trnval = Y_train(valindx,:);
            Y_trnval_temp = Y_trnval+1;
            data_val = [X_trnval'; Y_trnval_temp'];

            if L == 2 % L2 regularization
                [w,conv] = kyu_LR_L2opt(X_trntrn,Y_trntrn,lambda(i),maxiter);
                h_trntrn = h_log(w(1:end-1),w(end),X_trntrn);
                h_trnval = h_log(w(1:end-1),w(end),X_trnval);
                [perf_temp(k,:),~,~,p_RP_temp,~] = kyu_LR_perftest(data_val,npar,w,ID,cutoff);
                p_RP = [p_RP p_RP_temp];
                err_trntrn(i,k) = (h_trntrn - Y_trntrn)'*(h_trntrn - Y_trntrn)/(2*(l-CVsize));
                err_trnval(i,k) = (h_trnval - Y_trnval)'*(h_trnval - Y_trnval)/(2*CVsize); 
                auc_trnval(i,k) = perf_temp(k,6);

            elseif L == 1 % L1 (lasso)
               [w,FitInfo] = lassoglm(X_trntrn,Y_trntrn,'binomial','Lambda',lambda(i),'Alpha',alpha); 
               intercept = FitInfo.Intercept;
               err = FitInfo.Deviance;
               dim(i,k) = numel(find(w));
               %err_trntrn(i,k) = FitInfo.MSE;
               h_trntrn = h_log(w,intercept,X_trntrn);
               h_trnval = h_log(w,intercept,X_trnval);
               [perf_temp(k,:),~,~,p_RP_temp,~] = kyu_LR_perftest(data_val,npar,[w; intercept],ID,cutoff);
               p_RP = [p_RP p_RP_temp];
               err_trntrn(i,k) = (h_trntrn - Y_trntrn)'*(h_trntrn - Y_trntrn)/(2*numel(Y_trntrn));
               err_trnval(i,k) = (h_trnval - Y_trnval)'*(h_trnval - Y_trnval)/(2*CVsize);
               auc_trnval(i,k) = perf_temp(k,6);
            end


        end

        ROCcurves{i} = p_RP;

        perf(i,:) = nanmean(perf_temp,1);
        perf_stdev(i,:) = nanstd(perf_temp,1)/sqrt(CVsize);
        dim_n = mean(dim,2);


    end

    err_trntrn_n = nanmean(err_trntrn,2);
    err_trnval_n = nanmean(err_trnval,2);
    err_trnval_s = 2*nanstd(err_trnval,1,2)/sqrt(CVrep);
    auc_trnval_n = nanmean(auc_trnval,2);
    auc_trnval_s = 2*nanstd(auc_trnval,1,2)/sqrt(CVrep);


    disp(['tested lambda: ',num2str(lambda)]);

    if strcmp(metric,'MSE') 
        disp(['Cross-validation MSE: ',num2str(err_trnval_n')]);
        disp(['95% confidence inverval: ',num2str(err_trnval_s')]);
        [bestperf,bestindx] = min(err_trnval_n);
    else
        disp(['Cross-validation AUC: ',num2str(auc_trnval_n')]);
        disp(['95% confidence inverval: ',num2str(auc_trnval_s')]);
        [bestperf,bestindx] = max(auc_trnval_n); 
    end
    lambda_opt = lambda(bestindx);

end


if L==2
    [w_opt,conv] = kyu_LR_L2opt(X_train,Y_train,lambda_opt,maxiter);
elseif L==1
    [w_opt,FitInfo] = lasso(X_train,Y_train,'Lambda',lambda_opt,'Alpha',alpha);
    w_opt = [w_opt; FitInfo.Intercept];
else    
    [mu,w_opt,se, LR_chi2,convState] = drxlr_logistic_regression(X_train,Y_train,maxiter,1e-6);
    lambda_opt = 0;
    bestperf = 0;
end

% performance of the best model (best lambda)


    

function hval = h_log(w,intercept,X)
% logit link function
Xplus = [ones(size(X,1),1) X];
expo = Xplus * [intercept; w];
hval = 1./(1+exp(-expo));


function Ldel = L_del(err,X)
% computes gradient of log likelihood
l = size(X,1);
d = size(X,2) ;   
Ldel = zeros(1,d);    
for k = 1:l
    Ldel = Ldel + err(k)*X(k,:);
end

function [trnindx,valindx] = randdraw(N,k)
    % divides the list 1:N into 2 sets with length k and N-k
    valindx = randsample(1:N,k);
    trnindx = setdiff(1:N,valindx);
