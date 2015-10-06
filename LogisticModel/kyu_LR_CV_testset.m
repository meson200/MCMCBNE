function [w_opt,lambda_opt,err,auc] = kyu_LR_CV_testset(X_train,Y_train,X_test,Y_test,L,lambda,varargin)

% regularized logistic regression
% inputs:
% L: order of regularization
% lambda: a list of different lambda values to test
% difference with kyu_LR_CV.m: it takes a dedicated test set to tune a
% lambda (a parameter for a regularization term)

args = varargin;
nargs = length(args);
metric = 'MSE'; % MSE is the default evaluation metric
for i=1:2:nargs
 switch args{i},
  case 'metric',  metric = args{i+1};
  case 'verbose', verbose = args{i+1};
  end
end

% initialize parameters
npar = size(X_train,2);
maxiter = 500; % maximum iteration for logistic regression
alpha = 1;
err_trntrn = zeros(numel(lambda),1);
ROCcurves = cell(numel(lambda),1);

dim = zeros(numel(lambda),1);
auc = zeros(numel(lambda),1);
auc_trntrn = zeros(numel(lambda),1);
auc_test = zeros(numel(lambda),1);
err = zeros(numel(lambda),1);
err_trntrn = zeros(numel(lambda),1);
err_test = zeros(numel(lambda),1);

for i = 1:numel(lambda)

    perf_temp = zeros(7,1);
    X_trntrn = X_train; 
    Y_trntrn = Y_train;
    Y_test_temp = Y_test+1;
    Y_train_temp = Y_train+1;
    data_val = [X_test'; Y_test_temp'];
    data_trn = [X_train'; Y_train_temp'];

    alpha = 2-L; % determines ridge(L=2,a=0) or LASSO(L=1,a=1)   

    [w,FitInfo] = lassoglm(X_trntrn,Y_trntrn,'binomial','Lambda',lambda(i),'Alpha',alpha); 
    intercept = FitInfo.Intercept;
    dim(i) = numel(find(w));
    h_test = h_log(w,intercept,X_test);
    h_trntrn = h_log(w,intercept,X_trntrn);
    % training error
    err_trntrn(i) = (h_trntrn - Y_trntrn)'*(h_trntrn - Y_trntrn)/(2*numel(Y_trntrn));
    ID = zeros(numel(Y_trntrn),1);
    RPr = zeros(numel(Y_trntrn),1);
    [perf_temp_trn,~,~,~,~] = kyu_LR_perftest(data_trn,npar,[w; intercept],ID,RPr);
    auc_trntrn(i) = perf_temp_trn(6);
    % testing error
    err_test(i) = (h_test - Y_test)'*(h_test - Y_test)/(2*numel(Y_test)); 
    ID = zeros(numel(Y_test),1);
    RPr = zeros(numel(Y_test),1);
    [perf_temp_test,~,~,~] = kyu_LR_perftest(data_val,npar,[w; intercept],ID,RPr);
    auc_test(i) = perf_temp_test(6);
    % estimate .632 bootstrap error  
    err(i) = 0.368*err_trntrn(i) + 0.632*err_test(i); 
    auc(i) = 0.368*auc_trntrn(i) + 0.632*auc_test(i);
   
end


% optimize the lambda
if strcmp(metric,'MSE')
    if verbose==1
        disp(['tested lambda: ',num2str(lambda)]);
        disp(['Cross-validation MSE: ',num2str(err')]);
    end
    [~,bestindx] = min(err); 
else
    if verbose==1
        disp(['tested lambda: ',num2str(lambda)]);
        disp(['Cross-validation AUC: ',num2str(auc')]);
    end
    [~,bestindx] = max(auc); 
end

% refit the model with the optimized parameter to get the model coefficients
[w_opt,FitInfo] = lassoglm(X_train,Y_train,'binomial','Lambda',lambda(bestindx),'Alpha',alpha);
w_opt = [w_opt; FitInfo.Intercept];
lambda_opt = lambda(bestindx);

% -------------------------subfunctions------------------------------------

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
