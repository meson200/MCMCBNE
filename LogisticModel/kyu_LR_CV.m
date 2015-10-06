function [w_opt,lambda_opt,bestperf] = kyu_LR_CV(X_train,Y_train,L,CVrep,CVsize,ID,varargin)

% regularized logistic regression
% inputs:
% L: order of regularization
% lambda: a list of different lambda values to test
% outputs:
% w_opt: logistic model coefficients (the last element is an intercept)
% lambda_opt: the optimized lambda value 
% bestperf: 

args = varargin;
nargs = length(args);
metric = 'MSE'; % MSE is the default evaluation metric
for i=1:2:nargs
 switch args{i},
  case 'metric'
      metric = args{i+1};
  case 'cutoff'
      cutoff = args{i+1};
  end
end

l = size(X_train,1);
%X_train = [X_train ones(l,1)];
% initialize parameters
npar = size(X_train,2);
maxiter = 500; % maximum iteration for logistic regression

% range of values for RP model
%lambda = logspace(-4,3,10);

%lambda = [0.005 0.01 0.02 0.04 0.1 0.2 0.4];
lambda = [0.2 0.5 1 2 3 4 6];
alpha = 0.5;

% range of values for texture
%lambda = [0.001 0.002 0.004 0.008 0.016];
%lambda = [0.001];

err_trntrn = zeros(numel(lambda),CVrep);
ROCcurves = cell(numel(lambda),1);
lst = zeros(numel(lambda),4);

dim = zeros(numel(lambda),CVrep);
w_temp = zeros(numel(lambda),npar);
perf = zeros(numel(lambda),7);
perf_stdev = zeros(numel(lambda),7);


for i = 1:numel(lambda)
    

    acc = 0;
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
            
            % added on July 23, for the combined dataset only
            [perf_temp(k,:),~,~,p_RP_temp,~] = kyu_LR_perftest(data_val,npar,w,ID,cutoff);
%            p_RP_temp = [];
%            lsttmp = [0 0 0 0];
%             for j = 1:CVsize
%                [acc_inc,lst_inc,pRPseg] = kyu_LR_classify(data(:,j),w); 
%                %if acc_inc ==0 wrongpts = [wrongpts j]; end
%                p_RP_temp = [p_RP_temp pRPseg];
%                acc = acc + acc_inc;   
%                lsttmp = lsttmp + lst_inc;
%             end
%            lst(i,:) = lst(i,:) + lsttmp;
            p_RP = [p_RP p_RP_temp];
            err_trntrn(i,k) = (h_trntrn - Y_trntrn)'*(h_trntrn - Y_trntrn)/(2*(l-CVsize));
            err_trnval(i,k) = (h_trnval - Y_trnval)'*(h_trnval - Y_trnval)/(2*CVsize); 
            auc_trnval(i,k) = perf_temp(k,6);
            
        else % L1 (lasso)
           [w,FitInfo] = lassoglm(X_trntrn,Y_trntrn,'binomial','Lambda',lambda(i),'Alpha',alpha); 
           intercept = FitInfo.Intercept;
           err = FitInfo.Deviance;
           dim(i,k) = numel(find(w));
           %err_trntrn(i,k) = FitInfo.MSE;
           h_trntrn = h_log(w,intercept,X_trntrn);
           err_trntrn(i,k) = (h_trntrn - Y_trntrn)'*(h_trntrn - Y_trntrn)/(2*numel(Y_trntrn));
           h_trnval = h_log(w,intercept,X_trnval);
           p_RP = [];
           err_trnval(i,k) = (h_trnval - Y_trnval)'*(h_trnval - Y_trnval)/(2*CVsize); 
           [perf_temp,~,~,p_RP_temp,~] = kyu_LR_perftest(data_val,npar,w,cutoff);
           auc_trnval(i,k) = perf_temp(6);
        end  
    end
    
    ROCcurves{i} = p_RP;
    
    perf(i,:) = nanmean(perf_temp,1);
    perf_stdev(i,:) = nanstd(perf_temp,1)/sqrt(CVsize);
    dim_n = mean(dim,2);
    
    
end


err_trntrn_n = mean(err_trntrn,2);
err_trnval_n = mean(err_trnval,2);
err_trnval_s = 2*nanstd(err_trnval,1,2)/sqrt(CVrep);
auc_trnval_n = mean(auc_trnval,2);
auc_trnval_s = 2*nanstd(auc_trnval,1,2)/sqrt(CVrep);


disp(['tested lambda: ',num2str(lambda)]);

if strcmp(metric,'MSE') 
    disp(['Cross-validation MSE: ',num2str(err_trnval_n')]);
    disp(['95% confidence inverval: ',num2str(err_trnval_s')]);
    [~,minindx] = min(err_trnval_n); 
else
    disp(['Cross-validation AUC: ',num2str(auc_trnval_n')]);
    disp(['95% confidence inverval: ',num2str(auc_trnval_s')]);
    [~,minindx] = max(auc_trnval_n); 
end
if L==2
    [w_opt,conv] = kyu_LR_L2opt(X_train,Y_train,lambda(minindx),maxiter);
else
    [w_opt,FitInfo] = lasso(X_train,Y_train,'Lambda',lambda(minindx),'Alpha',alpha);
    w_opt = [w_opt; FitInfo.Intercept];
end

% performance of the best model (best lambda)
bestperf = perf(minindx,:);
bestperf_std = perf_stdev(minindx,:);
if L==2
    [ROC_X,ROC_Y,~,auc] = perfcurve(ROCcurves{minindx}(1,:),ROCcurves{minindx}(2,:),2);
    ROCXY = [ROC_X ROC_Y];
    perf(6) = auc;
end 
lambda_opt = lambda(minindx);

    

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
