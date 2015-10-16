function Probs = kyu_LogisticModel(data)

% build a multivariate logistic prediction model
% using the variable set saved in a dataset 'data'
% the model has a L2 regularization term which coefficient will be tuned 
% in cross-validation fashion within each training dataset
%
% input:
% data: the bootstrapped/CVed dataset (see kyu_BN_GeneratePartition.m)
%
% output:
% Probs: class probability for every training/testing instances in every
% BS/CV folds
% .Prob_trn: P(class=2) for the patients (IDs in the 1st column) encountered during training,
% contains duplicates. 2nd column is the true value of a class
% .Prob_test: P(class=2) for the patients in the bootstrap OOB samples
% .ProbByPat_(std)_trn: Average/95%CI of training P(class) for each patient
% .ProbByPat_(std)_test: Average/95%CI of testing P(class) for each patient


% all the logistic model hyperparameters saved here
args_logistic = struct(); 
args_logistic.CVrep = 100; % number of repetition that shuffles CV partitioning
args_logistic.CVsize = 10; % size of cross validation set
args_logistic.L = 2; % order of a regularization term
args_logistic.lambda = [0.2 0.5 1 2 3 4 6]; % choices of L2 reg. coefficient (lambda) to try
args_logistic.alpha = 0.5; % elastic regression params. used only when L=1
args_logistic.metric = 'MSE'; % which metric to use for choosing the best L2 coefficient
args_logistic.cutoff = 0.5; % decision making boundary on logistic probability
args_logistic.maxiter = 500; % maximum iteration for logistic regression

data_in = data.KM; 
labels = data_in.Labels';
data_c = data_in.data_c_orig;
pts_test = data_in.patients_test;
pts_train = data_in.patients_train;
num_nodes = numel(labels);
trn_cases = size(data_c,2);
RPratio = numel(find(data_c(num_nodes,:)==2))/trn_cases;
Nrand = numel(pts_train);
lambda_opt = zeros(Nrand,1);
Prob_trn = [];
Prob_test = [];

for i = 1:Nrand
    % training
    % train LR (L2 regularization)
    data_trn_c = data_c(:,pts_train{i});
    outcome_trn = data_trn_c(num_nodes,:)-1;
    data_test_c = data_c(:,pts_test{i});
    
    [w_opt,lambda_opt(i),~] = kyu_LR_CV(data_trn_c(1:num_nodes-1,:)',outcome_trn',pts_train{i},args_logistic);
    [perf_trn,ROCXY_trn,Pth_trn,~,plist_trn] = kyu_LR_perftest(data_trn_c,num_nodes-1,w_opt,pts_train{i},RPratio);
    % test
    [perf_test,ROCXY_test,Pth_test,~,plist_test] = kyu_LR_perftest(data_test_c,num_nodes-1,w_opt,pts_test{i},RPratio);
    
    Prob_trn = [Prob_trn;plist_trn];
    Prob_test = [Prob_test;plist_test];
    
end
Probs.Prob_trn = Prob_trn;
Probs.Prob_test = Prob_test;
Probs.EnsembleSizes = 1; % we don't use an ensemble
Probs = AverageProbability_PatbyPat(Probs,Nrand);


