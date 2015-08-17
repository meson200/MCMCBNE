function [Probs,Mreturn] = kyu_BN_BMA(gs_top,posterior,data_trn,data_test,Nrand,PtsTrn,PtsTest,ClassRate,EnsembleSizes,alpha_d)

% Apply the BN ensemble to bootstrap replicates to obtain class probability
% in bayesian model averaging (BMA) fashion
%
% list of input arguments:
% gs_top: a cell array of BN DAG ensemble, size N_f * N_e 
% (N_f:number of validation folds, N_e: max. number of DAGs sampled
% from each fold)
% posterior: posterior prob. of each graph in gs_top
% data_trn: datasets for training parameters (the same data used for graph learning) 
% data_test: datasets for prediction testing
% Nrand: number of folds (numel(data_trn))
% PtsTrn: instance identifier for training data folds (number 1~Number_of_instances), 
% (not the same as patient IDs!)
% PtsTest: instance identifier for testing data folds
% ClassRate: prior probability P(Class) to be used as a default operating pt a classifier
% EnsembleSizes: a list of ensemble sizes at which performance is evaluated
%
% the following are returned as output:
% Probs: class probability 
% Prob{i}.Prob_trn: P(class) for training instances in a fold i
% Prob{i}.Prob_test: P(class) for testing instances in a fold i
% Mreturn: polarity of dependency between nodes (Nnodes*Nnodes)
% Mreturn(i,j)>0: P(j=2|i=2) > P(j=2|i=1) (positive correlation)

ns = 2*ones(1,size(data_trn{1},1));
nvar = size(data_trn{1},1);
MBig = cell(Nrand,1);
occurBig = cell(Nrand,1);
Prob_trn = cell(Nrand,1);
Prob_test = cell(Nrand,1);
NoTopGraphs = max(EnsembleSizes);

parfor i = 1:Nrand
   trn = data_trn{i};
   test = data_test{i};
   post = posterior(i,:);
   bnet1 = cell(NoTopGraphs,1);
   eg1 = cell(NoTopGraphs,1);
   gs = gs_top(i,:);
   M = zeros(nvar,nvar);
   dagtemp = zeros(nvar,nvar);
   for j = 1:NoTopGraphs
        [bnet1{j},eg1{j}] = kyu_BN_paramlearn(trn,gs{j},ns,alpha_d,2);
        dagtemp = dagtemp + bnet1{j}.dag; 
        M = M + dag_polarity(bnet1{j},eg1{j});
   end 
   M = M./dagtemp;
   M(isnan(M)) = 0;
   MBig{i} = M;
   occur = zeros(nvar,nvar);
   for p = 1:nvar
       for q = 1:nvar
           if M(p,q) ~=0 occur(p,q) = 1; end
       end
   end
   occurBig{i} = occur;
   Ps_test = kyu_BN_ClassPs_Bayesian(test,bnet1,eg1,NoTopGraphs);
   Ps_trn = kyu_BN_ClassPs_Bayesian(trn,bnet1,eg1,NoTopGraphs);
   [Prob_trn{i},Prob_test{i}] = PerfTest_group(test,trn,Ps_test,Ps_trn,PtsTest{i},PtsTrn{i},post,ClassRate,EnsembleSizes);
   
end

Probs = struct('Prob_trn',Prob_trn,'Prob_test',Prob_test,'EnsembleSizes',EnsembleSizes);
MSum = 0; 
occurSum = 0;
for i = 1:Nrand
    MSum = MSum + MBig{i};
    occurSum = occurSum + occurBig{i};
end
Mreturn = MSum./occurSum;
Mreturn(isnan(Mreturn)) = 0;


% ---------------------------subfunctions----------------------------------


function [Prob_trn,Prob_test] = PerfTest_group(test,trn,Ps_test,Ps_trn,ID_test,ID_trn,post,ClassRate,EnsembleSizes)

Prob_trn = [];
Prob_test = [];
for p = 1:numel(EnsembleSizes)
    [~,~,~,plist_p_test] = kyu_BN_perftest_Bayesian(test,Ps_test,post,EnsembleSizes(p),ID_test,ClassRate);
    [~,~,~,plist_p_train] = kyu_BN_perftest_Bayesian(trn,Ps_trn,post,EnsembleSizes(p),ID_trn,ClassRate);
    Prob_trn = [Prob_trn plist_p_train(:,3)];
    Prob_test = [Prob_test plist_p_test(:,3)];
end
Prob_trn = [plist_p_train(:,1:2) Prob_trn];
Prob_test = [plist_p_test(:,1:2) Prob_test];




