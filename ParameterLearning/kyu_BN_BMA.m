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

useEM = 2; % use EM by default
ns = 2*ones(1,size(data_trn{1},1));
nvar = size(data_trn{1},1);
Prob_trn = cell(Nrand,1);
Prob_test = cell(Nrand,1);
NoTopGraphs = max(EnsembleSizes);
Pos_big = cell(Nrand,numel(EnsembleSizes));
Neg_big = cell(Nrand,numel(EnsembleSizes));

parfor i = 1:Nrand
   trn = data_trn{i};
   test = data_test{i};
   post = posterior(i,:);
   bnet1 = cell(NoTopGraphs,1);
   eg1 = cell(NoTopGraphs,1);
   gs = gs_top(i,:);
   dagtemp = zeros(nvar,nvar);
   Pos_row = cell(1,numel(EnsembleSizes));
   Neg_row = cell(1,numel(EnsembleSizes));
   pos = 0;
   neg = 0;
   for j = 1:NoTopGraphs
        [bnet1{j},eg1{j}] = kyu_BN_paramlearn(trn,gs{j},ns,alpha_d,useEM);
        dagtemp = dagtemp + bnet1{j}.dag; 
        % pm: links with positive correlation, nm: negative
        [pm,nm] = dag_polarity(bnet1{j},eg1{j});  
        pos = pos + pm;
        neg = neg + nm;
        [ism,ind] = ismember(j,EnsembleSizes);
        if ism
            Pos_row{ind} = pos;
            Neg_row{ind} = neg;
        end
   end 
   Pos_big(i,:) = Pos_row;
   Neg_big(i,:) = Neg_row;
   Ps_test = kyu_BN_ClassPs_Bayesian(test,bnet1,eg1,NoTopGraphs);
   Ps_trn = kyu_BN_ClassPs_Bayesian(trn,bnet1,eg1,NoTopGraphs);
   [Prob_trn{i},Prob_test{i}] = PerfTest_group(test,trn,Ps_test,Ps_trn,PtsTest{i},PtsTrn{i},post,ClassRate,EnsembleSizes);
   
end

Probs = struct('Prob_trn',Prob_trn,'Prob_test',Prob_test,'EnsembleSizes',EnsembleSizes);
Mreturn = cell(1,numel(EnsembleSizes));
for j = 1:numel(EnsembleSizes)
    pos_sum = 0;
    neg_sum = 0;
    for i = 1:Nrand
        pos_sum = pos_sum + Pos_big{i,j};
        neg_sum = neg_sum + Neg_big{i,j};
    end
    MM = (pos_sum - neg_sum)./(pos_sum + neg_sum);
    MM(isnan(MM)) = 0;
    Mreturn{j} = MM;
end


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




