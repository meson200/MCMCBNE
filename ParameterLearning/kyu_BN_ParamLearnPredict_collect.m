function [Probs,polarity] = kyu_BN_ParamLearnPredict_collect(data,ClusterOutput)

% combine batched jobs for classification test into a single struct  
% input: 
% data: the bootstrapped/CVed dataset used for training/testing
% ClusterOutput: output of kyu_BN_ParamLearnPredict_submit
%
% output:
% Probs: class probability for every training/testing instances in every
% BS/CV folds
% .Prob_trn: P(class=2) for the patients (IDs in the 1st column) encountered during training,
% contains duplicates. 2nd column is the true value of a class
% .Prob_test: P(class=2) for the patients in the bootstrap OOB samples
% .ProbByPat_(std)_trn: Average/95%CI of training P(class) for each patient
% .ProbByPat_(std)_test: Average/95%CI of testing P(class) for each patient
% polarity: direction of correlation between nodes (Nnodes*Nnodes) which was
% obtained from the trained BN parameters (CPTs), as a function of ensemble
% size
% polarity(i,j)>0: P(j=2|i=2) > P(j=2|i=1) (positive correlation)

Nfolds = 0;
NJobs = numel(ClusterOutput);
Npatients = size(data.KM.data_orig_missing,2);
ProbTrnTmp = [];
ProbTestTmp = [];
TestSizes = ClusterOutput{1}{1}(1).EnsembleSizes; 
polarity = repmat({0},[1 numel(TestSizes)]);

for i = 1:NJobs
    result_temp = ClusterOutput{i};
    Pb = result_temp{1};
    TestSizes = result_temp{1}(1).EnsembleSizes;
    for p = 1:numel(TestSizes)
        polarity{p} =  polarity{p} + result_temp{2}{p};
    end
    Nfolds = Nfolds + numel(Pb);
    for j = 1:numel(Pb)
        ProbTrnTmp = cat(1,ProbTrnTmp,Pb(j).Prob_trn);
        ProbTestTmp = cat(1,ProbTestTmp,Pb(j).Prob_test);
    end
end
for p = 1:numel(TestSizes)
    polarity{p} =  polarity{p}/NJobs;
end

Probs.Prob_trn = ProbTrnTmp;
Probs.Prob_test = ProbTestTmp;
Probs.EnsembleSizes = TestSizes;

% average out the predicted probabilities to obtain 
% bootstrap average class probability for each instance
Probs = AverageProbability_PatbyPat(Probs,Nfolds);
Probs.polarity = polarity;









