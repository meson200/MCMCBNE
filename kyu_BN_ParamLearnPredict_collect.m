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
NoEnsembleSizes = numel(TestSizes);
% instance-by-instance prediction result
Probs.ProbByPat_trn = zeros(Npatients,NoEnsembleSizes+2);
Probs.ProbByPat_test = zeros(Npatients,NoEnsembleSizes+2);
Probs.ProbByPat_std_trn = zeros(Npatients,NoEnsembleSizes+2);
Probs.ProbByPat_std_test = zeros(Npatients,NoEnsembleSizes+2);
Mtrn = Probs.Prob_trn(:,3:NoEnsembleSizes+2);
Mtest = Probs.Prob_test(:,3:NoEnsembleSizes+2);
for i = 1:Npatients
    
    pt_row_trn = find(Probs.Prob_trn(:,1)==i);
    pt_row_test = find(Probs.Prob_test(:,1)==i);
    Pmat_trn = Mtrn(pt_row_trn,:);
    Pmat_test = Mtest(pt_row_test,:);
    Probs.ProbByPat_trn(i,1:2) = Probs.Prob_trn(pt_row_trn(1),1:2);
    Probs.ProbByPat_test(i,1:2) = Probs.Prob_test(pt_row_test(1),1:2);
    Probs.ProbByPat_std_trn(i,1:2) = Probs.Prob_trn(pt_row_trn(1),1:2);
    Probs.ProbByPat_std_test(i,1:2) = Probs.Prob_test(pt_row_test(1),1:2);
    
    Probs.ProbByPat_trn(i,3:end) = mean(Pmat_trn,1);
    Probs.ProbByPat_std_trn(i,3:end) = std(Pmat_trn,1,1)/sqrt(Nfolds)*2;
    Probs.ProbByPat_test(i,3:end) = mean(Pmat_test,1);
    Probs.ProbByPat_std_test(i,3:end) = std(Pmat_test,1,1)/sqrt(Nfolds)*2;

end
Probs.polarity = polarity;
Probs.EnsembleSizes = TestSizes;








