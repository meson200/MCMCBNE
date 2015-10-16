function Probs = AverageProbability_PatbyPat(Probs,Nfolds)

% average out the predicted probabilities in an object Probs 
% to obtain bootstrap average class probability for each instance
%
% input: the Probs object that saves class probability for every training/testing instances in every
% BS/CV folds
% .Prob_trn: P(class=2) for the patients (IDs in the 1st column) encountered during training,
% contains duplicates. 2nd column is the true value of a class
% .Prob_test: P(class=2) for the patients in the bootstrap OOB samples
% Nfolds: number of bootstrap/CV folds
%
% output: the Probs object updated with the following additional fields
% .ProbByPat_(std)_trn: Average/95%CI of training P(class) for each patient
% .ProbByPat_(std)_test: Average/95%CI of testing P(class) for each patient

NoEnsembleSizes = numel(Probs.EnsembleSizes);
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