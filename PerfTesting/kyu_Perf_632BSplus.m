function [Perfs_632p,Perfs_632pM] = kyu_Perf_632BSplus(Probs,Probs_fitting,data)

% obtain the bootstrap Probsormance of a model by combining bootstrap testing result (Probs) and
% fitting (Probs_fitting) according to the .632+ bootstrap method
% inputs:
% Probs: P(class) for every training/testing instances
% Probs_fitting: fitted values of P(class) when the entire data is used for
% building a model
% data: the bootstrapped/CVed dataset (see kyu_BN_GeneratePartition.m)
%
% outputs:
% Perfs_632p: a struct containing the ROC metrics listed below. 
% Each column of the metrics represents a result from a correspnding ensemble size in 
% Probs.EnsembleSizes.
% .auc: area under the curve. 
% .se_opt: sensitivity at an optimal operating point (max Youden's index)
% .sp_opt: specificity at an optimal operating point (max Youden's index)
% (for the above 3, 1st row: average, 2nd row: 95% confidence interaval)
% .ROCY: Y-axis points (TPR) of ROC curves. X-axis points are fixed at
% regular intervals (0:0.01:1)
% .ROCY_ste: 95% CE of ROC y values
% Perfs_632pM: stores the raw samples of the 5 performance metrics 


% curves (first row) and 95% confidence intervals (second row)
% columns are used for different ensemble sizes
% Perfs_632pM: .632+ adjusted Probsormance for each bootstrap replicate

Nmodel = numel(Probs.EnsembleSizes); % number of different ensemble sizes tested
data_in = data.MI;
Nrand = numel(data_in.train);
Perfs_632p = struct('auc',zeros(2,Nmodel),'se_opt',zeros(2,Nmodel),'sp_opt',zeros(2,Nmodel),'ROCY',zeros(101,1));
Perfs_632pM = struct('auc',zeros(Nrand,Nmodel),'se_opt',zeros(Nrand,Nmodel),'sp_opt',zeros(Nrand,Nmodel));
ROCY_test = cell(Nmodel,1);

for i = 1:Nmodel
    
    % get 0.632+ AUC
    [~,TPR_trn,t_trn,auc_trn] = perfcurve(Probs_fitting.Prob_trn(:,2),Probs_fitting.Prob_trn(:,i+2),2,'XVals',0:0.01:1,'UseNearest','off');
    row_start = 1;
    ROCY_test{i} = zeros(101,Nrand);
    t_test = zeros(101,Nrand);
    TPR_test = zeros(101,Nrand);
    auc_test = zeros(Nrand,1);
    for m = 1:Nrand
       setsize = numel(data_in.patients_test{m});
       row_end = row_start + setsize-1;
       Pmat = Probs.Prob_test(row_start:row_end,:);
       try
            [~,TPR_test(:,m),t_test(:,m),auc_test(m)] = perfcurve(Pmat(:,2),Pmat(:,i+2),2,'XVals',0:0.01:1,'UseNearest','off');
       catch % if a testing set is filled with one class
            TPR_test(:,m) = NaN;
            t_test(:,m) = NaN;
            auc_test(m) = NaN;
       end
       row_start = row_end+1;
    end
    [Perfs_632p.auc(1,i),Perfs_632p.auc(2,i),Perfs_632pM.auc(:,i)] = kyu_632plusbootstrap(auc_trn,auc_test,0.5);

    % get 0.632+ ROC curve
    TPR_list = zeros(101,Nrand);
    TPR_632p = zeros(101,1);
    TPR_632p_ste = zeros(101,1);
    for j = 1:101
        % calculate no-information error rate (q) (Adler & Lausen 2009)
        ss = numel(Probs_fitting.Prob_trn(:,i+2));
        q = sum(Probs_fitting.Prob_trn(:,i+2)>t_trn(j));
        q = q/ss;
        disp(q)
        [TPR_632p(j),TPR_632p_ste(j),TPR_list(j,:)] = kyu_632plusbootstrap(TPR_trn(j),TPR_test(j,:),q);
    end
    FPR = 0:0.01:1;
    FPR = FPR';
    se_opt_list = zeros(Nrand,1);
    sp_opt_list = zeros(Nrand,1);
    for m =1:Nrand
        youden = TPR_list(:,m) + ones(101,1) - FPR;    
        [~,maxindx] = max(youden);
        se_opt_list(m) = TPR_list(maxindx,m);   
        sp_opt_list(m) = 1-FPR(maxindx);
    end
    
    eq_sample_size = numel(find(~isnan(se_opt_list)));
    Perfs_632p.se_opt(1,i) = nanmean(se_opt_list);
    Perfs_632p.se_opt(2,i) = 2*nanstd(se_opt_list)/sqrt(eq_sample_size);
    Perfs_632p.sp_opt(1,i) = nanmean(sp_opt_list);
    Perfs_632p.sp_opt(2,i) = 2*nanstd(sp_opt_list)/sqrt(eq_sample_size);
    Perfs_632p.ROCY(:,i) = TPR_632p;
    Perfs_632p.ROCY_ste(:,i) = TPR_632p_ste;
    Perfs_632pM.se_opt(:,i) = se_opt_list';
    Perfs_632pM.sp_opt(:,i) = sp_opt_list';
    
end


