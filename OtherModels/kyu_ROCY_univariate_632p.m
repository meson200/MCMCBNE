function [auc_list,se_opt_list,sp_opt_list,ROCY_632p,ROCY_632p_ste] = kyu_ROCY_univariate_632p(data_orig,outcome_orig,pts_test)

% evaluate bootstrap 632 plus ROC metrics

% training TPR
pol = 1;
Nrand = length(pts_test);
[~,TPR_trn,t_trn,auc_trn] = perfcurve(outcome_orig,data_orig,2,'XVals',0:0.01:1,'UseNearest','off');
if auc_trn<0.5
    pol = -1; % if predicting negative, flip the predictor values
end
[~,TPR_trn,t_trn,auc_trn] = perfcurve(outcome_orig,data_orig*pol,2,'XVals',0:0.01:1,'UseNearest','off');

% testing TPR
auc_632p = zeros(2,1);
auc_test = zeros(Nrand,1);
TPR_test = zeros(101,Nrand);
TPR_632p = zeros(101,1);
TPR_632p_ste = zeros(101,1);
TPR_list = zeros(101,Nrand);

for i = 1:Nrand
    data_oob = data_orig(pts_test{i})*pol;
    outcome_oob = outcome_orig(pts_test{i});
    try
        [~,TPR_test(:,i),~,auc_test(i)] = perfcurve(outcome_oob,data_oob,2,'XVals',0:0.01:1,'UseNearest','off');
    catch
        TPR_test(:,i) = NaN;
        auc_test(i) = NaN;
    end
end

% 632 bootstrap average trn+test

[auc_632p(1),auc_632p(2),auc_list] = kyu_632plusbootstrap(auc_trn,auc_test,0.5);

for j = 1:101
    q = sum(data_orig*pol>t_trn(j));
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
se_opt = [mean(se_opt_list);2*std(se_opt_list)/sqrt(Nrand)];
sp_opt = [mean(sp_opt_list);std(sp_opt_list)/sqrt(Nrand)];
ROCY_632p = TPR_632p;
ROCY_632p_ste = TPR_632p_ste;
