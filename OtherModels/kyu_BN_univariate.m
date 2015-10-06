function [Perfs_632p_uni,Perfs_632pM_uni] = kyu_BN_univariate(data,SBRTfilter,useNTCP)

% this function evaluates univariate AUC for the user chosen variables 
% in order to compare with multivariate models 
%
% input variables: 
%
% data: the data struct containing bootstrap testing sets (kyu_BN_GeneratePartition.m)
%
% SBRTfilter = 2: include only conventioanl fx (dose per fraction <= 2)
% SBRTfilter = 3: include only SBRT patients (fraction size <= 5)
% SBRTfilter = 1: SBRT+Conventional
%
% useNTCP: include physical model (TCP/NTCP) outputs into univariate models
% to compare (2) 
%
% output variables:
%
% Perfs_632p_uni: a struct bundle of performance metrics (auc/optimal
% SE/SP/ROC curves)
% Perfs_632pM_uni: an array of performance values from bootstrap sets

patients_test = data.MI.patients_test;

[~,data_trn_c,~,~,labels,~,~,~,~] = kyu_BN_readdata(SBRTfilter,3);

Nvar = numel(labels);
Npat = size(data_trn_c,2);
Nvar_NTCP = 4; % number of NTCP metrics, user defined
               % see the variable PhyModels in kyu_BN_GeneratePartition.m 
Nrand = numel(data.MI.patients_test);

Perfs_632p_uni = struct('auc',zeros(2,Nvar),'se_opt',zeros(2,Nvar),'sp_opt',zeros(2,Nvar),'ROCY',zeros(101,1));
Perfs_632pM_uni = struct('auc',zeros(Nrand,Nvar),'se_opt',zeros(Nrand,Nvar),'sp_opt',zeros(Nrand,Nvar));


auc_NTCP = zeros(2,4);
se_NTCP = zeros(2,4);
sp_NTCP = zeros(2,4);
ROCY_NTCP = zeros(101,Nvar_NTCP);
ROCY_NTCP_ste = zeros(101,Nvar_NTCP);
outcome_orig = data_trn_c(end,:);
for i = 1:Nvar
    data_orig = data_trn_c(i,:);
    [auc_list,se_list,sp_list,ROCY_temp,ROCY_ste_temp] = kyu_ROCY_univariate_632p(data_orig,outcome_orig,patients_test); 

    eq_sample_size = numel(find(~isnan(auc_list)));
    Perfs_632p_uni.auc(1,i) = nanmean(auc_list);
    Perfs_632p_uni.auc(2,i) = 2*nanstd(auc_list)/sqrt(eq_sample_size);
    Perfs_632p_uni.se_opt(1,i) = nanmean(se_list);
    Perfs_632p_uni.se_opt(2,i) = 2*nanstd(se_list)/sqrt(eq_sample_size);
    Perfs_632p_uni.sp_opt(1,i) = nanmean(sp_list);
    Perfs_632p_uni.sp_opt(2,i) = 2*nanstd(sp_list)/sqrt(eq_sample_size);
    Perfs_632p_uni.ROCY(:,i) = ROCY_temp;
    Perfs_632p_uni.ROCY_ste(:,i) = ROCY_ste_temp;

    Perfs_632pM_uni.auc(:,i) = auc_list';
    Perfs_632pM_uni.se_opt(:,i) = se_list';
    Perfs_632pM_uni.sp_opt(:,i) = sp_list';
   
end

Perfs_632p_uni.labels = labels;

% when NTCP option is on
% generate NTCP parameters
if useNTCP == 2
    NTCPs = zeros(Npat,4);
    NTCPs = medianimpute(NTCPs);
    NTCPs = NTCPs';
    for i = 1:Nvar_NTCP
        data_orig = NTCPs(i,:);
        [auc_NTCP(:,i),se_NTCP(:,i),sp_NTCP(:,i),ROCY_NTCP(:,i),ROCY_NTCP_ste(:,i)] = kyu_ROCY_univariate_632p(data_orig,outcome_orig,patients_test); 
    end
end










