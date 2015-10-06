function [perf,ROCXY,P_th,p_list] = kyu_BN_perftest_Bayesian(data,Ps,posterior,N,patientID,ClassRate)

% prediction by Bayesian averaging of BN ensemble
% input: 
% data: testing data in a BNT format (rows for variables, cols for
% instances)
% Ps: an array of P(class) computed by the ensemble
% posterior: an array of posteriors estimated by MCMC 
% N: the number of highest posterior graphs to average for prediction
%
% output
% perf: a list of classification accuracy metrics with the following order
% likelihood/0/accuracy/PPV/NPV/AUC/MCC
% ROCXY: 101X2 matrix for points of an ROI curve 
% col 1:SP(fixed at regular inteval 0.01), col 2:SE (classification res.)
% P_th: a list of predicted class probability for all the instances
% p_list: patient ID/positive prediction P/negative prediction P 

wrongpts = [];
acc = 0;
lst = [0 0 0 0];
p_RP = [];
num_val = size(data,2);
p_list = zeros(num_val,3);
for j = 1:num_val
   %[acc_inc,lst_inc,pRP] = kyu_BN_classify_Bayesian(data(:,j),bnetlist,posterior,N);
   [acc_inc,lst_inc,pRP] = kyu_BN_classify_Bayesian(data(:,j),posterior,N,Ps(:,j),ClassRate);
   if acc_inc ==0 wrongpts = [wrongpts j]; end
   p_RP = [p_RP pRP];
   p_list(j,:) = [patientID(j) pRP(1) pRP(2)];
   acc = acc + acc_inc;   
   lst = lst + lst_inc;
end
acc = acc/num_val;
ppv = lst(1)/(lst(1)+lst(2));
npv = lst(3)/(lst(3)+lst(4));
if lst(1:2) == [0 0] ppv = NaN; end
if lst(3:4) == [0 0] npv = NaN; end
ROCXY = nan(101,2);
if numel(unique(p_RP(1,:))) > 1
   [ROCXY(:,1),ROCXY(:,2),P_th,auc] = perfcurve(p_RP(1,:),p_RP(2,:),2,'XVals',[0:0.01:1],'UseNearest','off');
   %Xvals=0:0.01:1;
   %ROC=interp1(ROC_X,ROC_Y,Xvals);
else
   auc = NaN;
   P_th = NaN;
end

mcc_no = lst(1)*lst(3)-lst(2)*lst(4);
mcc_de = (lst(1)+lst(2))*(lst(1)+lst(4))*(lst(3)+lst(2))*(lst(3)+lst(4));
mcc_de = sqrt(mcc_de);
mcc = mcc_no/mcc_de;

%dag = bnet.dag;
%score = score_dags(data,ns,{dag},'scoring_fn', 'bayesian');
score = 999; 
perf = [score,0,acc,ppv,npv,auc,mcc];
