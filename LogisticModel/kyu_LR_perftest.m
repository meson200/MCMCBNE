function [perf,ROCXY,th,p_RP,p_list] = kyu_LR_perftest(data,num_nodes,b,patientID,RPratio)

wrongpts = [];
acc = 0;
lst = [0 0 0 0];
p_RP = [];
num_val = size(data,2);
p_list = zeros(num_val,3);
for j = 1:num_val
    [acc_inc,lst_inc,pRP] = kyu_LR_classify(data(:,j),b,RPratio); 
    if acc_inc ==0 wrongpts = [wrongpts patientID(j)]; end
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
try
    [ROC_X,ROC_Y,th,auc] = perfcurve(p_RP(1,:),p_RP(2,:),2,'XVals',0:0.01:1,'UseNearest','off');
    ROCXY = [ROC_X ROC_Y];
catch
    ROCXY = [];
    th = NaN;
    auc = NaN;
end

mcc_no = lst(1)*lst(3)-lst(2)*lst(4);
mcc_de = (lst(1)+lst(2))*(lst(1)+lst(4))*(lst(3)+lst(2))*(lst(3)+lst(4));
mcc_de = sqrt(mcc_de);
mcc = mcc_no/mcc_de;


xdata = data(1:num_nodes,:)';
y = data(num_nodes,:)';
x = [xdata ones(size(xdata,1),1)]; % constant factor
eta = x * b;
mu = drxlr_invlogit(eta);

seps = sqrt(eps);
score=sum(y.*log(mu+seps)+(1-y).*log(1-mu+seps)); % log likelihood of the LR model


perf = [score,0,acc,ppv,npv,auc,mcc];