function [acc,lst,p] = kyu_LR_classify(data,b,RPratio)

num_nodes = size(data,1);
evidence = cell(1,num_nodes);
outcome = data(num_nodes);
for k = 1:num_nodes-1
    evidence{k} = data(k);
end
evidence2 = cell2mat(evidence);
eta = [evidence2 1] * b;
mu = drxlr_invlogit(eta);

p = [outcome; mu];
if mu<RPratio RP_predict = 1; else RP_predict = 2; end
if RP_predict==outcome acc = 1; else acc = 0; end

tp = 0; fp = 0; tn = 0; fn = 0;
if outcome == 2 
    if RP_predict == 2 tp=1; else fn = 1; end
else        
    if RP_predict == 2 fp=1; else tn = 1; end
end

lst = [tp fp tn fn];