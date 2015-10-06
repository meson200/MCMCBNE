function [acc,lst,p] = kyu_BN_classify_Bayesian(data,posterior,N,Ps,ClassRate)

% calculate marginal probabilities P(Class=1) and P(Class=2) for a single new
% instance 
% via weighted averaging of prediction from each individual model in the
% ensemble (Ps)
% choose the value for Class that maximizes the joint P
% requires the Ps vector kyu_BN_classify_Bayesian_Ps
% for acc calculation,
% Class is classified as 2 when P(Class) > ClassRate
% acc returns 1 for correct classification

h = 0;

% for i = 1:N
%     bnet = bnetlist{i};
%     engine = jtree_inf_engine(bnet);
%     num_nodes = size(data,1);
%     evidence = cell(1,num_nodes);
%     
%     for k = 1:num_nodes-1
%         if isnan(data(k))
%             evidence{k} = [];
%         else
%             evidence{k} = data(k);
%         end
%     end
%     [engine, ll] = enter_evidence(engine, evidence);
%     marg = marginal_nodes(engine,num_nodes);
%     weight = posterior(k)/sum(posterior(1:N));
%     h = h + weight*marg.T(2);
% end


for i = 1:N
    weight = posterior(i)/sum(posterior(1:N));
    h = h + weight*Ps(i);
end
num_nodes = size(data,1);
outcome = data(num_nodes);    
p = [outcome; h];
if h<ClassRate Class_predict = 1; else Class_predict = 2; end
if Class_predict==outcome acc = 1; else acc = 0; end

tp = 0; fp = 0; tn = 0; fn = 0;
if acc == 1 
    if outcome == 2 
        tp = 1; 
    elseif outcome == 1
        tn = 1; 
    end
elseif acc == 0        
    if outcome == 1 
        fp = 1;
    elseif outcome == 2
        fn = 1; 
    end
end

lst = [tp fp tn fn];
