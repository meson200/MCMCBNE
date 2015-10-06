function Ps = kyu_BN_classify_Bayesian_Ps(data,bnetlist,enginelist,N)

% Derive class probability from each graphs in the ensemble bnetlist
% and save it in a vector Ps
% the Ps will be weighted by the posterior of the graph
% before being used for the final prediction
% inputs
% data: binary data with missing values (probabilistic imputation)
% bnetlist: a list of BNT bnet objects
% N: the number of graphs in the ensemble
%
% outputs
% Ps(i,j): P(class=2) for a patient j estimated by the i-th model in
% bnetlist

num_patients = size(data,2);
num_nodes = size(data,1); 
Ps = zeros(N,num_patients);
for j = 1:num_patients
    datasmall = data(:,j);
    for i = 1:N
        bnet = bnetlist{i};
        engine = enginelist{i};
        if ~isempty(bnet)
            evidence = cell(1,num_nodes);
            for k = 1:num_nodes-1
                if isnan(datasmall(k))
                    evidence{k} = [];
                else
                    evidence{k} = datasmall(k);
                end
            end
            [engine, ~] = enter_evidence(engine, evidence);
            marg = marginal_nodes(engine,num_nodes);
            Ps(i,j)= marg.T(2);
        else
            Ps(i,j)=0; 
        end
    end
end
