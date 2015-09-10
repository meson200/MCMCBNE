function [bnet,engine] = kyu_BN_paramlearn(data,dag,ns,alpha_d,option)

% Construct a full Bayesian Network model by learning parameters (CPTs)
% given the graph and data
% input
% data: training data with a BNT format (rows:nodes,cols:patients)
% dag: a DAG matrix 
% ns: node sizes
% alpha_d: equivalent sample size for dirichlet prior
% option:
% 2: EM algorithm (in case of missing data)
% 1: max. likelihood (filled data)

num_nodes = numel(ns);
bnet = mk_bnet(dag, ns);
for i=1:num_nodes
    bnet.CPD{i} = tabular_CPD(bnet, i,'prior_type','dirichlet','dirichlet_type','BDeu','dirichlet_weight',alpha_d);
end

cases = num2cell(data);
for i = 1:size(cases,1)
    for j = 1:size(cases,2)
        if isnan(cases{i,j}) cases{i,j} = [];
        end
    end
end
if option == 1
    bnet = learn_params_ml(bnet,data);
    engine = jtree_inf_engine(bnet);
else
     max_iter = 20;
     engine = jtree_inf_engine(bnet);
    [bnet, LLtrace, engine] = learn_params_em(engine, cases , max_iter);
end
