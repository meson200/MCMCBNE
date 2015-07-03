function [R,params] = kyu_BN_MCMC_convtest(data,ChainLength,whichprior,clustername)

% Run a MCMC graph search on a single training dataset 
% the purpose is to monitor the convergence of posterior distribution
% and determine the necessary parameters to be used for 
% time-consuming bootstrap training/validation
%
% input:
% data : a struct array of bootstrapped training data .
% (created by kyu_BN_RP_CombinedModel_generateBSdata.m)
% ChainLength: MCMC iterations, should be long enough to capture the
% convergence (~10^5 is acceptable for 8 nodes)
% clustername: name of the cluster profile installed in your Matlab
%
% output
% *(T): evaluated at the interations specified in SampledT_large
% *(t): evaluated at the interations specified in SampledT_small
% R: a struct that contains the following fields:
% -EdgePosterior: a matrix representing the frequency of connections between
%  nodes that appeared in the entire MCMC samples 
% -TopGraph(T): the MAP DAG (the graph that was sampled most)
% -DAGhistogram(T): Posterior distribution of DAGs
% -histpeaked(T): peakedness of the posterior
% -ACR(t): acceptance rate 
% -L1track(t): L1 measure of similarity between tracks (not implemented)
% -Beta_track(t): change in sampled hyperparameter beta (bio. prior only)
% -Beta_hist: Posterior distribution of the beta (bio. prior only)
% -Agree_caus(t): similarity between sampled graphs and a causality prior
% -Agree_ipa(t): similarity between sampled graphs and a biological prior
% -score(t): graph score (args_MCMC.ScoringFn)
% -score_upper(t): upper bound of the score (top 10% quantile)
% -score_lower(t): lower bound of the score (bottom 10% quantile)
% params: a struct of hyperparameters used for the simulation, has 2
% fields:
% -graph: params related to graphs (nodes/priors/...)
% -MCMC: params related to MCMC (burnin/initialization/...)

data_in = data.KM; % k-means discretized data for graph learning
labels = data_in.Labels';
num_nodes = numel(labels);
ns = 2*ones(num_nodes,1); % every node has a binary state
data_trn = data_in.data_orig;  
[dag_benchmark_caus,mask_caus,mask_ipa] = make_dagprior(labels,1);
[dag_benchmark_ipa,~,~] = make_dagprior(labels,2);
dag_benchmark_no = [];
if strcmp(whichprior,'causal')
    prior = dag_benchmark_no;
    benchmark = dag_benchmark_caus;
elseif strcmp(whichprior,'biological')
    prior = dag_benchmark_ipa;
    benchmark = dag_benchmark_ipa;
else
    prior = dag_benchmark_no;
    benchmark = [];
end

disp('-------------------Bayesian Network----------------------------');

% set simulation hyperparameters
% hyperparams on structures (number of nodes, priors,...)
args_graph = struct();
args_graph.labels = labels; % a cell array of variable names
args_graph.ns = ns; % number of distinct states a node can have
args_graph.num_nodes = num_nodes; % number of nodes in the graph
args_graph.dag_benchmark = benchmark; % benchmark graph
args_graph.prior = prior; % a prior DAG (see make_dagprior.m)
args_graph.priortype = whichprior;
args_graph.mask_caus = mask_caus; % a mask matrix used for MCMC/causality (see make_dagprior.m)
args_graph.mask_ipa = mask_ipa; % a DAG mask used for MCMC/biological prior (see make_dagprior.m)
% hyperparams for graph MCMC 
args_MCMC = struct();
args_MCMC.Burnin = 0; % burn in 
args_MCMC.ChainLength = ChainLength; % MCMC chain length (includes burn-in)
args_MCMC.ScoringFn = 'bayesian'; % graph scoring function
args_MCMC.alpha_d = 2; % equivalent sample size for dirichlet prior
args_MCMC.NoInit = 25; % number of random initialization on graphs
args_MCMC.InitDensity = 0.4; % degree of sparsity in DAG initialization 
args_MCMC.InitBeta = 10; % biological prior hyperparameter beta at t=0
args_MCMC.NoTopGraphs = 1000; % number of highest-posterior graphs to keep from MCMC samples (to save memory)
args_MCMC.thinning = 10; % store every 1 out of k samples in a chain. controls SampledT_small
args_MCMC.gsfreq = 10000; % frequency (in iteration) at which posterior is computed from samples. controls SampledT_large
args_MCMC.MaxParents = 3; % the maximum number of parents a node can have

params = struct('graph',args_graph,'MCMC',args_MCMC);

NoJobtoReq = args_MCMC.NoInit+1;
cluster = parcluster(clustername)
inputs = {data_trn,args_graphs,args_MCMC};
j = batch(cluster,'kyu_BN_trainMCMC_L1_guillimin_convtest',1,inputs,'matlabpool',NoJobtoReq,'CurrentDirectory', '.');
ok = wait(j,'finished',100000); % a job is forced to crash after spending 10^5 seconds in queue
if ok
    R = fetchOutputs(j);
    R = R{1,1};
else
    disp('the job crashed');   
    R = [];
end


