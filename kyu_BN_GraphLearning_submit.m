function result = kyu_BN_GraphLearning_submit(data,ChainLength,whichprior,clustername,validation,disc,index,result)

% Bayesian Network graph learning part 1
% this is where training data and parameters are set, jobs for MCMC
% simulation are created and sent to clusters  
%
% inputs:
% data: a struct array of bootstrap training data 
% (created by kyu_BN_RP_CombinedModel_generateBSdata.m)
% ChainLength: the MCMC iterations required for the posterior to converge
% (need to do the convergence test beforehand)
% whichprior: 'biological' or 'prior'(default)
% clustername: name of the cluster profile
% validation: 'bootstrap','cv'(cross validation), or 'fitting'(for 632+)
% disc: discretization method to choose, 'KM' or 'MI'
% index: required to specify which bootstrap replicates to send to a
% cluster ex) 1:100 (should be left [] for validation='cv')
% result: a cell storing the results from a previous run. first time run,
% set it []
%
% outputs:
% j: used only when cluster='westgrid'
% handed to kyu_BN_RP_CombinedModel_BS_part2 to retrieve results from a cluster
% outputs: the trained graph samples, used only when cluster='guillimin'
% results are retrieved while the jobs are being submitted because multiple
% jobs have to be run serially 

switch disc
    case 'KM'
        data_in = data.KM; % k-means discretized data for graph learning
    case 'MI'
        data_in = data.MI; % MMI based discretized data for graph learning
end
labels = data_in.Labels';
num_nodes = numel(labels);
ns = 2*ones(num_nodes,1);
% set priors
[dag_benchmark_caus,mask_caus,mask_ipa] = make_dagprior(labels,1);
[dag_benchmark_bio,~,~] = make_dagprior(labels,2);
dag_benchmark_no = [];
switch validation
    case 'bootstrap'
        trainmat = data_in.train;
    case 'cv'
        Nrand = numel(data_in.train); 
        trainmat = data_in.train;
        index = 1:Nrand;
    case 'fitting'
        Nrand = 1;
        index = 1;
        trainmat = {};
        trainmat{1} = data.KM.data_orig; 
end

disp('-------------------Bayesian Network----------------------------');

if strcmp(whichprior,'causal')
    prior = dag_benchmark_no;
elseif strcmp(whichprior,'biological')
    prior = dag_benchmark_bio;
else
    prior = dag_benchmark_no;
end

j = [];
outputs = [];
% set simulation hyperparameters
% hyperparams on structures (number of nodes, priors,...)
args_graph = struct();
args_graph.labels = labels; % a cell array of variable names
args_graph.ns = ns; % number of distinct states a node can have
args_graph.num_nodes = num_nodes; % number of nodes in the graph
args_graph.prior = prior; % a prior DAG (see make_dagprior.m)
args_graph.priortype = whichprior;
args_graph.mask_caus = mask_caus; % a mask matrix used for MCMC/causality (see make_dagprior.m)
args_graph.mask_ipa = mask_ipa; % a DAG mask used for MCMC/biological prior (see make_dagprior.m)
% hyperparams for graph MCMC 
args_MCMC = struct();
args_MCMC.Burnin = 10000; % burn in 
args_MCMC.ChainLength = ChainLength; % MCMC chain length (includes burn-in)
args_MCMC.ScoringFn = 'bayesian'; % graph scoring function
args_MCMC.alpha_d = 2; % equivalent sample size for dirichlet prior
args_MCMC.NoInit = 25; % number of random initialization on graphs
args_MCMC.InitDensity = 0.4; % density of sparse connections at DAG initialization 
args_MCMC.InitBeta = 10; % beta(t=0). (biological prior only)
args_MCMC.NoTopGraphs = 1000; % number of highest-posterior graphs to keep from MCMC samples (to save memory)
args_MCMC.thinning = 10; % store every 1 out of k samples in a chain (to save memory, little change in posterior)
args_MCMC.MaxParents = 3; % the maximum number of parents a node can have
% cluster setting
args_cluster = struct();
args_cluster.cluster = clustername;
args_cluster.index = index; % an array of indices for bootstrap samples used for graph training
% to run on the entire BS set of size 200, index should be 1:200

index = args_cluster.index;
NoInit = args_MCMC.NoInit;
Nrand = numel(index);
% determines the size of the job to be submitted each time
% the task is split into segments of jobs which will run serially
% size of each segment is FoldsPerJob X number of initialization 
FoldsPerJob = 2; 
st = 1:FoldsPerJob:Nrand;
Njobs = numel(st);
NoWorkerstoReq = NoInit*FoldsPerJob+1;
cluster = parcluster(clustername)
if isempty(result)
	result = cell(Njobs,1);
end
for i = 1:Njobs
   start = st(i);
   finish = min(st(i)+FoldsPerJob-1,numel(index));
   data_in = trainmat(index(start:finish));
   inputs = {data_in,args_graph,args_MCMC};
   j = batch(cluster,'kyu_BN_trainMCMC_L1_BS',3,inputs,'matlabpool',NoWorkerstoReq,'CurrentDirectory', '.');
   saveindex = floor(index(start)/FoldsPerJob)+1;
   ok = wait(j,'finished',100000);
   % the job will crash if it hasn't finished within 27 hours
   if ok
       result{saveindex} = fetchOutputs(j);
   else
       disp(['set ',num2str(index(start:finish)),' crashed']);   
       result{saveindex} = [];
   end
end













