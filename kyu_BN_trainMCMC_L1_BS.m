function [gs_top,post,betas_hist] = kyu_BN_trainMCMC_L1_BS(data,args_graph,args_MCMC)

% a main MCMC routine with parallelization

MCMCresults = struct('Label',0,'Prior',0,'EdgePosterior',0,'TopGraph',0,'DAGhistogram',0,'L1track',0,'Beta_track',0,'Beta_hist',0,'Agree_caus',0,'Agree_ipa',0);
NoFolds = numel(data);
num_nodes = args_graph.num_nodes;
ns = args_graph.ns;
dag_benchmark = args_graph.prior;
mask_caus = args_graph.mask_caus;
MCMC_maxlength = args_MCMC.ChainLength;
NoInit = args_MCMC.NoInit;
scoring_fn = args_MCMC.scoring_fn;
NoTopGraphs = args_MCMC.NoTopGraphs;
MCMC_burnin = args_MCMC.Burnin;
thinning = args_MCMC.thinning;
gsfreq = args_MCMC.gsfreq;
maxfanin = args_MCMC.MaxParents;
InitDensity = args_MCMC.InitDensity;
beta_i = args_MCMC.InitBeta;
% here you can set ESS differently for each node
alpha_d_nodes = alpha_d*ones(num_nodes,1);
alpha_d_nodes(end) = 1;
nodetypes = cell(num_nodes,1);
nodetypes(:) = {'tabular'};
dag_caus = ones(num_nodes,num_nodes) - mask_caus;
for i = 1:num_nodes
    dag_caus(i,i) = 0;
end

NoWorkers= NoInit*NoFolds;

% MCMC parameters
disp('Running MCMC...')
MCMC_n0 = 100;
L1_checkpts = MCMC_maxlength/MCMC_n0;
dim = [num_nodes num_nodes];
% L1 distance metrics: P(density), rrate (rejection rate)
% reference: SP Brooks et al., An Approachg to Diagnosing Total Variation
% Convergence of MCMC algoriths (1997)
%P = zeros(num_of_ss,L1_checkpts);
Bl = zeros(L1_checkpts,1);
%P_d = distributed.zeros(num_of_ss,L1_checkpts);
%accept_ratio = [];

len = (MCMC_maxlength - MCMC_burnin)/thinning;    
gs = cell(NoWorkers,len);
gs_top = cell(NoFolds,NoTopGraphs);
post = zeros(NoFolds,NoTopGraphs);
betas = zeros(NoWorkers,len);
betas_hist = zeros(NoFolds,201);

parfor i=1:NoWorkers

    % initialization of graphs
    datalrn = data{mod(NoWorkers,NoFolds)+1};
    density = InitDensity*rand(1);
    MCMC_init_dag = kyu_generateDAG(dim,density,dag_caus);
    [sampled_graphs, bs_w, acr_w,~,~,lik_w,post_w] = learn_struct_mcmc_L1(datalrn, ns, dag_caus,dag_benchmark, 'nsamples', MCMC_maxlength,  ...
                                                    'burnin',MCMC_burnin,'scoring_fn',scoring_fn,'init_dag',MCMC_init_dag,'init_beta',beta_i,...
                                                    'maxfanin',maxfanin,'type',nodetypes,'dirichlet',alpha_d_nodes);
    %sampled_graphs_thinned = sampled_graphs(1:thinning:end); % thinning of sampled graphs  
    gs(i,:) = sampled_graphs(1:thinning:end);
    bs_thinned = bs_w(1:thinning:end);
    betas(i,:) = bs_thinned;

end  

parfor j = 1:NoFolds

    start = NoInit*(j-1)+1;
    finish = NoInit*j;
    graphs = gs(start:finish,:);
    graphs = graphs';
    [DAGhistogram,~,~] = matrix_hist(graphs,gsfreq/thinning); 
    histo = DAGhistogram(end);
    posterior = histo.frequency;
    ngraphs = length(posterior);
    if ngraphs < NoTopGraphs
        for p = ngraphs+1:NoTopGraphs
            histo.matrix{p} = zeros(num_nodes);
            posterior(p) = 0;
        end
    end
    gs_top(j,:) = histo.matrix(1:NoTopGraphs);
    post(j,:) = posterior(1:NoTopGraphs);
    betas_fold = betas(start:finish,:);
    betas_hist_t = hist(betas_fold(:),linspace(0,20,201));
    betas_hist(j,:) = betas_hist_t/sum(betas_hist_t);
    
end





