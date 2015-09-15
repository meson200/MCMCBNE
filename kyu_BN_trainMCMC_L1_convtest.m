function MCMCresults = kyu_BN_trainMCMC_L1_convtest(data,args_graph,args_MCMC)

% a main MCMC routine with parallelization

MCMCresults = struct('Label',0,'Prior',0,'EdgePosterior',0,'TopGraph',0,'DAGhistogram',0,'L1track',0,'Beta_track',0,'Beta_hist',0,'Agree_caus',0,'Agree_ipa',0);
labels = args_graph.labels;
num_nodes = args_graph.num_nodes;
ns = args_graph.ns;
dag_benchmark = args_graph.prior;
mask_caus = args_graph.mask_caus;
mask_ipa = args_graph.mask_ipa;
MCMC_maxlength = args_MCMC.ChainLength;
NoInit = args_MCMC.NoInit;
scoring_fn = args_MCMC.ScoringFn;
NoTopGraphs = args_MCMC.NoTopGraphs;
MCMC_burnin = args_MCMC.Burnin;
thinning = args_MCMC.thinning;
gsfreq = args_MCMC.gsfreq;
maxfanin = args_MCMC.MaxParents;
InitDensity = args_MCMC.InitDensity;
alpha_d = args_MCMC.alpha_d;
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
NoWorkers=NoInit;

disp('Running MCMC...')
timepts_small = thinning:thinning:(MCMC_maxlength-MCMC_burnin);
timepts_large = gsfreq:gsfreq:MCMC_maxlength-MCMC_burnin;
MCMC_n0 = 100; % graph statistics will be evaluated for every MCMC_n0 graphs
MCMC_checkpts = MCMC_maxlength/MCMC_n0; 
dim = [num_nodes num_nodes];
% L1 distance metrics: P(density), rrate (rejection rate)
% reference: SP Brooks et al., An Approachg to Diagnosing Total Variation
% Convergence of MCMC algoriths (1997)
%P = zeros(num_of_ss,MCMC_checkpts);
Bl = zeros(MCMC_checkpts,1);
%P_d = distributed.zeros(num_of_ss,MCMC_checkpts);
accept_ratio = [];

len = (MCMC_maxlength - MCMC_burnin)/thinning;    
gs = cell(NoWorkers,len);
lik = zeros(NoWorkers,len);
post = zeros(NoWorkers,len);
acr = zeros(NoWorkers,MCMC_maxlength);
betas = zeros(NoWorkers,len);


parfor i=1:NoWorkers

    % initialization of graphs
    density = InitDensity*rand(1);
    MCMC_init_dag = kyu_generateDAG(dim,density,dag_caus);
    [sampled_graphs, bs_w, acr_w,~,~,lik_w,post_w] = learn_struct_mcmc_L1(data, ns, dag_caus,dag_benchmark, 'nsamples', MCMC_maxlength,  ...
                                                    'burnin',MCMC_burnin,'scoring_fn',scoring_fn,'init_dag',MCMC_init_dag,'init_beta',beta_i, ...
                                                    'maxfanin',maxfanin,'type',nodetypes,'dirichlet',alpha_d_nodes);
    sampled_graphs_thinned = sampled_graphs(1:thinning:end); % thinning of sampled graphs    
    bs_thinned = bs_w(1:thinning:end);
    lik_thinned = lik_w(1:thinning:end);
    post_thinned = post_w(1:thinning:end);
    betas(i,:) = bs_thinned(MCMC_burnin/thinning+1:end);
    lik(i,:) = lik_thinned(MCMC_burnin/thinning+1:end);
    post(i,:) = post_thinned(MCMC_burnin/thinning+1:end);
    gs(i,:) = sampled_graphs_thinned(MCMC_burnin/thinning+1:end);
    acr(i,:) = acr_w;

end  

% put things together
% convert composite variables to numeric

gs = gs';
betas = betas';
%gs = gs(:);
% calculate L1 statistics
%     for t = 1:L1_checkpts
%         rrate = 0;
%         for i = 1:num_of_ss 
%             for j = 1:num_of_ss
%                 rrate = rrate + (1-min(1,P(i,t)/P(j,t)));
%             end
%         end
%         Bl(t) = rrate/num_of_ss/(num_of_ss-1);
%     end

% find the maximum posterior graph
[betas_trend,betas_hist,mat_post,agree_caus,agree_ipa] = gs_stats(betas,gs,MCMC_checkpts,dag_benchmark,mask_caus,mask_ipa);
[mat_freq_out,hist_peaked,top_dag] = matrix_hist(gs,gsfreq/thinning);
likmax = quantile(lik,0.9,1);
likmin = quantile(lik,0.1,1);
lik = mean(lik,1);
post = mean(post,1);
acr = mean(acr,1);
energy = post-lik;
delta = -energy./betas_trend';

% plot the graph scores


MCMCresults.Beta_track = betas_trend;
MCMCresults.Beta_hist = betas_hist;
MCMCresults.Agree_caus = agree_caus;
MCMCresults.Agree_ipa = agree_ipa;
MCMCresults.score_upper = likmax;
MCMCresults.score_lower = likmin;
MCMCresults.score = lik;

% save the results into the struct

MCMCresults.Label = labels;
MCMCresults.Prior = dag_benchmark;
MCMCresults.EdgePosterior = mat_post;
MCMCresults.TopGraph = top_dag;
MCMCresults.DAGhistogram = mat_freq_out;
MCMCresults.histpeaked = hist_peaked;    
MCMCresults.ACR = acr;
MCMCresults.L1track = Bl;
MCMCresults.SampledT_small = timepts_small;
MCMCresults.SampledT_large = timepts_large;
    

