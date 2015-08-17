function output = kyu_BN_GraphLearning_collect(result)
% collect the results of graph learning from multiple workers
% and assemble them into a single struct

output = struct('gs_top',0,'posterior',0,'beta_hist',0);
gs = [];
post = [];
beta = [];

for i = 1:numel(result)
    res = fetchOutputs(result{i});
    gs_t = res{1};
    post_t = res{2};
    beta_t = res{3};
    gs = [gs; gs_t];
    post = [post;post_t]; 
    beta = [beta;beta_t];
end
output.gs_top = gs;
output.posterior = post;
output.beta_hist = beta;
