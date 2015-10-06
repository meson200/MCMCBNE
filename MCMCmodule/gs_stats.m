function [betas_trend,betas_hist,mat_post,agree_caus,agree_ipa] = gs_stats(betas,mats,MCMC_checkpts,dag_benchmark,mask_caus,mask_ipa)

% Compute the statistics of the sampled DAGs
% (t): the stat is computed as a function of iteration  
% betas_trend(t): change of biological hyperparameter (beta) over time
% betas_hist: histogram of sampled beta
% mat_post: a matrix representing the frequency of connections between
% nodes during the entire simulation
% agree_caus(t): similarity between sampled graphs and a causality prior
% agree_ipa(t): similarity between sampled graphs and a biological prior
% each element of the outputs is an average stat of n0 graphs
% (n0 = chain length/MCMC_checkpts)
%
% inputs:
% betas: a trace of a hyperparameter beta
% mats: a trace of the sampled graphs
% MCMC_checkpts: the number of times the stats are calculated 
% MCMC_checkpts = chainlength/n0 
% dag_benchmark: a prior DAG (see make_dagprior.m)
% mask_dag: a mask matrix used for MCMC/causality (see make_dagprior.m)
% mask_ipa: a DAG mask used for MCMC/biological prior (see make_dagprior.m)

%n = length(mats);
[n1,n2] = size(mats);
s1 = n1/MCMC_checkpts;
dim = size(mats{1,1});
mat_post = zeros(dim);
top_dag = zeros(dim);
agree_caus = zeros(MCMC_checkpts,1);
agree_ipa = zeros(MCMC_checkpts,1);
prob_cum = cell(1, MCMC_checkpts);
if isempty(dag_benchmark)
    dag_benchmark = zeros(dim);
end


for i = 1:MCMC_checkpts
    matssub = mats((i-1)*s1+1:i*s1,:);
    matssub = matssub(:);
    diff_caus = 0;
    diff_ipa  = 0;
    agr_cr = 0;
    agr_br = 0;
    score = 0;
    posterior = 0;
    for j = 1:numel(matssub)
        newmat = matssub{j};
        mat_post = mat_post+newmat;
        % agreement with prior
        agr_c = newmat & ~mask_caus ;
        agr_cr = agr_cr + numel(find(agr_c))/numel(find(newmat));
        agr_b = newmat & mask_ipa;
        agr_br = agr_br + numel(find(agr_b))/numel(find(mask_ipa));
        % likelihood (P(D|G))
        % posterior (P(D|G)*P(G|beta))
        %temp1 = abs(dag_benchmark-newmat).*mask_caus;
        %temp2 = abs(dag_benchmark-newmat).*mask_ipa;
        %diff_caus = diff_caus + sum(temp1(:));
        %diff_ipa = diff_ipa + sum(temp2(:));
    end
    prob_cum{i} = mat_post/i/numel(matssub);
    %agree_caus(i) = diff_caus/numel(matssub)/numel(find(mask_caus));
    %agree_ipa(i) = diff_ipa/numel(matssub)/numel(find(mask_ipa));
    temp1 = abs(dag_benchmark-prob_cum{i}).*mask_caus;
    temp2 = abs(dag_benchmark-prob_cum{i}).*mask_ipa;
    
    agree_caus(i) = agr_cr/numel(matssub);
    agree_ipa(i) = agr_br/numel(matssub);
    %agree_caus(i) = sum(temp1(:))/numel(find(mask_caus));
    %agree_ipa(i) = sum(temp2(:))/numel(find(mask_ipa));
    
end
mat_post = mat_post/(n1*n2);

% for i = 1:dim(1)
%     for j=1:dim(2)
%         if mat_post(i,j)>conf
%             top_dag(i,j) = 1;
%         end
%     end
% end
%disp(mat_post);
% bh = kyu_biograph(top_dag,labels);
% set(bh.Nodes,'FontSize',20);
% view(bh)
betas_trend = mean(betas,2);
betas_hist = hist(betas(:),linspace(1,30,30));
betas_hist = betas_hist/sum(betas_hist);





