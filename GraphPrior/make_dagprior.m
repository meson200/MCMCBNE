function [dag_prior,mask_caus,mask_ipa,category] = make_dagprior(labels,prior_type)

% Create a prior DAG from a given list of variables
% that enforces causality or previously reported relationships
% in between the chosen variables
% prior_type = 1: causality prior (non-causal relationships are assigned with
% a zero weight, which is otherwise all equal)
% prior_type = 2: biological prior (causality prior + previously reported
% biological relationships. Direct and indirect relationships given different
% weights)

% inputs
% labels: cell array of node names, including both RP and RP covariates
% prior_type: 1 for causality, 2 for biological
% outputs
% dag_prior: a matrix with size Nnodes*Nnodes saving weights between the nodes 
% dag_prior(i,j): weight from labels{i} -> labels{j} 
% mask_caus = a binary matrix for MCMC under causality . 
% 1: forbidden by causality, 0:allowed 
% mask_ipa: a binary matrix for MCMC under biological prior. 
% 1: the connection is allowed and posterior for the prior agreement is computed 
% 0: the connection is forbidden

% first categorize the variables into 5 categories
% 1: physical 
% 2: biological-baseline + PTV volume
% 3: biological-post radiation
% 4: clincal factors 
% 5: RP

nvar = numel(labels);
category = kyu_BN_RP_CategorizeVariables(labels);

% enforced causality
% Causality relationships will be set in between a pair of categories
% forbid the following non-causal links:
% post -> pre (3->2)
% biological_post -> physical (3->1)
% physical -> biological_pre (1->2)
% any variables -> clinical (*->4)
% clinical -> physical (4->1)
% all the links emanating from RP (5->1v2v3v4)
links_forbidden = {[3,2],[1 2],[3,1],[5,1],[5,2],[5,3],[5,4],[2 5],...
                   [3 4],[2,4],[1,4],[4 1]};
% IPA prior (collected from IPA database in 08/14/2013)
direct = {{'TGF','a2m'},{'IL6','a2m'},{'a2m','a2m'},{'IL6','IL6'},{'TGF','TGF'}};
indirect = {{'TGF','OPN'},{'TGF','ACE'},{'TGF','IL6'},{'OPN','IL6'},{'ACE','ACE'},{'OPN','OPN'},{'IL6','PTVvol'},{'a2m','RP'},{'ACE','RP'},{'IL6','RP'},{'MLD','RP'},{'V30','RP'},{'V20','RP'}};
forbid_sp = {{'V30','MLD'},{'V20','MLD'},{'V30','MLD_BED'},{'V20','MLD_BED'},{'MHD','PTVCOMSI'}};


if prior_type == 2 % Werhli & Husmeier
    direct_weight = 0.9;
    indirect_weight = 0.6;
else
    direct_weight = 0.5;
    indirect_weight = 0.5;
end    

dag_prior = zeros(nvar,nvar);
mask_caus = zeros(nvar,nvar);
mask_ipa = zeros(nvar,nvar);

for i=1:nvar
    for j=1:nvar
        link_cate = [category(i) category(j)];
        forbidden_cat = 0;
        forbidden_spec = 0;
        for k = 1:numel(links_forbidden)
            if isequal(link_cate,links_forbidden{k})
               forbidden_cat = 1;
               mask_caus(i,j) = 1;
            end
        end
        for m = 1:numel(forbid_sp)
            if strcmp(labels(i),forbid_sp{m}(1)) && strcmp(labels(j),forbid_sp{m}(2))
               forbidden_spec = 1;
               mask_caus(i,j) = 1;
            end
        end
        forbidden = forbidden_cat || forbidden_spec;
            
        if forbidden==1 || i==j
           dag_prior(i,j) = 0;
        else
           [bm1,~] = strtok(labels(i),'_');
           [bm2,~] = strtok(labels(j),'_');
           same_biomarker = strcmp(bm1,bm2);
           found_d = 0;
           found_i = 0;
           for p = 1:numel(direct)
              found_d = found_d || (occurrence(labels(i),direct{p}) && occurrence(labels(j),direct{p}));  
           end
           for u = 1:numel(indirect)
              found_i = found_i || (occurrence(labels(i),indirect{u}) && occurrence(labels(j),indirect{u}));
           end
           %if same_biomarker
           %    dag_prior(i,j) = direct_weight;
           %    mask_ipa(i,j) = 1;
           if found_d
               dag_prior(i,j) = direct_weight;
               mask_ipa(i,j) = 1; 
           elseif found_i
               dag_prior(i,j) = indirect_weight;
               mask_ipa(i,j) = 1;
           elseif occurrence(labels(j),{'V20','V30'}) && strcmp(labels(i),'MLD') 
               dag_prior(i,j) = 1; 
           else
               dag_prior(i,j) = 0.5;
           end
        end
    end
end




