function [selected,model_coal,lambda_hist,err,auc] = Stability_LASSO(Nrand,metric,verbose) 

% bootstrap test on LASSO variable selection
% inputs:
% Nrand: number of bootstrap repetition
% metric: performance measure to which a lambda is tuned to. 'AUC' or 'MSE'  
% verbose: on(1) off(0)
% outputs:
% selected: a binary matrix that stores selec- tion results for every bootstrap runs 
% (1: selected, 0: not selected)
% model_coal: models in the matrix 'selected' coalesced into much smaller
% number
% lambda_hist: histogram of lambda values chosen during bootstrap repetitions
% err: mean square error of models with different lambdas at each BS runs
% auc: auc values with different lambdas at each BS runs

normalize = 3; % K-means normalization
WhichFraction = 2; % include only conventioanl fx (fraction size < 3)

[~,~,data_trn_missing,data_trn_c_missing,labels,~] = kyu_BN_readdata_combined(WhichFraction,normalize);
% eliminate missing data (for variable selection)
colstodelete = [];
for i = 1:size(data_trn_missing,2)
   if ~isempty(find(isnan(data_trn_missing(:,i)))) 
      colstodelete = [colstodelete i]; 
   end
end
% remove missing data
data_trn_missing(:,colstodelete) = [];
data_trn_c_missing(:,colstodelete) = [];
data_trn = data_trn_missing;
data_trn_c = data_trn_c_missing;


labels = labels';
num_nodes = numel(labels)+1;
trn_cases = size(data_trn,2);
selected = zeros(Nrand,num_nodes-1);
%lambda = [0.02 0.04 0.06 0.08 0.1 0.12];
lambda = [0.002 0.005 0.01 0.02 0.04 0.06 0.08 0.12 0.2 0.4];
lambda_opts = zeros(Nrand,1);
err = zeros(Nrand,length(lambda));
auc = zeros(Nrand,length(lambda));

for i=1:Nrand

    % generate BS samples
    disp(['bootstrap sample #',num2str(i)]);
    pot = 1:trn_cases;
    if WhichFraction ~=3
        % to ensure the fraction of samples between two institution is
        % preserved in every bootstrap sets
        SampleWashU = 22;
        SampleWashU_missing = numel(colstodelete<=SampleWashU);
        pot1 = 1:SampleWashU-SampleWashU_missing;
        pot2 = SampleWashU-SampleWashU_missing+1:trn_cases;
        % sampling with replacement
        sb1_cases = randsample(pot1,numel(pot1),true);
        sb2_cases = randsample(pot2,numel(pot2),true);
        sb_cases = [sb1_cases sb2_cases];
    else
        sb_cases = randsample(pot,trn_cases,true);
    end   
    leftout = setdiff(pot,sb_cases);
    train = data_trn_c(:,sb_cases);
    test = data_trn_c(:,leftout);
    outcome_trn = train(num_nodes,:)-1;
    outcome_test = test(num_nodes,:)-1;
    [w_opt,lambda_opts(i),err(i,:),auc(i,:)] = kyu_LR_CV_testset(train(1:num_nodes-1,:)',outcome_trn',...
                                               test(1:num_nodes-1,:)',outcome_test',1,lambda,'metric',metric,'verbose',verbose);
    nonzeros = find(w_opt<0.001);   
    nonzeros = nonzeros(1:end-1);
    selected(i,nonzeros) = 1;
    varstring = '';
    for p = 1:numel(nonzeros)
        if p>1
            varstring = strcat(varstring,', ');
        end
        varstring = strcat(varstring,labels{nonzeros(p)});
    end
    if verbose==1
        disp(['Selected variables: ',varstring]);
    end
end

% generate lambda histogram (selection frequency of discrete lambda values)
lambda_hist = zeros(numel(lambda),1);
for i = 1:numel(lambda)
    lambda_hist(i) = numel(find(lambda_opts==lambda(i)));
end
lambda_hist = lambda_hist/Nrand;


% model coalescing
% pre-model selection: the models should contain significant prognistic
% factors 
% rowstoremove = [];
% sig = find(sp_p<0.05);
% for i = 1:Nrand
%     if selected(i,sig) ~= 1 
%        rowstoremove = [rowstoremove i]; 
%     end
% end
% selected(rowstoremove,:) = [];
pp = findcommonrows(selected);
model_coal = [];
num_models = 10;
corr_cutoff = 0.7; %cutoff for similarity between models in terms of a spearman coefficient
[num_models_new, model_list, model_freq] = drxlr_apply_tree_reduction_kyu(data_trn_c', pp.frequency, pp.row, num_models, 0.70);
[freq_sorted,sortedindex] = sort(model_freq,'descend');
model_coal.frequency = freq_sorted;
for l = 1:numel(sortedindex)
   model_coal.row{l} = model_list{l}; 
end





