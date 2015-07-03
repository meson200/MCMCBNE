function [selected,CElist,CEvar_avg,CErand,blanket,labels] = Stability_KS(blanketsize,Nvarmin,Nrand,verbose)

% bootstrap test on KS filtering result
% in order to obtain statistics on filtered variable

% inputs:
% blanketsize: the number of variables in the Markov blanket (size)
% Nvarmin: number of variables to select (desired data dimensionality)
% Nrand: number of bootstrap replicates
% verbose: display messages (1) or not (0)
% outputs:
% selected: a binary matrix that stores selec- tion results for every bootstrap runs 
% (1: selected, 0: not selected)
% CElist: avg cross entropy of removed variables at each round of elimination
% CEvar_avg: average CE of each variable when it is removed
% CErand: CE of a random variable
% blanket: a Nvar X Nvar matrix with the (i,j)th element denoting 
% a frequency that a variable i is used as a blanket to remove a variable j
% labels: a struct of the names of candidate variables to be filtered

normalize = 3; % K-means normalization
WhichFraction = 2; % include only conventioanl fx (fraction size < 3)

[~,~,data_trn_missing,data_trn_c_missing,labels,mi] = kyu_BN_readdata(WhichFraction,normalize);

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

Nvar = size(data_trn,1)-1;
blanket = zeros(Nvar,Nvar);
labels = labels';
trn_cases = size(data_trn,2);

selected = zeros(Nrand,Nvar);
removed = ones(Nrand,Nvar);
CElist = zeros(Nrand,Nvar-Nvarmin);
CEvar = zeros(Nrand,Nvar);
CErand = zeros(Nrand,1);
CEvar_avg = zeros(1,Nvar);
CEvar_std = zeros(1,Nvar);

i = 1;
while i < Nrand + 1

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
    train_c = data_trn_c(:,sb_cases);
    train = data_trn(:,sb_cases);
    train_x_cont = train_c(1:Nvar,:);
    [loc,CElist_temp,CEvar_temp,CErand_temp,blanket_t] = KSfilter(train,labels,Nvarmin,blanketsize,verbose);
    blanket = blanket + blanket_t;
    if sum(CElist_temp) < 10
        CElist(i,:) = CElist_temp; 
        CEvar(i,:) = CEvar_temp;
        selected(i,loc) = 1;
        removed(i,loc) = 0;
        CErand(i) = CErand_temp;
        i = i+1;
    end
    
end

% find the distinct combinations of variables that was selected together
uniquerows  = unique(selected,'rows');
freq = zeros(size(uniquerows,1),1);
for i = 1:size(uniquerows,1)
   for j = 1:Nrand
       if isequal(uniquerows(i,:),selected(j,:))
            freq(i) = freq(i) + 1;
       end
   end
end
[freq_sorted,sortedindex] = sort(freq,'descend');
% this is the variable combination that was selected most
uniquerows_sorted = uniquerows(sortedindex, :);

for i = 1:Nvar
    sample = CEvar(find(removed(:,i)),i);
    CEvar_avg(i) = mean(sample);
    CEvar_std(i) = 2*std(sample)/sqrt(numel(sample));
end
