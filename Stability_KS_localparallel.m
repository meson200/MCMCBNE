function [selected,CElist,CEvar_avg,CErand,blanket,labels] = Stability_KS_localparallel(blanketsize,Nvarmin,Nrand,verbose)

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

disc = 2; % K-means discretization
WhichFraction = 3; % include only conventioanl fx (fraction size < 3)
BS_institution = 2; % when creating bootstrap samples, 
%                     force the sampling to preserve the fraction of samples between source institution (1) or not (0)

[~,~,data_trn_missing,data_trn_c_missing,labels,mi,studyid,~,~] = kyu_BN_readdata_forme(WhichFraction,disc);

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
studyid(colstodelete) = [];
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

Nworker = 8;
Nrand_pw = Nrand/Nworker;
matlabpool open 8

CElisttemp = repmat({zeros(Nrand_pw,Nvar-Nvarmin)},[1 Nworker]);
CEvartemp = repmat({zeros(Nrand_pw,Nvar)},[1 Nworker]);
selectedtemp = repmat({zeros(Nrand_pw,Nvar)},[1 Nworker]);
removedtemp = repmat({ones(Nrand_pw,Nvar)},[1 Nworker]);
CErandtemp = repmat({zeros(Nrand_pw,1)},[1 Nworker]);
blankettemp = repmat({zeros(Nvar,Nvar)},[1 Nworker]);

parfor u =1:Nworker
    i = 1;
    while i < Nrand_pw + 1

         % generate BS samples
        disp(['bootstrap sample #',num2str(i)]);
        pot = 1:trn_cases;
        if BS_institution == 1
            % to ensure the fraction of sample size between institutions is
            % preserved in every bootstrap sets
            %instit_tags = 
            [inst,group] = kyu_GroupbyInstitution(studyid);
            sb_cases = [];
            for p = 1:length(inst)
                pot = group{p};
                sb_cases_temp = randsample(pot,numel(pot),true);
                sb_cases = [sb_cases sb_cases_temp];
            end
        else
            sb_cases = randsample(pot,trn_cases,true);
        end
        train_c = data_trn_c(:,sb_cases);
        train = data_trn(:,sb_cases);
        train_x_cont = train_c(1:Nvar,:);
        [loc,CElist_temp,CEvar_temp,CErand_temp,blanket_t] = KSfilter(train,labels,Nvarmin,blanketsize,verbose);
        disp(trace(blanket_t))
        blankettemp{u} = blankettemp{u} + blanket_t;
        if sum(CElist_temp) < 10
            CElisttemp{u}(i,:) = CElist_temp; 
            CEvartemp{u}(i,:) = CEvar_temp;
            selectedtemp{u}(i,loc) = 1;
            removedtemp{u}(i,loc) = 0;
            CErandtemp{u}(i) = CErand_temp;
            i = i+1;
        end

    end
end
matlabpool close

for u = 1:Nworker
    start = Nrand_pw*(u-1)+1;
    fin = Nrand_pw*u;
    CElist(start:fin,:) = CElisttemp{u};
    CEvar(start:fin,:) = CEvartemp{u};
    selected(start:fin,:) = selectedtemp{u};
    removed(start:fin,:) = removedtemp{u};
    CErand(start:fin) = CErandtemp{u};
    blanket = blanket + blankettemp{u};
    
end


% find the distinct combinations of variables that was selected together
% uniquerows  = unique(selected,'rows');
% freq = zeros(size(uniquerows,1),1);
% for i = 1:size(uniquerows,1)
%    for j = 1:Nrand
%        if isequal(uniquerows(i,:),selected(j,:))
%             freq(i) = freq(i) + 1;
%        end
%    end
% end
% [freq_sorted,sortedindex] = sort(freq,'descend');
% % this is the variable combination that was selected most
% uniquerows_sorted = uniquerows(sortedindex, :);
% 
for i = 1:Nvar
    sample = CEvar(find(removed(:,i)),i);
    CEvar_avg(i) = nanmean(sample);
    CEvar_std(i) = 2*nanstd(sample)/sqrt(numel(sample));
end
