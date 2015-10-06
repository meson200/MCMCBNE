function [selected,CEmin,CEvar,CE_rand,blanket] = KSfilter(data,labels,N,k,verbose)

% Variable selection algorithm by Koller & Sahami (1996)
% data: discretized data matrix. rows: variables, columns: instances (BNT format) 
%       last row of data should be a class
% N: desired number of variables to be selected
% k: size of markov blanket
% verbose: display messages (1) or not (0)
% outputs:
% selected: a binary list of length=N(variable) indicating if the 
% CEmin: CE of the eliminated variable at each round
% CEvar: CE of each variable when it was eliminated
% CErand: CE of a variable with random values i.e. the cutoff for
% usefulness
% blanket: IxJ matrix denoting frequency of a variable (i) used as a blanket for removing a variable (j) 

data_X = data(1:end-1,:);
class=  data(end,:)';
ndata = size(data_X,2);
nvar = size(data_X,1);
nvar_original = nvar;
blanket = zeros(nvar,nvar);
nelim = nvar-N;
CEmin = zeros(1,nelim);
CEvar = zeros(1,nvar);
% initiate correlation matrix (cross-entropy)
%rho = abs(corrcoef(data_X_c'));
rho = zeros(nvar,nvar);
for i = 1:nvar
    for j = 1:nvar
        rho(i,j) = AvgCrossEnt(i,j,data_X,class);
        %rho(i,j) = 
    end    
end
elim_order = [];
selected = 1:nvar;
% backward selection: remove variables one by one
for i = 1:nelim
    nvar = numel(selected);
    CE = zeros(nvar,1);
    if k>0
        % B is a matrix that stores Markov blanket variables for a variable
        % VAR at each row B(VAR,:)
        B = zeros(nvar,k); 
    else
        B = [];
    end
    for j = 1:nvar
        [~,lst] = sort(rho(selected(j),:));
        for m = 1:numel(elim_order)
            lst(lst==elim_order(m)) = []; % the eliminated variables can't form a blanket
        end
        lst(lst==selected(j)) = []; % one can't eliminate oneself
        %disp(lst);
        if k>0
            for m = 1:k
                B(j,m) = lst(m); % choose the highest correlated variable(s) as a blanket
            end
            CE(j) = AvgCrossEnt(selected(j),B(j,:),data_X,class);
        else
            CE(j) = AvgCrossEnt(selected(j),[],data_X,class);
        end    
    end
    [CEmin(i),E_ind] = min(CE);
    MinIndice = find(CE == CEmin(i)); 
    if numel(MinIndice)>1 % when multiple variables record minimum CE
        RandInt = floor(1 + numel(MinIndice)*rand(1)); % pick one of them randomly
        E_ind = MinIndice(RandInt);
    end
    E = selected(E_ind); % eliminate the variable with the lowest cross-entropy (most redundant)
    CEvar(E) = CEmin(i);
    elim_order = [elim_order E];
    selected(selected==E) = [];
    Bstr = [];
    for m = 1:k
        Bstr = strcat(Bstr,',',labels{B(E_ind,m)});
        blanket(B(E_ind,m),E) = 1;
    end
    str1 = sprintf('Removed %s (CE=%f), blanket: %s',labels{E},CEmin(i),Bstr);
    str2 = sprintf('# of variables left: %d',nvar-1);
    if verbose>0
        disp(str1);
        disp(str2);
    end
    
end

%cross-entropy of a random variable wrt the class
CE_rand_temp = zeros(100,1);
for pp = 1:100
    kk = zeros(nvar_original,1);
    randvec = round(rand(ndata,1));
    randvec = randvec + 1;
    CE_rand_temp(pp) = AvgCrossEnt(nvar_original+1,[],[data_X;randvec'],class);
    for uu = 1:nvar_original
        kk(uu) = AvgCrossEnt(nvar_original+1,uu,[data_X;randvec'],class);
    end
    CE_rand_temp(pp) = median(kk);
end
CE_rand = median(CE_rand_temp);

varstr = [];
for k = 1:numel(selected)
    varstr = strcat(varstr,',',labels{selected(k)});
end
if verbose>0
    disp(['variables selected: ',varstr]);
end


function CE = AvgCrossEnt(e,B,data_X,class)

% evaluates an average cross entropy between the two distribution
% P(C|B) and P(C|e,B) for measuring information gain with a variable e
% under the presence of a blanket B
% inputs: 
% e: variable for which the information gain is evaluated
% B: markov blanket of a class
% data_X: data matrix, excludes a class
% class: a 1D array for a class variable

nB = numel(B);
ncases = 2*(1+nB);
count_B_e = compute_counts([data_X(B,:); data_X(e,:)], 2*ones(1,numel(B)+1));
size_B_e = size(count_B_e);
count_B_e = count_B_e(:);
count_c_B = compute_counts([class'; data_X(B,:)], 2*ones(1,numel(B)+1));
size_c_B = size(count_c_B);
count_c_B = count_c_B(:);
count_c_B_e = compute_counts([class'; data_X(B,:); data_X(e,:);], 2*ones(1,numel(B)+2));
size_c_B_e = size(count_c_B_e);
count_c_B_e = count_c_B_e(:);
N = size(data_X,2); % sample size
CE = 0;
for i = 1:ncases
    Sub = cell(1,nB+1);
    [Sub{:}] = ind2sub(size_B_e,i);
    Sub2 = cell2num(Sub);
    % compute p1 = P(C|B,e)
    temp1 = num2cell([1 Sub2]);
    temp2 = num2cell([2 Sub2]);
    ind_hi = sub2ind(size_c_B_e,temp1{:}); 
    ind_lo = sub2ind(size_c_B_e,temp2{:}); 
    p1 = [count_c_B_e(ind_hi) count_c_B_e(ind_lo)];
    if sum(p1) == 0
        p1 = [0.5 0.5];
    else
        p1 = p1/sum(p1);
    end
    % compute p2 = P(C|B)
    temp1 = num2cell([1 Sub2(1:end-1)]);
    temp2 = num2cell([2 Sub2(1:end-1)]);
    ind_hi = sub2ind(size_c_B,temp1{:}); 
    ind_lo = sub2ind(size_c_B,temp2{:}); 
    p2 = [count_c_B(ind_hi) count_c_B(ind_lo)];
    if sum(p2) == 0
        p2 = [0.5 0.5];
    else
        p2 = p2/sum(p2);
    end
    count = count_B_e(i);
    CE = CE + cross_entropy(p1,p2)*count/N;
end


