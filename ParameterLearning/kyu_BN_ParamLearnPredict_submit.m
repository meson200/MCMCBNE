function result = kyu_BN_ParamLearnPredict_submit(data,MCMCresult,whichinput,clustername,validation,disc,EnsembleSizes,index,result)

% Learn parameters of MCMC sampled graphs and 
% obtain class probability under a bootstrap/CV validation setting
% computation will be distributed to a cluster specified at 'whichcluster'
% a single core will be assigned to param learning/testing for a single
% data fold
% due to current the maximum number of job restriction at the cluster, 
% the task is serialized into a sequence of smaller jobs of size =
% FoldsPerJob cores

% inputs:
% data: the bootstrapped/CVed dataset (see kyu_BN_GeneratePartition.m)
% MCMCresult: a struct from MCMC (see kyu_BN_GraphLearning_collect.m)
% whichinput: 
% -'missing': inference with original dataset containing missing values
% -'nointra': make inference without biomkr 'intra' or 'ratio' (so, only the baseline data)   
% clustername: name of the cluster 
% validation: 
% -'bootstrap'
% -'cv': cross validation
% -'fitting': parameters are learned from an original training dataset
% -'external': parameters are learned in a training dataset and inference
% is made in a validation set
% disc: discretization method to choose 
% EnsembleSizes: a list of ensemble sizes to test performance 
% (should be in an increasing order ex: [1 50 100])
% max(EnsembleSizes) should be < size(gs_top,2) (number of high posterior models
% obtained)
% index: required to specify which bootstrap replicates to send to a
% cluster ex) 1:100 (should be left [] for validation='cv')
% result: a cell storing the results from a previous run. first time run,
% set it []
%
% outputs:
% result: cell array of results passed from batch runs of kyu_BN_BMA.m
% one element saves the results from one batch run
% needs to be converted by kyu_BN_ParamLearnPredict_collect.m to be used
% for the next step

switch disc
    case 'KM'
        data_in = data.KM; % k-means discretized data for graph learning
    case 'MI'
        data_in = data.MI; % MMI based discretized data for graph learning
end
data_o = data_in.data_orig;
patno_train = data_in.patients_train;
patno_test = data_in.patients_test;
trainmat_m = data_in.train_missing;
if strcmp(whichinput,'missing')
    testmat_m = data_in.test_missing;
elseif strcmp(whichinput,'nointra')
    testmat_m = data_in.test_nointra;
end
alpha_d = 2; % equivalent sample size for dirichlet prior
labels = data_in.Labels;
npat = size(data_o,2);
nvar = size(data_o,1);
gs_top = MCMCresult.gs_top;
post = MCMCresult.posterior;
ClassRate = numel(find(data_o(nvar,:)==2))/npat;
Nrand = numel(index);
switch validation
    case 'bootstrap' % boostrap
        FoldsPerJob = 10; % number of workers requested per job
                           
    case 'cv' % cross validation (valid when sample size < 100)
        FoldsPerJob = Nrand;
    case 'fitting' % evaluate fitting performance (train data=test data)
        FoldsPerJob = 1;
        trainmat_m = {};
        trainmat_m{1} = data_in.data_orig_missing;
        testmat_m = {};
        testmat_m{1} = data_in.data_orig_missing;
        if trcmp(whichinput,'nointra')
            testmat_m{1} = remove_intra(testmat_m{1},labels);
        end
        patno_train = {};
        patno_train{1} = 1:npat;
        patno_test = {};
        patno_test{1} = 1:npat;
    case 'external'    
        FoldsPerJob = 1;
        trainmat_m = {};
        trainmat_m{1} = data_in.train_missing;
        testmat_m = {};
        testmat_m{1} = data_in.test_missing;
        if strcmp(whichinput,'nointra')
            testmat_m{1} = remove_intra(testmat_m{1},labels);
        end
        npat_train = size(data_in.train_missing,2);
        patno_train = {};
        patno_train{1} = 1:npat_train;
        npat_test = size(data_in.test_missing,2);
        patno_test = {};
        patno_test{1} = 1:npat_test;    

end
st = 1:FoldsPerJob:Nrand;
Njobs = numel(st);
NoWorkerstoReq = FoldsPerJob+1;
if isempty(result)
	result = cell(Njobs,1);
end

for p = 1:Njobs
    start = st(p);
    finish = min(st(p)+FoldsPerJob-1,numel(index));
    saveindex = floor(index(start)/FoldsPerJob)+1;
    graphs = gs_top(index(start:finish),:);
    posts = post(index(start:finish),:);
    data_tr = trainmat_m(index(start:finish));
    data_test = testmat_m(index(start:finish));
    pts_tr = patno_train(index(start:finish));
    pts_test = patno_test(index(start:finish));
    inputs = {graphs,posts,data_tr,data_test,FoldsPerJob,pts_tr,pts_test,ClassRate,EnsembleSizes,alpha_d};
    if strcmp(clustername,'serial')
        [Pr,Mt] = kyu_BN_BMA_serial(inputs{:});
        result{saveindex}{1} = Pr;
        result{saveindex}{2} = Mt;
    else 
        cluster = parcluster(clustername);
        j = batch(cluster,'kyu_BN_BMA',2,inputs,'matlabpool',NoWorkerstoReq,'CurrentDirectory', '.')
        ok = wait(j,'finished',100000); % wait 27 hrs maximum
        if ok
           result{saveindex} = fetchOutputs(j);
        else
           disp(['Jobs for folds',num2str([index(start) index(finish)]),'crashed']);   
           result{saveindex} = [];
        end
    end
end 





