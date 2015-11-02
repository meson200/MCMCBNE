function data = kyu_BN_GeneratePartition(Npart,val)

% generates bootstrap samples by sampling without replacement from the
% original dataset
%
% input
% Npart: number of partitions 
% val: validation method. 'CV' (cross validation) or 'BS' (bootstrap)
%
% output 
% data: consists of two bootstrap sample sets
% data.KM: discretized with the KM method (used for structure finding)
% data.MI: discretized with the maximum MI method (used for BN parameter
% learning)
% each bootstrap set is structured into the following fields:
% train: training set, imputed/discrete
% train_c: training set, imputed/continuous 
% train_missing: train set, not imputed/discrete 
% test_missing: test set, not imputed/discrete 
% test_nointra: biomarker intra- or ratio- removed from test_missing
% test_c: test set, imputed/continuous
% bio_*: the training/testing set containing only biological variables
% phy: the test set containing probability estimates from physical models, set to zero 
% *orig*: the source dataset 

SBRTfilter = 3; % fractionation filter
[data_trn_KM,data_trn_c_KM,data_trn_missing_KM,data_trn_c_missing_KM,labels,~,~,~,~] = kyu_BN_readdata(SBRTfilter,3);
[data_trn_MI,data_trn_c_MI,data_trn_missing_MI,data_trn_c_missing_MI,~,~,studyID,FracSize,COMSI] = kyu_BN_readdata(SBRTfilter,2);
labels = labels';
% find which variables are biomarkers
category = kyu_BN_RP_CategorizeVariables(labels);
% before changing, see the category assignment in kyu_BN_RP_CategorizeVariables.m
nbio = find(category==2 | category==3); 
num_nodes = numel(labels)+1;
trn_cases = size(data_trn_KM,2);
labels{num_nodes} = 'RP';

% physical model prediction currently set to zero
% make your own physics model here (TCP,NTCP)
% PhyModels = TCP(studyID,FracSize,...) 
PhyModels = zeros(trn_cases,4);

if strcmp(val,'BS') % bootstrap validation
    data_KM = PutData([],[],Npart,nbio,labels,data_trn_KM,data_trn_c_KM,data_trn_missing_KM,data_trn_c_missing_KM,PhyModels);
    pts_train = data_KM.patients_train; 
    pts_test = data_KM.patients_test; 
else % cross-validation
    % generate partitions
    c = cvpartition(1:trn_cases,'KFold',Npart);
    for count = 1:Npart
        sb_cases = 1:trn_cases;
        idx_end = sum(c.TestSize(1:count));
        idx_start = idx_end-c.TestSize(count)+1;
        sb_cases(idx_start:idx_end) = [];
        leftout = idx_start:idx_end;
        pts_train{count} = sb_cases; 
        pts_test{count} = leftout;
    end
    data_KM = PutData(pts_train,pts_test,Npart,nbio,labels,data_trn_KM,data_trn_c_KM,data_trn_missing_KM,data_trn_c_missing_KM,PhyModels);
end
data_MI = PutData(pts_train,pts_test,Npart,nbio,labels,data_trn_MI,data_trn_c_MI,data_trn_missing_MI,data_trn_c_missing_MI,PhyModels);
data = struct('KM',data_KM,'MI',data_MI);



function Ds = PutData(pts_train,pts_test,Npart,nbio,labels,data,data_c,data_missing,data_c_missing,PhyModels)

Ds = struct('train',0,'train_missing',0,'test_missing',0,'Labels',[]);
SampleSize = size(data,2);
Nvar_total = numel(labels);
Ds.Labels = labels;
Ds.train = cell(Npart,1);
Ds.train_c = cell(Npart,1);
Ds.train_missing = cell(Npart,1);
Ds.train_c_missing = cell(Npart,1);
Ds.test_missing = cell(Npart,1);
Ds.test_nointra = cell(Npart,1);
Ds.bio_labels = labels([nbio numel(labels)]);
Ds.bio_train = cell(Npart,1);
Ds.bio_train_missing = cell(Npart,1);
Ds.bio_test_missing = cell(Npart,1);
Ds.bio_test_nointra = cell(Npart,1);
Ds.data_orig_missing = data_missing;
Ds.data_c_orig_missing = data_c_missing;
Ds.data_orig = data;
Ds.data_c_orig = data_c;
Ds.phy = cell(Npart,1);
Ds.patients_train = [];
Ds.patients_test = [];

for i = 1:Npart
        if ~isempty(pts_train) && ~isempty(pts_test)
            sb_cases = pts_train{i}; 
            leftout = pts_test{i};
        else
            ok = false;
            while ~ok
                pot = 1:SampleSize;
                sb_cases = randsample(pot,SampleSize,true);
                leftout = setdiff(pot,sb_cases);
                % accept the random draw only when there are more than 2
                % classes in the testing set
                ok = numel(unique(data(end,leftout))) > 1; 
            end
        end
        Ds.patients_train{i} = sb_cases;
        Ds.patients_test{i} = leftout;
        Ds.train{i} = data(:,sb_cases);
        Ds.train_missing{i} = data_missing(:,sb_cases);
        Ds.train_c{i} = data_c(:,sb_cases);
        Ds.train_c_missing{i} = data_c_missing(:,sb_cases);
        Ds.test_missing{i} = data(:,leftout);
        Ds.test_nointra{i} = remove_intra(Ds.test_missing{i},labels);
        Ds.test_c{i} = data_c(:,leftout);
        Ds.test_c_missing{i} = data_c_missing(:,leftout);
        shortened = [nbio Nvar_total];
        Ds.bio_train{i} = Ds.train{i}(shortened,:);
        Ds.bio_train_missing{i} = Ds.train_missing{i}(shortened,:);
        Ds.bio_test_missing{i} = Ds.test_missing{i}(shortened,:);
        Ds.bio_test_nointra{i} = Ds.test_nointra{i}(shortened,:);
        PhyModeltemp = PhyModels(leftout,:);
        Ds.phy{i} = PhyModeltemp;
end


