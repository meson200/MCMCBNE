function [data,data_c,data_missing,data_c_missing,labels,mi,studyid,FracSize,COMSI] = kyu_BN_readdata_combined(SBRTfilter,disc)

% reads patient data from .xls files into a data matrix 
% data matrices are in the following dimensions:
% # rows: number of the variables selected at the input dialog + 1 (outcome)
% # cols: number of patients
% the last row is reserved for an outcome vector
%
% inputs
%
% SBRTfilter = 2: include only conventioanl fx (fraction size < 3)
% SBRTfilter = 3: include only SBRT patients (fraction size > 3)
% SBRTfilter = 1: SBRT+Conventional
% 
% disc = 2: supervised (maximum mutual information)
% disc = 3: unsupervised (k-means)
%
% outputs
%
% data: discretized data/ missing data filled in
% data_c: data in a z-scored continuous scale/ missing data filled in 
% data_missing: discretized data/ missing data left as NaN
% data_c_missing: continuous cata / missing data left as NaN
% labels: a cell array of variable names
% mi: mutual information between variables and RP
% studyid: ID of the imported patients
% FracSize: fraction size of the imported patients
% COMSI: superior-inferior location of a PTV centroid 
%
% normalization option
% norm = 0: don't normalize, data in raw scale
% norm = 1: standardization (z-score) 
% norm = 2: rescaling (min-max) 
% norm = 3: standardization + hyper tangent mapping 
norm = 1;
% number of bins for discretization 
bins = 2;
% imputation method
imp = 2;
% imp = 1: fill the missing values with the group median 
% imp = 2: k-nearest neighbor (fill with the values of the most similar
% entries)

% scaling factors between the 2 institutions for some variables
% due to differences in units(eg. Gy/cGy)/ calibration factors (biomkrs)
% McGill values will be multiplied by scaling_values
% calibration factor for a2M obtained by taking the ratio
% WashU/McGill(Conventional) from June 2014 data
scaling_markers = {'age','a2m','ACE','OPN','Tx','MLD','MHD','V20','V30'};
scaling_values = [365,1/3700,1000,1000,100,100,100,0.01,0.01];

% path to the xls files
% physical and biological data separated in two .xls files
dir_name_WashU = '/Users/kyu/Desktop/Patients/NSCLC_WUSTL/';
filename_bio_WashU = 'WashU_biomarkers.xls';
filename_phy_WashU = 'WashU_features_20150301.xls';
fullpath_bio_WashU = cat(2,dir_name_WashU,filename_bio_WashU);
fullpath_phy_WashU = cat(2,dir_name_WashU,filename_phy_WashU);
dir_name_McGill = '/Users/kyu/Desktop/Patients/NSCLC_McGill/';
filename_bio_McGill = 'McGill_biomarkers_20150213.xls';
filename_phy_McGill = 'McGill_features_20150225.xls';
fullpath_bio_McGill = cat(2,dir_name_McGill,filename_bio_McGill);
fullpath_phy_McGill = cat(2,dir_name_McGill,filename_phy_McGill);

% read the WashU data first and obtain the list of variables from the user prompt 
[data_raw_phy_1,data_raw_phy_missing_1,class_1,pts_to_include_1,studyid,FracSize_1] = kyu_readphysical(fullpath_phy_WashU,imp);
studyid_1 = studyid(pts_to_include_1);
[data_raw_bio_1,data_raw_bio_missing_1] = kyu_readbiomarkers(fullpath_bio_WashU,FracSize_1,imp);
data_raw = [data_raw_bio_1 data_raw_phy_1];
data_raw_missing = [data_raw_bio_missing_1 data_raw_phy_missing_1];
for i = 1:numel(data_raw_phy_1)
   if strcmp(data_raw_phy_1(1,i).name,'PTVCOMSI') 
       COMSI_1 = data_raw_phy_1(1,i).value(pts_to_include_1);
   end
end
FracSize_1 = FracSize_1(pts_to_include_1);

% input dialog
disp('The following variables are found: \n');
for k = 1:size(data_raw,2)
    disp(sprintf('%d. %s \n',k,data_raw(k).name{1,1}));
end
resp = input('choose the variables, separated by comma: ','s');
if ~isempty(resp)
    temp = textscan(resp,'%d','delimiter',',');
    vars = cell2mat(temp);
else
    vars = 1:size(data_raw,2);  
end
labels = {};
num_nodes = numel(vars);
for i=1:num_nodes
    ind = vars(i);
    labels{i} = data_raw(ind).name{1,1};
end 
data_X_c_1 = kyu_preprocess(data_raw,labels,pts_to_include_1,[],[]);
data_X_c_missing_1 = kyu_preprocess(data_raw_missing,labels,pts_to_include_1,[],[]);
[data_raw_phy_2,data_raw_phy_missing_2,class_2,pts_to_include_2,studyid,FracSize_2] = kyu_readphysical(fullpath_phy_McGill,2);
studyid_2 = studyid(pts_to_include_2)+100;
[data_raw_bio_2,data_raw_bio_missing_2] = kyu_readbiomarkers(fullpath_bio_McGill,FracSize_2,2);
for i = 1:numel(data_raw_phy_2)
   if strcmp(data_raw_phy_2(1,i).name,'PTVCOMSI') 
       COMSI_2 = data_raw_phy_2(1,i).value(pts_to_include_2);
   end
end

% merge the dataset from the two institutions
FracSize_2 = FracSize_2(pts_to_include_2);
studyid = [studyid_1; studyid_2];
COMSI = [COMSI_1; COMSI_2];
FracSize = [FracSize_1; FracSize_2];
data_raw = [data_raw_bio_2 data_raw_phy_2];
data_raw_missing = [data_raw_bio_missing_2 data_raw_phy_missing_2];
data_X_c_2 = kyu_preprocess(data_raw,labels,pts_to_include_2,scaling_markers,scaling_values,[]);                                                                          
data_X_c_missing_2 = kyu_preprocess(data_raw_missing,labels,pts_to_include_2,scaling_markers,scaling_values,[]);
data_X_c_raw = [data_X_c_1; data_X_c_2];
data_X_c_raw_missing = [data_X_c_missing_1; data_X_c_missing_2];
class_1 = class_1(pts_to_include_1);
class_2 = class_2(pts_to_include_2);
class = [class_1; class_2];

% apply a fraction size filter
switch SBRTfilter
    case 1
       FracFilter = find(FracSize>0);
    case 2
       FracFilter = find(FracSize<2.01);
    case 3
       FracFilter = find(FracSize>3);
end
data_X_c_raw = data_X_c_raw(FracFilter,:);
data_X_c_raw_missing = data_X_c_raw_missing(FracFilter,:);
class = class(FracFilter);
FracSize = FracSize(FracFilter);
COMSI = COMSI(FracFilter);
studyid = studyid(FracFilter);

% normalize 
[data_X_c,data_X_c_missing] = kyu_normalize(data_X_c_raw,data_X_c_raw_missing,norm);

% discretize
[data_X,data_X_missing,b_filled,b_missing,mi] = kyu_discretize(data_X_c,data_X_c_missing,class,disc,bins);

% plot discretization results
%scales_plot = [3700/10^6,1,1,1,1/1000,1,1/1000,1,1,1/100,1/100,1,100,100,1/365,1];
%mp2013_plotbins(data_X_c,data_X_c_raw,b_filled,labels,scales_plot);


% eliminate missing data (for variable selection)
data_X_forKS = data_X;
class_forKS = class;
rowtodelete = [];
for i = 1:size(data_X,1)
   if ~isempty(find(isnan(data_X_missing(i,:)))) 
      rowtodelete = [rowtodelete i]; 
   end
end
data_X_forKS(rowtodelete,:) = [];
class_forKS(rowtodelete) = [];
%mi = p;

%optional: Sparse PCA
% for j = 5:num_nodes
%     [all_cards,all_vars,all_Zs] = sparsePCA(data_X_c,j,j,10,1);
%     var_selected = find(all_Zs);  
%     varstr = [];
%     for k = var_selected'
%         varstr = strcat(varstr,',',labels{k});
%     end
%     disp(['variables selected: ',varstr,' variance:',num2str(all_vars)]);
% end

events = find(class==2);

% add a class to a data
data = reshapeforbn(data_X,class);
data_c = reshapeforbn(data_X_c,class);
data_missing = reshapeforbn(data_X_missing,class);
data_c_missing = reshapeforbn(data_X_c_missing,class);

disp('--------------odds ratio of the selected variables-----------------');
y_trn = class-1;
oddsratio = zeros(3,size(data_X,2));
for i=1:num_nodes
    x_trn = data_X_c(:,i);
    [~,b_trn,se_trn,~,~] = drxlr_logistic_regression(x_trn,y_trn,300,1e-4);
    or_trn = [b_trn(1),b_trn(1)-1.96*se_trn(1),b_trn(1)+1.96*se_trn(1)];
    or_trn = exp(or_trn(1,:));
    oddsratio(:,i) = or_trn';
    str = sprintf('%d.%s : %4.2f(%4.2f,%4.2f)',i,labels{i},or_trn(1),or_trn(2),or_trn(3));
    disp(str)      
end

% save the imported data into the .csv file
content = [data_X_c_missing'; class'];
content2 = num2cell(content);
header = labels';
header = cat(1,header,'RP');
DmatF = [header content2]; 
for i = 1:size(DmatF,1)
    for j = 1:size(DmatF,2)
        if isnan(DmatF{i,j})
           DmatF{i,j} = ''; 
        end
    end
end
cell2csv('ForMIS.csv',DmatF');

function data_out = reshapeforbn(data_in,class)
    
data_out = [data_in class];
data_out = data_out';


    








