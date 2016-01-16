function [data,data_c,data_missing,data_c_missing,labels,mi,studyid,FracSize,COMSI,obj_pp] = kyu_BN_readdata(SBRTfilter,Instfilter,disc,varargin)

% reads patient data from 2 .xls files (one clinical/dosimetric, one biological) into a data matrix 
% data matrices are in the following dimensions:
% # rows: number of the variables selected at the input dialog + 1 (outcome)
% # cols: number of patients
% the last row is reserved for an outcome vector
%
% inputs
% SBRTfilter: filter applied to fraction size
% -SBRTfilter = 2: include only conventioanl fx (dose per fraction <= 2)
% -SBRTfilter = 3: include only SBRT patients (fraction size <= 5)
% -SBRTfilter = 1: SBRT+Conventional
% InstFilter: filter applied to source institution, cell array containing
% letters standing for the institutions to import, which is the first letter
% of the variable 'studyID' from kyu_readphysical.m ex: {'W','C'}
% disc: discretization option
% -disc = 1: median
% -disc = 2: supervised (maximum mutual information)
% -disc = 3: unsupervised (k-means)
% -disc = 4: use bin boundaries for training set, should be accompanied by 
% an optional argument obj_pp
% varargin: (for testing set only) a struct for normalization params created from training 
%
% outputs
% data: discretized data/ missing data filled in
% data_c: data in a z-scored continuous scale/ missing data filled in 
% data_missing: discretized data/ missing data left as NaN
% data_c_missing: continuous cata / missing data left as NaN
% labels: a cell array of variable names
% mi: mutual information between variables and RP
% studyid: cell array of patient identifiers. one letter for an institution followed by a
% 3-digit number (ex:W001)
% FracSize: fraction size of the imported patients
% COMSI: superior-inferior location of a PTV centroid 
% obj_pp: a collection of parameters used for pre-processing the data
% can be handed as an input param for another run of kyu_BN_readdata.m
% .imp: imputation parameter
% -imp = 1: fill the missing values with the group median 
% -imp = 2: k-nearest neighbor (fill with the values of the most similar
% -entries)
% .disc: discretization params. see kyu_discretize.m 
% .norm: normalization params. see kyu_normalize.m
% .labels: cell array of variable names imported to the training set

class_name = 'RP'; % choose your endpoint here. will be saved to an array "class"

if ~isempty(varargin) % import post-processing parameters if given as an input argument
    obj_pp = varargin{1,1}; 
    nbins = obj_pp.disc.bins;
else
    obj_pp = struct(); % otherwise, set new parameters here
    % normalization option
    % norm = 0: don't normalize, data in raw scale
    % norm = 1: standardization (z-score) 
    % norm = 2: rescaling (min-max) 
    % norm = 3: standardization + hyper tangent mapping 
    norm = 1;
    % number of bins for discretization 
    nbins = 2;
    % imputation method
    obj_pp.imp = 2;
end

% path to the xls files
% physical and biological data separated in two .xls files
dir_name = '~/Box Sync/SKyu/lungdata/';
filename_bio = 'biomarkers_20150622.xls';
filename_phy = 'dosimetry_clinical_20150709.xls';
filename_meta = 'lung_metadata.csv';
fullpath_bio = cat(2,dir_name,filename_bio);
fullpath_phy = cat(2,dir_name,filename_phy);
fullpath_meta = cat(2,dir_name,filename_meta);
% read the clinical/physical data first and obtain the list of variables from the user prompt 
[data_raw_phy,data_raw_phy_missing,class,pts_to_include,studyid,FracSize] = kyu_readphysical(fullpath_phy,obj_pp.imp,class_name);
[data_raw_bio,data_raw_bio_missing] = kyu_readbiomarkers(fullpath_bio,studyid,FracSize,obj_pp.imp);
data_raw = [data_raw_bio data_raw_phy];
data_raw_missing = [data_raw_bio_missing data_raw_phy_missing];
for i = 1:numel(data_raw_phy)
   if strcmp(data_raw_phy(1,i).name,'PTVCOMSI') 
       COMSI = data_raw_phy(1,i).value;
   elseif strcmp(data_raw_phy(1,i).name,'NumFrac')    
       NumFrac = data_raw_phy(1,i).value;
   end
end

% apply a SBRT filter here
switch SBRTfilter
    case 1
       FracFilter = find(NumFrac>0);
    case 2
       FracFilter = find(FracSize<=2);
    case 3
       FracFilter = find(NumFrac<=5);
end
% apply institution filter here
instfilt = [];
if ~iscell(Instfilter)
    Instfilter = {Instfilter};
end
for i = 1:numel(studyid)
    if occurrence(studyid{i}(1),Instfilter)
        instfilt = [instfilt i];
    end
end
final_filter = intersect(FracFilter,pts_to_include);
final_filter = intersect(final_filter,instfilt);
class = class(final_filter);
FracSize = FracSize(final_filter);
studyid = studyid(final_filter);
COMSI = COMSI(final_filter);
NumFrac = NumFrac(final_filter);



% select the variables to import
if isempty(varargin) % training set: shows user prompt
% input dialog
    disp('The following variables are found: \n');
    namelist = {};
    for k = 1:size(data_raw,2)
        disp(sprintf('%d. %s \n',k,data_raw(k).name{1,1}));
        namelist = [namelist data_raw(k).name];
    end
    resp = input('choose the variables, separated by comma: ','s');
    if ~isempty(resp)
        temp = textscan(resp,'%d','delimiter',',');
        vars = cell2mat(temp);
    else
        vars = 1:size(data_raw,2);  
    end
    labels = namelist(vars);
    obj_pp.labels = labels;
else % testing set: import the same variables as training set, which are saved in obj_pp
    labels = obj_pp.labels;
end
num_nodes = numel(labels);

%[data_X_c_raw,data_X_c_raw_missing] = kyu_preprocess(data_raw,data_raw_missing,labels,pts_to_include,[],[]);
[data_X_c_raw,data_X_c_raw_missing] = kyu_preprocess(data_raw,data_raw_missing,labels,final_filter,[],[]);
units = kyu_readmetadata(fullpath_meta,labels);


if isempty(varargin) % training set: create new disc/norm and save them to obj_pp
    [data_X_c,data_X_c_missing,obj_pp.norm] = kyu_normalize(data_X_c_raw,data_X_c_raw_missing,norm); % normalize 
    [data_X,data_X_missing,mi,obj_pp.disc] = kyu_discretize(data_X_c,data_X_c_missing,class,disc,nbins); % discretize
else % testing set: apply the norm/disc from the training set specified in obj_pp
    [data_X_c,data_X_c_missing,~] = kyu_normalize(data_X_c_raw,data_X_c_raw_missing,4,obj_pp.norm); % normalize 
    [data_X,data_X_missing,mi,~] = kyu_discretize(data_X_c,data_X_c_missing,class,4,nbins,obj_pp.disc); % discretize
end
b_filled = obj_pp.disc.boundary_filled;
% plot discretization results
kyu_plotbins(data_X_c,data_X_c_raw,b_filled,labels,units);


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

% add a class to a data
data = reshapeforbn(data_X,class);
data_c = reshapeforbn(data_X_c,class);
data_missing = reshapeforbn(data_X_missing,class);
data_c_missing = reshapeforbn(data_X_c_missing,class);

disp('--------------odds ratio of the selected variables-----------------');
y_trn = class-1;
oddsratio = zeros(3,size(data_X_c_missing,2));
for i=1:num_nodes
    x_trn = data_X_c_missing(:,i);
    % remove NaN
    keeprows = ~isnan(x_trn) & ~isnan(y_trn);
    x_trn_2 = x_trn(keeprows);
    y_trn_2 = y_trn(keeprows);
    % continuous data - logistic fit
    if numel(unique(x_trn_2))>2
        [~,b_trn,se_trn,LR_chi2,~] = drxlr_logistic_regression(x_trn_2,y_trn_2,300,1e-4);
        or_trn = [b_trn(1),b_trn(1)-1.96*se_trn(1),b_trn(1)+1.96*se_trn(1)];
        or_trn = exp(or_trn(1,:));
        oddsratio(:,i) = or_trn';
        p = 1-chi2cdf(LR_chi2,1); % p-value from a likelihood ratio test
    else % binary data - contingency table
        x_disc = data_missing(i,:);
        x_disc_2 = x_disc(keeprows);
        [tbl,~,p] = crosstab(x_disc_2,y_trn_2);
        if ~isempty(find(tbl(:)==0)) % when there is a zero count
            tbl = tbl + 1; % add a pseudo count, this will give an inaccurate OR
        end
        or = (tbl(1)*tbl(4))/(tbl(2)*tbl(3));
        se = sqrt(1/tbl(1) + 1/tbl(2) + 1/tbl(3) + 1/tbl(4));
        or_trn = [or, or-1.96*se, or+1.96*se];     
        oddsratio(:,i) = or_trn';  
    end
    str = sprintf('%d.%s : %4.2f(%4.2f,%4.2f), p=%0.2f',i,labels{i},or_trn(1),or_trn(2),or_trn(3),p);
    disp(str)      
end
event = numel(find(y_trn==1))/numel(y_trn)*100;
disp(['(event rate: ',num2str(event),'%)']);

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


    








