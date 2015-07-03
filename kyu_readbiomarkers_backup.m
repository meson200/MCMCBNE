function [data_raw,data_raw_missing] = kyu_readbiomarkers_backup(NameoftheFile,FracSize,imp)

% this function imports biomarker measurement data from .xls and
% saves it as a matrix of marker concentration values in a raw scale
% each time point is marked by the following identifier:
% 1: pre-RT, 2: intra-RT, 3: end-RT, 4: post3months, 5: post6months
% for SBRT:
% 1: pre-RT, 2: end-RT, 3: post3months, 4: post6months
%
% inputs
% NameoftheFile: path to the .xls file
% FracSize: a vector of fraction size for the patients
% imp: imputation method (1: median, 2: KNN imputation) 
%
% outputs
%
% data_raw: a struct variable with the size = (No. of variable found in the file)
% and missing data filled-in
% data_raw_missing: same as data_raw, except that missing entries left NaN


data_raw = struct('name',{},'value',{});
[num,txt] = xlsread(NameoftheFile);
time = [];
id = [];
no_marker = (size(txt,2)-1)/2;

for k=2:size(txt,1)
    [tmp1,tmp2] = strtok(txt(k,1),'-');
    time = [time -1*str2num(tmp2{1,1})];    
    %tmp3 = strtok(tmp1,'L');  
    inst = tmp1{1}(1);
    tmp3 = tmp1{1}(2:end);
    
    id = [id str2num(tmp3)];
end
%idss = unique(id,'first');
[dummy,I]=unique(id,'first');
idss=id(sort(I));
SBRTids = idss(find(FracSize~=2));
samplesize = numel(idss);
%fraction = []; % binary array to denote fractionation (2:SBRT, 1:conventional)
%id2  = [];
% for i=1:samplesize
%     if numel(find(id==idss(i)))<5 % SBRT patients have less than 5 time points
%         fraction = [fraction 2];
%         id2 = [id2 i*ones(1,4)]; 
%         time = [time 1 2 4 5]; % ad-hoc assignment for SBRT pts only: use end-RT level as if it were intra-RT, 3-month as if end-RT 
%     else % conventional 
%         fraction = [fraction 1];
%         id2 = [id2 i*ones(1,5)];
%         time = [time 1 2 3 4 5];
%     end
% end   

    
num_of_var_per_marker = 4;
%3 time points (preRT,intraRT,endRT) 
num_var = no_marker*num_of_var_per_marker;

varname = cell(1,num_var);
for i=1:no_marker
    markername = txt{1,2*i};
    ind = num_of_var_per_marker*(i-1); 
    varname{1,ind+1} = [markername '_pre'];
    varname{1,ind+2} = [markername '_intra'];
    varname{1,ind+3} = [markername '_end'];
    varname{1,ind+4} = [markername '_ratio'];
end



temp_missing = zeros(samplesize,no_marker*(num_of_var_per_marker-1));
for i=1:no_marker
    %each biomarker
    temp_pre = zeros(samplesize,1);
    temp_intra = zeros(samplesize,1);
    temp_post = zeros(samplesize,1);
    for j = 1:size(num,1)
        row = find(id(j)==idss);
        if time(j) == 1
            temp_pre(row) = num(j,2*i-1);
        elseif time(j) == 2
            temp_intra(row) = num(j,2*i-1);
        elseif time(j) == 3 && isempty(find(id(j)==SBRTids)) % conventional patients
            temp_post(row) = num(j,2*i-1); 
        end
    end
    temp_missing(:,(num_of_var_per_marker-1)*(i-1)+1) = temp_pre; 
    temp_missing(:,(num_of_var_per_marker-1)*(i-1)+2) = temp_intra;
    temp_missing(:,(num_of_var_per_marker-1)*(i-1)+3) = temp_post;
    
end

% imputation 
temp = temp_missing;
switch imp
    case 1
        temp = medianimpute(temp_missing);
    case 2
        temp = knnimpute_kyu(temp_missing');
        temp = temp';
end

% put everything into the data_raw cell array
for i=1:no_marker
    for j = 1:4
       data_raw(4*(i-1)+j).name{1,1} = varname{1,4*(i-1)+j}; 
       data_raw_missing(4*(i-1)+j).name{1,1} = varname{1,4*(i-1)+j}; 
       if j == 4
           data_raw_missing(4*i).value = (temp_missing(:,3*i-1)-temp_missing(:,3*i-2))./temp_missing(:,3*i-2)*100;
           data_raw(4*i).value = (temp(:,3*i-1)-temp(:,3*i-2))./temp(:,3*i-2)*100;
       else
           data_raw_missing(4*(i-1)+j).value = temp_missing(:,3*(i-1)+j);
           data_raw(4*(i-1)+j).value = temp(:,3*(i-1)+j);
       end
    end
end