function [data_raw,data_raw_missing,class,pts_to_include,studyid,FracSize] = kyu_readphysical(NameoftheFile,imp,class_name)

% read the .xls file that contains physical/clinical variables
% input
%
% NameoftheFile: path to the .xls file
% imp: imputation method
% class_name: name of the variable that should be a class
%
% output
%
% data_raw: a struct for imputed variables with length = (No. of variable found in the file)
% data_raw_missing: same as data_raw, except that missing entries left NaN
% class: a vector of an outcome (RP in this case)
% pts_to_include: patients with value=1 is included into the training data
% studyid: cell array of patient identifiers. one letter for an institution followed by a
% 3-digit number (ex:W001)

data_raw = struct('name',{},'value',{});
[num,txt] = xlsread(NameoftheFile);
dim = size(txt);
samplesize = dim(1)-1;
if isnan(num(1,1))
    num = num(2:end,:);
end    
tmparray = zeros(samplesize,1);
pts_to_include = ones(samplesize,1);
FracSize = ones(samplesize,1);

% populate data_raw
% goes through column names found in the .xls file
% instructions to enumerate data depends on the type of variables
i=1;

% offset is needed for indexing numeric columns
% if the first column is string
instit = cell(samplesize,1);
offset = 1;
for k=1:size(txt,2)
    header = char(txt(1,k));
    switch header
        case 'Instit' % institution. under construction
            instit = txt(2:end,k);
        case 'Seq' % patient ID
            ptno = num(:,k-offset);
        case 'Rx'    
        case 'DOB' 
        case 'XRTStart'
        case 'RPdate'
        case 'RP_coded'    
        case 'RP_coded2'
        case 'failure'
            tmparray = txt(2:end,k);
            tmparray2 = zeros(samplesize,1);
            data_raw(i).name = txt(1,k);
            
            for u = 1:samplesize
               if strcmp('no disease',tmparray(u))
                   tmparray2(u) = 1;
               else
                   tmparray2(u) = 2;  
               end
            end    
            data_raw(i).value = tmparray2;
            i = i+1;
        case 'RPdate'    
        case 'LastFU'     
        case 'FUmonths'
        case 'RPmonths'
        case 'Expiration'    
        % skip the "diff" variables
        case 'a2m_diff'
        case 'IL6_diff'
        case 'ACE_diff'
        case 'TGF_diff'  
        % RP    
        case 'RP'
            for u=1:samplesize
                if num(u,k-offset) < 2 
                    tmparray(u) = 1; % define the event as RP grade >2 
                else
                    tmparray(u) = 2;
                end
            end
            class = tmparray; % class in this case is a binary for RP occurrence
            data_raw(i).name = {'Pneumonitis'};
            data_raw(i).value = num(:,k-offset);
            i = i+1;
        case 'age'
           tmparray = num(:,k-offset);
           data_raw(i).name = {'age'};
           data_raw(i).value = tmparray;
           i = i+1;
        case 'include'    
            pts_to_include = num(:,k-offset);
        case 'FracSize'
            data_raw(i).name = txt(1,k);
            data_raw(i).value = num(:,k-offset);
            i = i+1;            
            for u=1:samplesize
                %if num(u,k) ~= 2 
                %    SBRT(u) = 2; 
                %end
                FracSize(u) = num(u,k-offset);
            end  
        case {'chemo','tumorstage'}
            %data_raw(i).name = txt(1,k);
            %data_raw(i).value = txt(2:end,k);
            %i = i+1;
        case 'smoking'    
            tmparray = txt(2:end,k);
            tmparray2 = zeros(samplesize,1);
            data_raw(i).name = txt(1,k);
            for u = 1:samplesize
               switch tmparray{u}
                   case 'O'
                       tmparray2(u) = 2;
                   case 'ex' % currently,ex-smokers are treated as non-smokers
                       tmparray2(u) = 1;
                   case 'no'
                       tmparray2(u) = 1;
                   otherwise
                       tmparray2(u) = NaN;
               end
            end
            data_raw(i).value = tmparray2;
            i = i+1;
        case {'V20_BED','V30_BED'}
        case {'rpdate','fibdate'}    
        otherwise    
            data_raw(i).name = txt(1,k);
            data_raw(i).value = num(:,k-offset);
            i = i+1;
    end
end

% create instance idenfifier (source institution+ID number)
studyid = cell(samplesize,1);
for i = 1:samplesize
    if ptno(i) < 10
        studyid(i) = strcat(instit(i),'00',num2str(ptno(i)));
    elseif ptno(i) < 100
        studyid(i) = strcat(instit(i),'0',num2str(ptno(i)));
    else
        studyid(i) = strcat(instit(i),num2str(ptno(i)));
    end
end


temp_missing = [];
for i=1:length(data_raw)
   temp_missing = [temp_missing data_raw(i).value];
end

% imputation 
temp = temp_missing;
switch imp
    case 1
        temp = medianimpute(temp_missing);
    case 2
        temp = knnimpute_kyu(temp_missing',1);
        temp = temp';
end

data_raw_missing = data_raw;

for i=1:length(data_raw)
   data_raw(i).value = temp(:,i);
   data_raw_missing(i).value = temp_missing(:,i); 
   % identify a class variable
   if strcmp(data_raw(i).name,class_name)
       class = data_raw(i).value;
   end
end

pts_to_include = find(pts_to_include==1);








