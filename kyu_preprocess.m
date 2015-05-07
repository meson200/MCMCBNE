function data_X = kyu_preprocess(data_raw,labels,pts_to_include,varargin)

% creates a data matrix from the raw data imported from spreadsheets
% inputs
% data_raw: struct array from kyu_readphysical and kyu_readbiomarkers
% labels: array of variables names chosen by a user from the input dialog


% preprocess the raw data 
% 1. slice the raw data into the user-chosen variables & patients to include
% 2. Rescale some variables if necessary (demanded by varargin)
% varargin: varlable labels to be scaled (cell array), scaling factor (numeric array)

labels_to_scale = [];
if ~isempty(varargin) 
    labels_to_scale = varargin{1,1}; 
    sfacs = varargin{1,2}; 
end
num_cases = numel(pts_to_include);
num_nodes = numel(labels);
mi = zeros(num_nodes,1);
bound = zeros(num_nodes,1); 
data_X = zeros(num_cases,num_nodes);
% choose only the variables chosen by the user
for i = 1:size(data_raw,2)
    for j = 1:num_nodes 
        if strcmp(data_raw(1,i).name,labels{j})==1 
            if isempty(labels_to_scale)
                found = 0;
            else
                found = [];
                for k = 1:numel(labels_to_scale)
                    lbl = data_raw(1,i).name{1};
                    foundit = ~isempty(strfind(lbl,labels_to_scale{k})) && isempty(findstr(data_raw(1,i).name{1},'ratio'));
                    found = [found foundit];    
                end
            end
            if sum(found)>0
                rawdata = data_raw(i).value(pts_to_include)*sfacs(find(found));
            else    
                rawdata = data_raw(i).value(pts_to_include);
            end
            data_X(:,j) = rawdata;
        end
    end
end

