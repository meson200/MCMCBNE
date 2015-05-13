function [data_n,data_n_missing] = kyu_normalize(data,data_missing,mode)

% normalize the data
% the data is assumed to be columns = variables
% normalization option
% mode = 1: standardization (z-score) 
% mode = 2: rescaling (min-max) 
% mode = 3: standardization + hyper tangent mapping 

switch mode
    case 0
        data_n = data;
        data_n_missing = data_missing;
    case 1 
        data_n = zscore(data);
        data_n_missing = nanzscore(data_missing);
    case 2
        data_n = mapminmax(data',0,1);
        data_n = data_n';
        data_n_missing = mapminmax(data_missing',0,1);
        data_n_missing = data_n_missing';
    case 3
        % standardize the data first
        data_n = zscore(data);
        data_n_missing = nanzscore(data_missing);
        % map the data to [-1 1] via a hyperbolic tangent function 
        data_n = tansig(data_n);
        data_n_missing = tansig(data_n_missing);
end
