function [data_n,data_n_missing,obj_norm] = kyu_normalize(data,data_missing,mode,varargin)

% normalize the data
% the data is assumed to be columns = variables
% normalization option
% mode = 1: standardization (z-score) 
% mode = 2: rescaling (min-max) 
% mode = 3: standardization + hyper tangent mapping 
% mode = 4: apply a normalization method used for trainng set (parameters
% given in varargin)
% varargin: (for testing set only) a struct for normalization params created from training 
%
% outputs
% data_n/data_n_missing: data normalized from imputed/non-imputed data
% obj_norm: discretization field of the post processing object
% (obj_pp.disc)
% saves the following fields for later use:
% .mode: see above
% .param: parameters used for respective normalization mode (imputed data)
% .param_missing: same as above (data with nan)


if ~isempty(varargin) 
    obj_norm = varargin{1,1}; 
else
    obj_norm = struct();
    obj_norm.mode = mode;
end

switch mode
    case 0
        data_n = data;
        data_n_missing = data_missing;
    case 1 
        norm_args = struct();
        norm_args_missing = struct();
        [data_n,norm_args.mu,norm_args.sigma] = zscore(data);
        [data_n_missing,norm_args_missing.mu,norm_args_missing.sigma] = nanzscore(data_missing);
    case 2
        [data_n,norm_args] = mapminmax(data',0,1);
        data_n = data_n';
        [data_n_missing,norm_args_missing] = mapminmax(data_missing',0,1);
        data_n_missing = data_n_missing';
    case 3
        % standardize the data first
        norm_args = struct();
        norm_args_missing = struct();
        [data_n,norm_args.mu,norm_args.sigma] = zscore(data);
        [data_n_missing,norm_args_missing.mu,norm_args_missing.sigma] = nanzscore(data_missing);
        % map the data to [-1 1] via a hyperbolic tangent function 
        data_n = tansig(data_n);
        data_n_missing = tansig(data_n_missing);
    case 4
        mode_p = obj_norm.mode;
        data_n = apply_normalization(data,mode_p,obj_norm.param);
        data_n_missing = apply_normalization(data,mode_p,obj_norm.param_missing);
end
if isempty(varargin) % create a struct for normalization params if not created
    obj_norm.mode = mode;
    obj_norm.param = norm_args;
    obj_norm.param_missing = norm_args_missing;
end

function data_out = apply_normalization(data_in,mode,norm_args)
data_out = zeros(size(data_in));
switch mode
    case 1
        for u = 1:size(data_in,2)
            data_out(:,u) = (data_in(:,u)-norm_args.mu(u))/norm_args.sigma(u);
        end
    case 2
        data_out = mapminmax(data',norm_args);
        data_out = data_out';
    case 3
        data_out = (data_in-norm_args.mu)./norm_args.sigma;
        data_out = tansig(data_out);
end
