function [num_models, model_list, model_freq] = drxlr_apply_tree_reduction_kyu(x, freq, xmodel, numTopModels, correlation_cutoff)
%  multiple models tree reduction
%
% Originally written by Issam El Naqa, 2005
% Extracted for general use by AJH, 2005
% Modified by Kyu, 2014 to generalize to different model orders
%
%
%  Usage: [num_models, model_list, model_freq] = drxlr_apply_tree_reduction(x, model_order, freq,xmodel, variables)
%  NOTE:  model_freq is UNIMPLEMENTED 

%  This function should have the ultimate frequencies of the models that have been
%  reduced in the following variable, but this is currently unimplemented
model_freq = [];

[dummy, inds]=sort(freq); inds=inds(end:-1:1);
freqs=freq(inds);
%xmodels=xmodel(:,inds);
%xmodels = xmodel{inds};
%model_list=xmodels(:,1:numTopModels);
model_list=xmodel(1:numTopModels);
num_models=size(model_list,2);
model_freq=freqs(1:numTopModels); % just truncate the rest!
i=1;
while i < num_models % top-> bottom
    disim_ij=[1:i];  % collect dissimilar i, j models
    for j=i+1:num_models
        %test_model=model_list(:,j);
        test_model = model_list{j};
        test_model_order = numel(find(test_model));
        test_corr=zeros(1,test_model_order);
        ref_model = model_list{i};        
        ref_model_order = numel(find(ref_model));
        
        if test_model_order == ref_model_order
           test_corr=zeros(1,test_model_order);
           vars_test = find(test_model);
           vars_ref = find(ref_model);
           for m = 1:ref_model_order 
               for n = 1:test_model_order
                   xsp=spearman(x(:,vars_ref(m)),x(:,vars_test(n)));
                    if abs(xsp) > correlation_cutoff
                        test_corr(m)=1;
                        %break;
                    end
                end
           end
           if sum(test_corr)<test_model_order % not similar model
               disim_ij=[disim_ij,j];
           else
               model_freq(i)=model_freq(i)+model_freq(j); % same model has been repeated
           end
        else
           disim_ij=[disim_ij,j];
        end
    end
    model_list=model_list(disim_ij);
    model_freq=model_freq(disim_ij);
    num_models=size(model_list,2);
    i=i+1;
end

return