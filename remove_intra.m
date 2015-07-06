function data_out = remove_intra(data,labels)
% simulate the absence of intra samples 
% by artifically dropping out all the intra-values

indx_intra = [];
data_out = data;
for i = 1:numel(labels)
   if ~isempty(regexpi(labels{i},'intra')) || ~isempty(regexpi(labels{i},'ratio'))
      indx_intra = [indx_intra i];
   end
end
if ~isempty(indx_intra)
    data_out(indx_intra,:) = NaN;
end