function data_filled = medianimpute(data_missing)
% fill in the NaN by the median of the population
data_filled = data_missing;
[SizeR,SizeC] = size(data_missing);

% identify missing vals
nanVals = isnan(data_missing);
[rows,cols] = find(nanVals);
for i = 1:numel(cols)
    data_filled(rows,cols(i)) = nanmedian(data_missing(:,cols(i)));
end


