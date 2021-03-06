function [found,loc] = occurrence(label,list)

% finds the occurrence of a label in a list of names
% returns the boolean for the occurrence (found) and the location of the
% label in the list (i)

found = 0;
loc = 0;
for i = 1:numel(list)
    aaa = regexpi(label,list{i});
    if iscell(aaa)
        aaa = aaa{1};
    end
    found = found || ~isempty(aaa);
    if ~isempty(aaa)
        loc = i;
    end
end