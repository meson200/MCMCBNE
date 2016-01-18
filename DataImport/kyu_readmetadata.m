function units = kyu_readmetadata(filepath,labels)

fid = fopen(filepath,'r');
% scan the header file
% skip the first line
line = fgetl(fid);
txt = textscan(fid,'%s%s%s', 'delimiter',',');
units = cell(numel(labels),1);
%[num,txt] = xlsread(NameoftheFile);
names = txt{1};
units_temp = txt{2};
for i = 1:numel(labels)
    [found,loc] = occurrence(labels{i},names);
    if found
        units{i} = units_temp{loc};
    end
end