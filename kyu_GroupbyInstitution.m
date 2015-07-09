function [inst,group] = kyu_GroupbyInstitution(studyid)

% group the instances by a source institution
% that should be indicated by the first letter of "studyid"

p = [studyid{:}];
p = p(1:4:end);
histo = zeros(1,26);
for n=1:length(p)
    currentLetter=p(n);
    histo(currentLetter-64)=histo(currentLetter-64)+1;
end
lett = find(histo);
inst = cell(length(lett),1);
group = cell(length(lett),1);
for i = 1:length(lett)
   inst{i} = char(lett(i)+64); 
   group{i} = strfind(p,inst{i}); 
end